import enum

from dassp.core.structures import SegmentCopyNumberRecord, AdjacencyType, Haplotype, Phasing
from dassp.core.graph import IntervalAdjacencyGraph
import gurobipy as g
import more_itertools as miter

from simulation.parts import reverse_segment


class SCNTStrategy(enum.Enum):
    CLONAL = 0
    DISTINC = 1
    MmMIXTURE = 2


class DASSPModel(object):
    def __init__(self,
                 segments,  # list of segments
                 adjacencies,  # list of adjacencies (i.e., pairs of positions defining segments extremities)
                 scnr,  # list of segment copy number records (can be greater than the list of segments, but can not be smaller)
                 n,  # number of clones for the problem
                 allele_specific=True,  # whether to treat the input in an allele specific manner, or not (if not, all segment copy number records are translated into a
                 #   non allele-specific ones, and one then proceeds as if one of two alleles was lost for every segment
                 na_false_positive_rate=0.0,
                 adjacencies_groups=None):  # a list of list of adjacencies (groups), where each adjacency group must be fully present in at least one of the clones
        self.config = {
            "homo_na_only": True,
            "SCN_strategy": None,
            "telomere_max_allele_specific": True
        }
        if adjacencies_groups is None:
            adjacencies_groups = []
        self.allele_specific = allele_specific
        self.segments = segments
        self.adjacencies = adjacencies
        self.segments_cn_records = scnr
        self.clone_cnt = n
        self.adjacencies_groups = adjacencies_groups
        self.na_false_positive_rate = na_false_positive_rate
        self.internal_check()
        self.iag = IntervalAdjacencyGraph(segments=self.segments, adjacencies=self.adjacencies)
        self.positions_idxs = [data["object"].idx for n, data in self.iag.nodes(data=True)]
        self.scnr_by_sids = {scnr.segment.idx: scnr for scnr in self.segments_cn_records}
        self.clone_ids = list(range(self.clone_cnt))
        self.big_ms_by_a_ids = self.compute_big_ms()
        self.optimized = False

        self.variables = {
            # mixing proportions of the clones
            "Phi": {i: None for i in self.clone_ids},
            # variables specifying whether the segment has majority/minority copy number in a specific clone
            "B": {i: {s.idx: None for s in self.segments} for i in self.clone_ids},
            # copy numbers of reference adjacencies (phased versions) in each clone
            "Yr": {i: {a.idx: {Phasing.AA: None,
                               Phasing.AB: None,
                               Phasing.BA: None,
                               Phasing.BB: None}
                       for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.adjacencies)}
                   for i in self.clone_ids},
            # copy numbers of novel adjacencies in each clone. Single value is present as at most one phasing of a novel adjacency is possible
            "Yn": {i: {a.idx: None for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.adjacencies)} for i in self.clone_ids},
            # adjacency phasing indicators (clone specific).
            #   As adjacencies copy numbers will be bound to be positive integers, these indicators would ensure positivity/zero values
            "P": {i: {a.idx: {Phasing.AA: None,
                              Phasing.AB: None,
                              Phasing.BA: None,
                              Phasing.BB: None}
                      for a in self.adjacencies}
                  for i in self.clone_ids},
            # auxiliary variables representing a product of binary values (i.e., phasing indicators) and adjacency copy numbers
            "Prod": {
                # clone-specific adjacency copy numbers (as a product of adjacency copy number variables Yr/Yn and respective binary indicators
                "PY": {
                    i: {a.idx: {
                        Phasing.AA: None, Phasing.AB: None, Phasing.BA: None, Phasing.BB: None}
                        for a in self.adjacencies}
                    for i in self.clone_ids
                },
                # adjacency specific indicator to signal presence/absence across all clones
                #   only novel adjacencies have constraints on this one
                "Pc": {a.idx: None for a in self.adjacencies},
                # clone specific adjacency presence-absence indicator
                "CsA": {i: {a.idx: None for a in self.adjacencies} for i in self.clone_ids},
                # clone-specific adjacency_group presence/absence indicator
                "CsGr": {i:
                             {group.idx: None for group in self.adjacencies_groups}
                         for i in self.clone_ids},
                # indicator of an adjacency group being present across all clones
                "Grc": {group.idx: None for group in self.adjacencies_groups}
            },
            # auxiliary variables representing t_i(v) function values on each segment extremity in a clone specific manner).
            "T": {i: {p_idx: {Haplotype.A: None, Haplotype.B: None} for p_idx in self.positions_idxs} for i in self.clone_ids},
            # auxiliary variable representing a maximum of t_i(v) (over il i.e., all clones) for each segment extremity
            "M": {p_idx: {Haplotype.A: None, Haplotype.B: None} for p_idx in self.positions_idxs},
            # auxiliary variables representing presence/absence of a particular phasing for a given adjacency across clones
            "Pc": {a.idx: {Phasing.AA: None,
                           Phasing.AB: None,
                           Phasing.BA: None,
                           Phasing.BB: None}
                   for a in self.adjacencies}
        }

    def internal_check(self):
        if len(self.segments) == 0:
            raise ValueError("No segments provided")
        segments_set = {s for s in self.segments}
        scnr_segments_set = {scnr.segment for scnr in self.segments_cn_records}
        if len(segments_set.symmetric_difference(scnr_segments_set)) != 0:
            raise ValueError("Segments copy number records and list of segments do not agree")
        if self.clone_cnt <= 0:
            raise ValueError("Clone count has to be positive, {n} was provided".format(n=self.clone_cnt))
        adjacencies_set = {a.idx for a in self.adjacencies}
        for adjacency_group in self.adjacencies_groups:
            for adjacency in adjacency_group.adjacencies:
                if adjacency.idx not in adjacencies_set:
                    raise ValueError("Adjacency {a} specified in one of the adjacencies groups is not present in the list of supplied adjacencies".format(a=str(adjacency)))
        if self.na_false_positive_rate < 0 or self.na_false_positive_rate > 1:
            raise ValueError("Novel Adjacencies False Positive rate has to be a positive value between 0 and 1. Supplied: {nafpr}"
                             "".format(nafpr=str(self.na_false_positive_rate)))

    def transform_copy_number_records_to_non_allele_specific(self):
        new_cnr = []
        for scnr in self.segments_cn_records:
            new_cnr.append(SegmentCopyNumberRecord.to_non_allele_specific(scnr))
        self.segments_cn_records = new_cnr

    def build_gurobi_model(self):
        self.model = g.Model("DASSP")
        self.gm = self.model
        #####
        #
        # creating variables
        #
        ####
        for i in list(self.variables["Phi"].keys()):
            self.variables["Phi"][i] = self.gm.addVar(lb=min((0.05, 1.0 / self.clone_cnt)), ub=1.0, vtype=g.GRB.CONTINUOUS, name="phi_{{{clone_id}}}".format(clone_id=str(i)))
        for i in list(self.variables["B"].keys()):
            for segment_id in list(self.variables["B"][i].keys()):
                self.variables["B"][i][segment_id] = self.gm.addVar(vtype=g.GRB.BINARY, name="b_{{{clone_id},{j}}}".format(clone_id=str(i), j=str(segment_id)))
        for i in list(self.variables["Yr"].keys()):
            for adjacency_id in list(self.variables["Yr"][i].keys()):
                self.variables["Yr"][i][adjacency_id][Phasing.AA] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER,
                                                                                   name="y_{{r,{clone_id},{aid},AA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id][Phasing.AB] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER,
                                                                                   name="y_{{r,{clone_id},{aid},AB}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id][Phasing.BA] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER,
                                                                                   name="y_{{r,{clone_id},{aid},BA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id][Phasing.BB] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER,
                                                                                   name="y_{{r,{clone_id},{aid},BB}}".format(clone_id=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for adjacency_id in list(self.variables["Yn"][i].keys()):
                self.variables["Yn"][i][adjacency_id] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER, name="y_{{n,{i},{aid}}}".format(i=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for adjacency_id in list(self.variables["P"][i].keys()):
                self.variables["P"][i][adjacency_id][Phasing.AA] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                  name="p_{{{clone_id},{aid},AA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id][Phasing.AB] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                  name="p_{{{clone_id},{aid},AB}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id][Phasing.BA] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                  name="p_{{{clone_id},{aid},BA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id][Phasing.BB] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                  name="p_{{{clone_id},{aid},BB}}".format(clone_id=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for vertex_id in list(self.variables["T"][i].keys()):
                self.variables["T"][i][vertex_id][Haplotype.A] = self.gm.addVar(vtype=g.GRB.INTEGER, name="t_{{{clone_id},{p_id},A}}".format(clone_id=str(i), p_id=str(vertex_id)))
                self.variables["T"][i][vertex_id][Haplotype.B] = self.gm.addVar(vtype=g.GRB.INTEGER, name="t_{{{clone_id},{p_id},B}}".format(clone_id=str(i), p_id=str(vertex_id)))
        for vertex_id in list(self.variables["M"].keys()):
            self.variables["M"][vertex_id][Haplotype.A] = self.gm.addVar(vtype=g.GRB.INTEGER, name="m_{{{p_id},A}}".format(p_id=str(vertex_id)))
            self.variables["M"][vertex_id][Haplotype.B] = self.gm.addVar(vtype=g.GRB.INTEGER, name="m_{{{p_id},B}}".format(p_id=str(vertex_id)))
        for adjacency_id in list(self.variables["Pc"].keys()):
            self.variables["Pc"][adjacency_id][Phasing.AA] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},AA}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id][Phasing.AB] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},AB}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id][Phasing.BA] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},BA}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id][Phasing.BB] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},BB}}".format(aid=str(adjacency_id)))
        for i in self.clone_ids:
            for a_id in list(self.variables["Prod"]["PY"][i].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["Prod"]["PY"][i][a_id][ph] = self.gm.addVar(lb=0,
                                                                               vtype=g.GRB.INTEGER,
                                                                               name="py_{{{clone_id},{a_id},{phasing}}}"
                                                                                    "".format(clone_id=str(i), phasing=ph, a_id=str(a_id)))
        for a in self.adjacencies:
            self.variables["Prod"]["Pc"][a.idx] = self.gm.addVar(vtype=g.GRB.BINARY, name="pc_{{{a_id}}}".format(a_id=str(a.idx)))
        for i in self.clone_ids:
            for adjacency_id in self.variables["Prod"]["CsA"][i]:
                self.variables["Prod"]["CsA"][i][adjacency_id] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                name="CsA_{{{clone_id},{a_id}}}"
                                                                                     "".format(clone_id=str(i), a_id=adjacency_id))
            for adj_group in self.adjacencies_groups:
                self.variables["Prod"]["CsGr"][i][adj_group.idx] = self.gm.addVar(vtype=g.GRB.BINARY,
                                                                                  name="CsGr_{{{clone_id},{ag_id}}}"
                                                                                       "".format(clone_id=str(i), ag_id=str(adj_group.idx)))
        for adj_group in self.adjacencies_groups:
            self.variables["Prod"]["Grc"][adj_group.idx] = self.gm.addVar(vtype=g.GRB.BINARY, name="Grc_{{{ag_id}}}".format(ag_id=str(adj_group.id)))
        #####
        #
        # creating constraints
        #
        ####
        # bigM constraints to linearize integer-binary product of adjacency phasing and respective adjacency copy numbers
        for i in self.clone_ids:
            for adjacency in self.adjacencies:
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    big_m = max(1, self.big_ms_by_a_ids[adjacency.idx])
                    self.gm.addConstr(big_m * self.variables["P"][i][adjacency.idx][ph],
                                      g.GRB.GREATER_EQUAL,
                                      self.variables["Prod"]["PY"][i][adjacency.idx][ph])
                    self.gm.addConstr(self.variables["Prod"]["PY"][i][adjacency.idx][ph], g.GRB.GREATER_EQUAL, 0.0, name="lb_py")
                    if adjacency.adjacency_type == AdjacencyType.NOVEL:
                        a_cn_var = self.variables["Yn"][i][adjacency.idx]
                    else:
                        a_cn_var = self.variables["Yr"][i][adjacency.idx][ph]
                    self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][i][adjacency.idx][ph])
                    self.gm.addConstr(a_cn_var - big_m * (1 - self.variables["P"][i][adjacency.idx][ph]),
                                      g.GRB.LESS_EQUAL,
                                      self.variables["Prod"]["PY"][i][adjacency.idx][ph])
                self.gm.addGenConstrOr(self.variables["Prod"]["CsA"][i][adjacency.idx],
                                       [self.variables["P"][i][adjacency.idx][Phasing.AA],
                                        self.variables["P"][i][adjacency.idx][Phasing.AB],
                                        self.variables["P"][i][adjacency.idx][Phasing.BA],
                                        self.variables["P"][i][adjacency.idx][Phasing.BB]])


        self.gm.addConstr(g.quicksum(self.variables["Phi"][i] for i in self.variables["Phi"]), g.GRB.EQUAL, 1.0, name="convex_combination")

        for c_p, c_f in zip(self.clone_ids[:-1], self.clone_ids[1:]):
            self.gm.addConstr(self.variables["Phi"][c_p] >= self.variables["Phi"][c_f], name="clonal_ordering_{cp}_{cf}".format(cp=str(c_p), cf=str(c_f)))
        for u, v, data in self.iag.segment_edges(data=True):
            segment = data["object"]
            if segment.start_position.idx != u:
                u, v = v, u
            ####
            # processing "u" position, which is always "t" segment's extremity
            ####
            novel_adjacency_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=u))
            novel_adjacencies = [data["object"] for _, __, data in novel_adjacency_edges_w_data]
            reference_adjacency_edges_w_data = list(self.iag.ref_adjacency_edges(data=True, nbunch=u))
            reference_adjacencies = [data["object"] for _, __, data in reference_adjacency_edges_w_data]
            if len(novel_adjacencies) + len(reference_adjacencies) != 0:
                for i in self.clone_ids:
                    ####
                    # processing allele Haplotype.A
                    ####
                    t_ta_lin_exp = g.LinExpr()
                    t_ta_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_a_segmentcn.cn))
                    t_ta_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_a_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.AA]), mult=-1.0)
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.BA]), mult=-1.0)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.AA] * self.variables["P"][i][ref_adjacency.idx][Phasing.AA]), mult=-1.0)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.BA] * self.variables["P"][i][ref_adjacency.idx][Phasing.BA]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if u == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AA]), mult=na_multiplier)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.AA]), mult=na_multiplier)
                        if beg:
                            t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                            # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                        else:
                            t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                            # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                    self.gm.addConstr(self.variables["T"][i][u][Haplotype.A], g.GRB.EQUAL, t_ta_lin_exp)
                    self.gm.addConstr(t_ta_lin_exp, g.GRB.GREATER_EQUAL, 0.0, name="non_negative_telomere_score_{{{clone_id},{vertex},A}}".format(clone_id=str(i), vertex=str(u)))
                    ####
                    # processing allele Haplotype.B
                    ####
                    t_tb_lin_exp = g.LinExpr()
                    t_tb_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_b_segmentcn.cn))
                    t_tb_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_b_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.AB]), mult=-1.0)
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.BB]), mult=-1.0)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.AB] * self.variables["P"][i][ref_adjacency.idx][Phasing.AB]), mult=-1.0)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.BB] * self.variables["P"][i][ref_adjacency.idx][Phasing.BB]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if u == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BB]), mult=na_multiplier)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.BB]), mult=na_multiplier)
                        if beg:
                            t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                            # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                        else:
                            t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                            # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                    self.gm.addConstr(t_tb_lin_exp, g.GRB.GREATER_EQUAL, 0, name="non_negative_telomere_score_{{{clone_id},{vertex},B}}".format(clone_id=str(i), vertex=str(u)))
                    self.gm.addConstr(self.variables["T"][i][u][Haplotype.B], g.GRB.EQUAL, t_tb_lin_exp)
                    # self.gm.addConstr(t_tb_lin_exp >= 0, )
            ####
            # processing "v" position, which is always "h" segment's extremity
            ####
            novel_adjacency_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=v))
            novel_adjacencies = [data["object"] for _, __, data in novel_adjacency_edges_w_data]
            reference_adjacency_edges_w_data = list(self.iag.ref_adjacency_edges(data=True, nbunch=v))
            reference_adjacencies = [data["object"] for _, __, data in reference_adjacency_edges_w_data]
            if len(novel_adjacencies) + len(reference_adjacencies) != 0:
                for i in self.clone_ids:
                    ####
                    # processing allele Haplotype.A
                    ###
                    t_ha_lin_exp = g.LinExpr()
                    t_ha_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_a_segmentcn.cn))
                    t_ha_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_a_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.AA]), mult=-1.0)
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.AB]), mult=-1.0)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.AA] * self.variables["P"][i][ref_adjacency.idx][Phasing.AA]), mult=-1.0)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.AB] * self.variables["P"][i][ref_adjacency.idx][Phasing.AB]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if v == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AA]), mult=na_multiplier)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.AA]), mult=na_multiplier)
                        if beg:
                            t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                            # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                        else:
                            t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                            # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                    self.gm.addConstr(t_ha_lin_exp, g.GRB.GREATER_EQUAL, 0, name="non_negative_telomere_score_{{{clone_id},{vertex},A}}".format(clone_id=str(i), vertex=str(v)))
                    self.gm.addConstr(self.variables["T"][i][v][Haplotype.A], g.GRB.EQUAL, t_ha_lin_exp)

                    ####
                    # processing allele Haplotype.B
                    ####
                    t_hb_lin_exp = g.LinExpr()
                    t_hb_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_b_segmentcn.cn))
                    t_hb_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_b_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.BA]), mult=-1.0)
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx][Phasing.BB]), mult=-1.0)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.BA] * self.variables["P"][i][ref_adjacency.idx][Phasing.BA]), mult=-1.0)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx][Phasing.BB] * self.variables["P"][i][ref_adjacency.idx][Phasing.BB]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if v == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 < novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BB]), mult=na_multiplier)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][Phasing.BB]), mult=-1.0)
                        if beg:
                            t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.BA]), mult=na_multiplier)
                            # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][Phasing.BA]))
                        else:
                            t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx][Phasing.AB]), mult=na_multiplier)
                            # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][Phasing.AB]))
                    self.gm.addConstr(t_hb_lin_exp, g.GRB.GREATER_EQUAL, 0, name="non_negative_telomere_score_{{{clone_id},{vertex},B}}".format(clone_id=str(i), vertex=str(v)))
                    self.gm.addConstr(self.variables["T"][i][v][Haplotype.B], g.GRB.EQUAL, t_hb_lin_exp)
        novel_adjacencies = [data["object"] for u, v, data in self.iag.nov_adjacency_edges(data=True)]
        for na in novel_adjacencies:
            for i in self.clone_ids:
                self.gm.addConstr(g.quicksum([self.variables["P"][i][na.idx][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]), g.GRB.LESS_EQUAL, 1.0)
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx][Phasing.AA], [self.variables["P"][i][na.idx][Phasing.AA] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx][Phasing.AB], [self.variables["P"][i][na.idx][Phasing.AB] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx][Phasing.BA], [self.variables["P"][i][na.idx][Phasing.BA] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx][Phasing.BB], [self.variables["P"][i][na.idx][Phasing.BB] for i in self.clone_ids])
            #####
            # requiring every NA to be present at least once in some of clones
            #####
            pc_ph_lin_exp = g.LinExpr(
                self.variables["Pc"][na.idx][Phasing.AA] + self.variables["Pc"][na.idx][Phasing.AB] + self.variables["Pc"][na.idx][Phasing.BA] + self.variables["Pc"][na.idx][
                    Phasing.BB])
            self.gm.addConstr(pc_ph_lin_exp,
                              g.GRB.LESS_EQUAL,
                              1)
            self.gm.addConstr(self.variables["Prod"]["Pc"][na.idx], g.GRB.EQUAL, pc_ph_lin_exp)
            # NA involving the same segments must be on the same alleles
            s1 = self.iag.get_segment_edge(node=na.position1.idx, data=True)[2]["object"]
            s2 = self.iag.get_segment_edge(node=na.position2.idx, data=True)[2]["object"]
            if s1 == s2:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][na.idx][Phasing.AB], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][na.idx][Phasing.BA], g.GRB.EQUAL, 0)
        # NA False positive rate
        fraction = 1.0 / len(novel_adjacencies)
        fp_lin_exp = g.LinExpr(g.quicksum([self.variables["Prod"]["Pc"][na.idx] * fraction for na in novel_adjacencies]))
        self.gm.addConstr(fp_lin_exp, g.GRB.GREATER_EQUAL, 1 - self.na_false_positive_rate)
        ###
        for u, v, data in self.iag.ref_adjacency_edges(data=True):
            ra = data["object"]

            segment1 = self.iag.get_segment_edge(data=True, node=u)[2]["object"]
            segment2 = self.iag.get_segment_edge(data=True, node=v)[2]["object"]
            if segment1.start_position > segment2.end_position:
                segment1, segment2 = segment2, segment1
                u, v = v, u
            s1cnr = self.scnr_by_sids[segment1.idx]
            s2cnr = self.scnr_by_sids[segment2.idx]
            if not s1cnr.allele_specific:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.BA], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.BB], g.GRB.EQUAL, 0)
            else:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.AA] + self.variables["P"][i][ra.idx][Phasing.AB], g.GRB.LESS_EQUAL, 1)
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.BA] + self.variables["P"][i][ra.idx][Phasing.BB], g.GRB.LESS_EQUAL, 1)
            if not s2cnr.allele_specific:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.AB], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.BB], g.GRB.EQUAL, 0)
            else:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.AA] + self.variables["P"][i][ra.idx][Phasing.BA], g.GRB.LESS_EQUAL, 1)
                    self.gm.addConstr(self.variables["P"][i][ra.idx][Phasing.AB] + self.variables["P"][i][ra.idx][Phasing.BB], g.GRB.LESS_EQUAL, 1)
            #####
            #
            # constraints enforcing consistent phasing of novel adjacencies across all clones
            #
            #####
            self.gm.addGenConstrOr(self.variables["Pc"][ra.idx][Phasing.AA], [self.variables["P"][i][ra.idx][Phasing.AA] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][ra.idx][Phasing.AB], [self.variables["P"][i][ra.idx][Phasing.AB] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][ra.idx][Phasing.BA], [self.variables["P"][i][ra.idx][Phasing.BA] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][ra.idx][Phasing.BB], [self.variables["P"][i][ra.idx][Phasing.BB] for i in self.clone_ids])
            self.gm.addConstr(self.variables["Pc"][ra.idx][Phasing.AA] + self.variables["Pc"][ra.idx][Phasing.AB] + self.variables["Pc"][ra.idx][Phasing.BA] + self.variables["Pc"][ra.idx][Phasing.BB],
                              g.GRB.LESS_EQUAL, 2)
        ####
        #
        # adding balancing equations on segments extremities that are reference telomeres without novel adjacencies
        #
        ###
        for p_id in self.positions_idxs:
            segment = self.iag.get_segment_edge(node=p_id, data=True)[2]["object"]
            if len(list(self.iag.adjacency_edges(nbunch=p_id))) == 0:
                for i in self.clone_ids:
                    t_va_lin_exp = g.LinExpr(self.variables["B"][i][segment.idx]*self.scnr_by_sids[segment.idx].maj_a_segmentcn.cn)
                    t_va_lin_exp.add(g.LinExpr((1-self.variables["B"][i][segment.idx])*self.scnr_by_sids[segment.idx].min_a_segmentcn.cn))
                    self.gm.addConstr(t_va_lin_exp, g.GRB.GREATER_EQUAL, 0, name="non_negative_telomere_score_{{rt,{clone_id},{vertex},A}}".format(clone_id=str(i), vertex=str(p_id)))
                    self.gm.addConstr(self.variables["T"][i][p_id][Haplotype.A], g.GRB.EQUAL, t_va_lin_exp)
                    t_vb_lin_exp = g.LinExpr(self.variables["B"][i][segment.idx]*self.scnr_by_sids[segment.idx].maj_b_segmentcn.cn)
                    t_vb_lin_exp.add(g.LinExpr((1-self.variables["B"][i][segment.idx])*self.scnr_by_sids[segment.idx].min_b_segmentcn.cn))
                    self.gm.addConstr(t_vb_lin_exp, g.GRB.GREATER_EQUAL, 0, name="non_negative_telomere_score_{{rt,{clone_id},{vertex},B}}".format(clone_id=str(i), vertex=str(p_id)))
                    self.gm.addConstr(self.variables["T"][i][p_id][Haplotype.B], g.GRB.EQUAL, t_vb_lin_exp)
        ####
        #
        # maximum over telomeres imbalance across clones
        #
        ####
        for p_id in self.positions_idxs:
            a_extremities = [self.variables["T"][i][p_id][Haplotype.A] for i in self.clone_ids]
            b_extremities = [self.variables["T"][i][p_id][Haplotype.B] for i in self.clone_ids]
            if self.config["telomere_max_allele_specific"]:
                self.gm.addGenConstrMax(self.variables["M"][p_id][Haplotype.A], a_extremities)
                self.gm.addGenConstrMax(self.variables["M"][p_id][Haplotype.B], b_extremities)
            else:
                self.gm.addGenConstrMax(self.variables["M"][p_id][Haplotype.A], a_extremities + b_extremities)
                self.gm.addConstr(self.variables["M"][p_id][Haplotype.B], g.GRB.EQUAL, self.variables["M"][p_id][Haplotype.A])

        ####
        #
        # adjacency group constrains
        #
        ###
        for i in self.clone_ids:
            for adj_group in self.adjacencies_groups:
                self.gm.addGenConstrAnd(self.variables["Prod"]["CsGr"][i][adj_group.idx],
                                        [self.variables["Prod"]["CsA"][i][a.idx] for a in adj_group.adjacencies])

        for adj_group in self.adjacencies_groups:
            self.gm.addGenConstrOr(self.variables["Prod"]["Grc"][adj_group.idx],
                                   [self.variables["Prod"]["CsGr"][i][adj_group.idx] for i in self.clone_ids])
            self.gm.addConstr(self.variables["Prod"]["Grc"]["Grc"][adj_group.idx], g.GRB.EQUAL, 1)
        ######
        #
        # ensuring that for every segment that has a distinct majority/minority copy number states both are present in at least some clones
        #
        ######
        for s in self.segments:
            scnr = self.scnr_by_sids[s.idx]
            if scnr.distinct_min_cn:
                self.gm.addConstr(g.quicksum([self.variables["B"][i][s.idx] for i in self.clone_ids]), g.GRB.LESS_EQUAL, len(self.clone_ids) - 1)
                self.gm.addConstr(g.quicksum([self.variables["B"][i][s.idx] for i in self.clone_ids]), g.GRB.GREATER_EQUAL, 1)
            else:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["B"][i][s.idx], g.GRB.EQUAL, 1)
        ######
        #
        # ensuring that first clone has all the majority copy number states, and the rest can be arranged in any fashion
        #
        ######
        if self.config["SCN_strategy"] == SCNTStrategy.CLONAL:
            maj_clone = self.clone_ids[0]
            for s in self.segments:
                self.gm.addConstr(self.variables["B"][maj_clone][s.idx], g.GRB.EQUAL, 1)
        elif self.config["SCN_strategy"] == SCNTStrategy.MmMIXTURE:
            for s in self.segments:
                self.gm.addConstr(g.quicksum([self.variables["Phi"][i] * self.variables["B"][i][s.idx] for i in self.clone_ids]), g.GRB.GREATER_EQUAL, 0.5)
        scalar_coefficient = 1 if self.config["telomere_max_allele_specific"] else .5
        objective_lin_exp = g.LinExpr()
        objective_lin_exp.add(g.quicksum([self.variables["M"][p_id][Haplotype.A] for p_id in self.positions_idxs]), mult=scalar_coefficient)
        objective_lin_exp.add(g.quicksum([self.variables["M"][p_id][Haplotype.B] for p_id in self.positions_idxs]), mult=scalar_coefficient)
        self.gm.setObjective(objective_lin_exp, g.GRB.MINIMIZE)

    def solve_gurobi_model(self):
        self.gm.optimize()
        self.optimized = True

    def compute_big_ms(self):
        result = {}
        for adjacency in self.adjacencies:
            s1 = self.iag.get_segment_edge(node=adjacency.position1.idx, data=True)[2]["object"]
            s2 = self.iag.get_segment_edge(node=adjacency.position2.idx, data=True)[2]["object"]
            s1_all_segment_copy_number_values = [self.scnr_by_sids[s1.idx].maj_a_segmentcn.cn, self.scnr_by_sids[s1.idx].maj_b_segmentcn.cn,
                                                 self.scnr_by_sids[s1.idx].min_a_segmentcn.cn, self.scnr_by_sids[s1.idx].min_b_segmentcn.cn]
            s2_all_segment_copy_number_values = [self.scnr_by_sids[s2.idx].maj_a_segmentcn.cn, self.scnr_by_sids[s2.idx].maj_b_segmentcn.cn,
                                                 self.scnr_by_sids[s2.idx].min_a_segmentcn.cn, self.scnr_by_sids[s2.idx].min_b_segmentcn.cn]
            big_m = max(s1_all_segment_copy_number_values + s2_all_segment_copy_number_values)
            result[adjacency.idx] = big_m
        return result

    def get_adjacency_copy_number(self, adjacency, clone_id):
        assert self.optimized
        result = {}
        for phasing in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
            result[phasing] = self.variables["Prod"]["PY"][clone_id][adjacency.idx][phasing].x
        return result

    def get_segment_copy_number(self, segment, clone_id):
        assert self.optimized
        if segment.is_reversed:
            segment_idx = reverse_segment(segment=segment).idx
            reverse_segment(segment=segment)
        else:
            segment_idx = segment.idx
        result = {}
        if self.variables["B"][clone_id][segment_idx].x == 1:
            result[Haplotype.A] = self.scnr_by_sids[segment_idx].maj_a_segmentcn.cn
            result[Haplotype.B] = self.scnr_by_sids[segment_idx].maj_b_segmentcn.cn
        else:
            result[Haplotype.A] = self.scnr_by_sids[segment_idx].min_a_segmentcn.cn
            result[Haplotype.B] = self.scnr_by_sids[segment_idx].min_b_segmentcn.cn
        return result
