import enum

from dassp.core.structures import SegmentCopyNumberRecord, AdjacencyType
from dassp.core.graph import IntervalAdjacencyGraph
import gurobipy as g
import more_itertools as miter


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
                 adjacencies_groups=None):  # a list of list of adjacencies (groups), where each adjacency group must be fully present in at least one of the clones
        self.config = {
            "homo_na_only": True,
            "SCN_strategy": SCNTStrategy.CLONAL,
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
        self.internal_check()
        self.iag = IntervalAdjacencyGraph(segments=self.segments, adjacencies=self.adjacencies)
        self.positions_idxs = [data["object"].idx for n, data in self.iag.nodes(data=True)]
        self.scnr_by_sids = {scnr.segment.idx: scnr for scnr in self.segments_cn_records}
        self.clone_ids = list(range(self.clone_cnt))
        self.big_ms_by_a_ids = self.compute_big_ms()
        self.variables = {
            # mixing proportions of the clones
            "Phi": {i: None for i in self.clone_ids},
            # variables specifying whether the segment has majority/minority copy number in a specific clone
            "B": {i: {s.idx: None for s in self.segments} for i in self.clone_ids},
            # copy numbers of reference adjacencies (phased versions) in each clone
            "Yr": {i: {a.idx: {"AA": None,
                               "AB": None,
                               "BA": None,
                               "BB": None}
                       for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.adjacencies)}
                   for i in self.clone_ids},
            # copy numbers of novel adjacencies in each clone. Single value is present as at most one phasing of a novel adjacency is possible
            "Yn": {i: {a.idx: None for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.adjacencies)} for i in self.clone_ids},
            # adjacency phasing indicators (clone specific).
            #   As adjacencies copy numbers will be bound to be positive integers, these indicators would ensure positivity/zero values
            "P": {i: {a.idx: {"AA": None,
                              "AB": None,
                              "BA": None,
                              "BB": None}
                      for a in self.adjacencies}
                  for i in self.clone_ids},
            # auxiliary variables representing a product of binary values (i.e., phasing indicators) and adjacency copy numbers
            "Prod": {
                "PY": {
                    i: {a.idx: {
                        "AA": None, "AB": None, "BA": None, "BB": None}
                        for a in self.adjacencies}
                    for i in self.clone_ids
                }
            },
            # auxiliary variables representing t_i(v) function values on each segment extremity in a clone specific manner).
            "T": {i: {p_idx: {"A": None, "B": None} for p_idx in self.positions_idxs} for i in self.clone_ids},
            # auxiliary variable representing a maximum of t_i(v) (over il i.e., all clones) for each segment extremity
            "M": {p_idx: {"A": None, "B": None} for p_idx in self.positions_idxs},
            # auxiliary variables representing presence/absence of a particular phasing for a given adjacency across clones
            "Pc": {a.idx: {"AA": None,
                           "AB": None,
                           "BA": None,
                           "BB": None}
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
            for adjacency in adjacency_group:
                if adjacency.idx not in adjacencies_set:
                    raise ValueError("Adjacency {a} specified in one of the adjacencies groups is not present in the list of supplied adjacencies".format(a=str(adjacency)))

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
                self.variables["Yr"][i][adjacency_id]["AA"] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER,
                                                                             name="y_{{{clone_id},{aid},AA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id]["AB"] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER,
                                                                             name="y_{{{clone_id},{aid},AB}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id]["BA"] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER,
                                                                             name="y_{{{clone_id},{aid},BA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["Yr"][i][adjacency_id]["BB"] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER,
                                                                             name="y_{{{clone_id},{aid},BB}}".format(clone_id=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for adjacency_id in list(self.variables["Yn"][i].keys()):
                self.variables["Yn"][i][adjacency_id] = self.gm.addVar(lb=1.0, vtype=g.GRB.INTEGER, name="y_{{{i},{aid}}}".format(i=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for adjacency_id in list(self.variables["P"][i].keys()):
                self.variables["P"][i][adjacency_id]["AA"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{clone_id},{aid},AA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id]["AB"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{clone_id},{aid},AB}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id]["BA"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{clone_id},{aid},BA}}".format(clone_id=str(i), aid=str(adjacency_id)))
                self.variables["P"][i][adjacency_id]["BB"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{clone_id},{aid},BB}}".format(clone_id=str(i), aid=str(adjacency_id)))
        for i in self.clone_ids:
            for vertex_id in list(self.variables["T"][i].keys()):
                self.variables["T"][i][vertex_id]["A"] = self.gm.addVar(vtype=g.GRB.INTEGER, name="t_{{{clone_id},{p_id},A}}".format(clone_id=str(i), p_id=str(vertex_id)))
                self.variables["T"][i][vertex_id]["B"] = self.gm.addVar(vtype=g.GRB.INTEGER, name="t_{{{clone_id},{p_id},B}}".format(clone_id=str(i), p_id=str(vertex_id)))
        for vertex_id in list(self.variables["M"].keys()):
            self.variables["M"][vertex_id]["A"] = self.gm.addVar(vtype=g.GRB.INTEGER, name="m_{{{p_id},A}}".format(p_id=str(vertex_id)))
            self.variables["M"][vertex_id]["B"] = self.gm.addVar(vtype=g.GRB.INTEGER, name="m_{{{p_id},B}}".format(p_id=str(vertex_id)))
        for adjacency_id in list(self.variables["Pc"].keys()):
            self.variables["Pc"][adjacency_id]["AA"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},AA}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id]["AB"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},AB}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id]["BA"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},BA}}".format(aid=str(adjacency_id)))
            self.variables["Pc"][adjacency_id]["BB"] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid},BB}}".format(aid=str(adjacency_id)))
        for i in self.clone_ids:
            for a_id in list(self.variables["Prod"]["PY"][i].keys()):
                for ph in ["AA", "AB", "BA", "BB"]:
                    self.variables["Prod"]["PY"][i][a_id][ph] = self.gm.addVar(vtype=g.GRB.INTEGER, name="py_{{{clone_id},{a_id},{phasing}}}".format(clone_id=str(i),
                                                                                                                                                     a_id=a_id,
                                                                                                                                                     phasing=ph))
        #####
        #
        # creating constraints
        #
        ####
        # bigM constraints to linearize integer-binary product of adjacency phasing and respective adjacency copy numbers
        for i in self.clone_ids:
            for adjacency in self.adjacencies:
                for ph in ["AA", "AB", "BA", "BB"]:
                    self.gm.addConstr(self.big_ms_by_a_ids[adjacency.idx] * self.variables["P"][i][adjacency.idx][ph],
                                      g.GRB.GREATER_EQUAL,
                                      self.variables["Prod"]["PY"][i][adjacency.idx][ph])
                    self.gm.addConstr(self.variables["Prod"]["PY"][i][adjacency.idx][ph], g.GRB.GREATER_EQUAL, 0.0)
                    if adjacency.adjacency_type == AdjacencyType.NOVEL:
                        a_cn_var = self.variables["Yn"][i][adjacency.idx]
                    else:
                        a_cn_var = self.variables["Yr"][i][adjacency.idx][ph]
                    self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][i][adjacency.idx][ph])
                    self.gm.addConstr(a_cn_var - self.big_ms_by_a_ids[adjacency.idx] * (1 - self.variables["P"][i][adjacency.idx][ph]),
                                      g.GRB.LESS_EQUAL,
                                      self.variables["Prod"]["PY"][i][adjacency.idx][ph])
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
                    # processing allele "A"
                    ####
                    t_ta_lin_exp = g.LinExpr()
                    t_ta_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_a_segmentcn.cn))
                    t_ta_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_a_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["AA"]), mult=-1.0)
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["BA"]), mult=-1.0)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["AA"] * self.variables["P"][i][ref_adjacency.idx]["AA"]), mult=-1.0)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["BA"] * self.variables["P"][i][ref_adjacency.idx]["BA"]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if u == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AA"]), mult=na_multiplier)
                        # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["AA"]), mult=na_multiplier)
                        if beg:
                            t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                            # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                        else:
                            t_ta_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                            # t_ta_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                    self.gm.addConstr(self.variables["T"][i][u]["A"], g.GRB.EQUAL, t_ta_lin_exp)
                    self.gm.addConstr(t_ta_lin_exp >= 0, name="non_negative_telomere_score_{{{clone_id},{vertex},A}}".format(clone_id=str(i), vertex=str(u)))
                    ####
                    # processing allele "B"
                    ####
                    t_tb_lin_exp = g.LinExpr()
                    t_tb_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_b_segmentcn.cn))
                    t_tb_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_b_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["AB"]), mult=-1.0)
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["BB"]), mult=-1.0)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["AB"] * self.variables["P"][i][ref_adjacency.idx]["AB"]), mult=-1.0)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["BB"] * self.variables["P"][i][ref_adjacency.idx]["BB"]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if u == novel_adjacency.position1:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BB"]), mult=na_multiplier)
                        # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["BB"]), mult=na_multiplier)
                        if beg:
                            t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                            # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                        else:
                            t_tb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                            # t_tb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                    self.gm.addConstr(self.variables["T"][i][u]["B"], g.GRB.EQUAL, t_tb_lin_exp)
                    self.gm.addConstr(t_tb_lin_exp >= 0, name="non_negative_telomere_score_{{{clone_id},{vertex},B}}".format(clone_id=str(i), vertex=str(u)))
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
                    # processing allele "A"
                    ###
                    t_ha_lin_exp = g.LinExpr()
                    t_ha_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_a_segmentcn.cn))
                    t_ha_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_a_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["AA"]), mult=-1.0)
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["AB"]), mult=-1.0)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["AA"] * self.variables["P"][i][ref_adjacency.idx]["AA"]), mult=-1.0)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["AB"] * self.variables["P"][i][ref_adjacency.idx]["AB"]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if v == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 <= novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AA"]), mult=na_multiplier)
                        # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["AA"]), mult=na_multiplier)
                        if beg:
                            t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                            # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                        else:
                            t_ha_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                            # t_ha_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                    self.gm.addConstr(self.variables["T"][i][v]["A"], g.GRB.EQUAL, t_ha_lin_exp)
                    self.gm.addConstr(t_ha_lin_exp >= 0, name="non_negative_telomere_score_{{{clone_id},{vertex},A}}".format(clone_id=str(i), vertex=str(v)))
                    ####
                    # processing allele "B"
                    ####
                    t_hb_lin_exp = g.LinExpr()
                    t_hb_lin_exp.add(g.LinExpr(self.variables["B"][i][segment.idx] * self.scnr_by_sids[segment.idx].maj_b_segmentcn.cn))
                    t_hb_lin_exp.add(g.LinExpr((1 - self.variables["B"][i][segment.idx]) * self.scnr_by_sids[segment.idx].min_b_segmentcn.cn))
                    for ref_adjacency in reference_adjacencies:
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["BA"]), mult=-1.0)
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][ref_adjacency.idx]["BB"]), mult=-1.0)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["BA"] * self.variables["P"][i][ref_adjacency.idx]["BA"]), mult=-1.0)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yr"][i][ref_adjacency.idx]["BB"] * self.variables["P"][i][ref_adjacency.idx]["BB"]), mult=-1.0)
                    for novel_adjacency in novel_adjacencies:
                        if v == novel_adjacency.position1.idx:
                            beg = novel_adjacency.position1 < novel_adjacency.position2
                        else:
                            beg = novel_adjacency.position2 < novel_adjacency.position1
                        na_multiplier = -1.0 if novel_adjacency.position1 != novel_adjacency.position2 else -2.0
                        t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BB"]), mult=na_multiplier)
                        # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i]["BB"]), mult=-1.0)
                        if beg:
                            t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["BA"]), mult=na_multiplier)
                            # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i]["BA"]))
                        else:
                            t_hb_lin_exp.add(g.LinExpr(self.variables["Prod"]["PY"][i][novel_adjacency.idx]["AB"]), mult=na_multiplier)
                            # t_hb_lin_exp.add(g.LinExpr(self.variables["Yn"][i][novel_adjacency.idx] * self.variables["P"][i]["AB"]))
                    self.gm.addConstr(self.variables["T"][i][v]["B"], g.GRB.EQUAL, t_hb_lin_exp)
                    self.gm.addConstr(t_hb_lin_exp >= 0, name="non_negative_telomere_score_{{{clone_id},{vertex},B}}".format(clone_id=str(i), vertex=str(v)))
        novel_adjacencies = [data["object"] for u, v, data in self.iag.nov_adjacency_edges(data=True)]
        for na in novel_adjacencies:
            for i in self.clone_ids:
                self.gm.addConstr(g.quicksum([self.variables["P"][i][na.idx][ph] for ph in ["AA", "AB", "BA", "BB"]]), g.GRB.LESS_EQUAL, 1.0)
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx]["AA"], [self.variables["P"][i][na.idx]["AA"] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx]["AB"], [self.variables["P"][i][na.idx]["AB"] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx]["BA"], [self.variables["P"][i][na.idx]["BA"] for i in self.clone_ids])
            self.gm.addGenConstrOr(self.variables["Pc"][na.idx]["BB"], [self.variables["P"][i][na.idx]["BB"] for i in self.clone_ids])
            #####
            # requiring every NA to be present at least once in some of clones
            #####
            self.gm.addConstr(self.variables["Pc"][na.idx]["AA"] + self.variables["Pc"][na.idx]["AB"] + self.variables["Pc"][na.idx]["BA"] + self.variables["Pc"][na.idx]["BB"],
                              g.GRB.EQUAL,
                              1)
            # NA involving the same segments must be on the same alleles
            s1 = self.iag.get_segment_edge(node=na.position1.idx, data=True)[2]["object"]
            s2 = self.iag.get_segment_edge(node=na.position2.idx, data=True)[2]["object"]
            if s1 == s2:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][na.idx]["AB"], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][na.idx]["BA"], g.GRB.EQUAL, 0)

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
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["BA"], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["BB"], g.GRB.EQUAL, 0)
            else:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["AA"] + self.variables["P"][i][ra.idx]["AB"], g.GRB.LESS_EQUAL, 1)
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["BA"] + self.variables["P"][i][ra.idx]["BB"], g.GRB.LESS_EQUAL, 1)
            if not s2cnr.allele_specific:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["AB"], g.GRB.EQUAL, 0)
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["BB"], g.GRB.EQUAL, 0)
            else:
                for i in self.clone_ids:
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["AA"] + self.variables["P"][i][ra.idx]["BA"], g.GRB.LESS_EQUAL, 1)
                    self.gm.addConstr(self.variables["P"][i][ra.idx]["AB"] + self.variables["P"][i][ra.idx]["BB"], g.GRB.LESS_EQUAL, 1)

        for p_id in self.positions_idxs:
            a_extremities = [self.variables["T"][i][p_id]["A"] for i in self.clone_ids]
            b_extremities = [self.variables["T"][i][p_id]["B"] for i in self.clone_ids]
            if self.config["telomere_max_allele_specific"]:
                self.gm.addGenConstrMax(self.variables["M"][p_id]["A"], a_extremities)
                self.gm.addGenConstrMax(self.variables["M"][p_id]["B"], b_extremities)
            else:
                self.gm.addGenConstrMax(self.variables["M"][p_id]["A"], a_extremities + b_extremities)
                self.gm.addConstr(self.variables["M"][p_id]["B"], g.GRB.EQUAL, self.variables["M"][p_id]["A"])

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
        objective_lin_exp.add(g.quicksum([self.variables["M"][p_id]["A"] for p_id in self.positions_idxs]), mult=scalar_coefficient)
        objective_lin_exp.add(g.quicksum([self.variables["M"][p_id]["B"] for p_id in self.positions_idxs]), mult=scalar_coefficient)
        self.gm.setObjective(objective_lin_exp, g.GRB.MINIMIZE)

    def solve_gurobi_model(self):
        self.gm.optimize()

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
