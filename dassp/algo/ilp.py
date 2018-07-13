from collections import defaultdict

from dassp.core.graph import IntervalAdjacencyGraph
from dassp.core.structures import SegmentCopyNumberProfile, AdjacencyCopyNumberProfile, check_and_fill_segments_to_fragments
from dassp.core.structures import Phasing, AdjacencyType, Haplotype, get_aabb_for_ra, get_abba_for_na_and_position
import gurobi as g


class SatModelSingleClone(object):
    def __init__(self,
                 hapl_segments,
                 hapl_adjacencies,
                 scnp,
                 hapl_telomeres):
        self.hapl_segments = hapl_segments
        self.hapl_adjacencies = hapl_adjacencies
        self.scnp = scnp
        self.hapl_telomeres = hapl_telomeres
        self.iag = IntervalAdjacencyGraph(segments=self.hapl_segments, adjacencies=self.hapl_adjacencies)
        self.iag.build_graph()
        self.big_Ms_by_a_ids = self.compute_big_ms()
        self.variables = {
            # binary variable for assigning A/B labeling for measured allele-specific segment copy number pairs (c', c'')
            "B": {s.stable_id_non_hap: None for s in self.hapl_segments},
            # integer copy number for diploid reference adjacencies
            "YR": {a.stable_id_non_phased: {Phasing.AA: None,
                                            Phasing.AB: None,
                                            Phasing.BA: None,
                                            Phasing.BB: None}
                   for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies)},
            # A/B assignment (phasing) for all adjacencies, both reference and novel. Reference will be constrained hard for AA, BB
            "P": {a.stable_id_non_phased: {Phasing.AA: None,
                                           Phasing.AB: None,
                                           Phasing.BA: None,
                                           Phasing.BB: None}
                  for a in self.hapl_adjacencies},
            # integer copy number for diploid counterpart for a single measured haploid novel adjacency
            #   since only one such counterpart can exist, a single variable is used
            "YN": {a.stable_id_non_phased: None
                   for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies)},
            #
            "Prod": {
                "PY": {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in self.hapl_adjacencies},
            },
        }
        self.model = g.Model("DASSP_sat_single_clone")
        self.gm = self.model

    def compute_big_ms(self):
        result = {}
        for adjacency in self.hapl_adjacencies:
            s1 = self.iag.get_segment_edge(node=adjacency.position1, data=True)[2]["object"]
            s2 = self.iag.get_segment_edge(node=adjacency.position2, data=True)[2]["object"]
            s1_all_segment_copy_number_values = [self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s1, haplotype=Haplotype.A),
                                                 self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s1, haplotype=Haplotype.B)]
            s2_all_segment_copy_number_values = [self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s2, haplotype=Haplotype.A),
                                                 self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s2, haplotype=Haplotype.B)]
            big_m = max(s1_all_segment_copy_number_values + s2_all_segment_copy_number_values)
            result[adjacency.stable_id_non_phased] = big_m
        return result

    def build_gurobi_model(self):
        self.define_variables()
        self.define_constraints()

    def define_variables(self):

        for segment_id in list(self.variables["B"].keys()):
            self.variables["B"][segment_id] = self.gm.addVar(vtype=g.GRB.BINARY, name="b_{{{j}}}".format(j=str(segment_id)))

        for adjacency_id in list(self.variables["YR"].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.variables["YR"][adjacency_id][ph] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="y_{{r,{aid},{phas}}}".format(aid=str(adjacency_id), phas=str(ph)))

        for adjacency_id in list(self.variables["YN"].keys()):
            self.variables["YN"][adjacency_id] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="y_{{n,{aid}}}".format(aid=str(adjacency_id)))

        for adjacency_id in list(self.variables["P"].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.variables["P"][adjacency_id][ph] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{aid}, {phas}}}".format(aid=str(adjacency_id), phas=str(ph)))

        for adjacency_id in list(self.variables["Prod"]["PY"].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.variables["Prod"]["PY"][adjacency_id][ph] = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER, name="py_{{{a_id},{phas}}}".format(a_id=str(adjacency_id), phas=str(ph)))

    def define_constraints(self):
        self.define_big_m_constraints_for_adjacencies_copy_numbers()
        self.define_constraints_for_phasing()
        self.define_constraints_on_nodes()

    def solve_model(self):
        self.gm.optimize()

    def define_big_m_constraints_for_adjacencies_copy_numbers(self):
        for adjacency in self.hapl_adjacencies:
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                a_id = adjacency.stable_id_non_phased
                big_m = max(1, self.big_Ms_by_a_ids[a_id])
                self.gm.addConstr(big_m * self.variables["P"][a_id][ph], g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][a_id][ph])
                self.gm.addConstr(self.variables["Prod"]["PY"][a_id][ph], g.GRB.GREATER_EQUAL, 0.0)
                if adjacency.adjacency_type == AdjacencyType.NOVEL:
                    a_cn_var = self.variables["YN"][a_id]
                else:
                    a_cn_var = self.variables["YR"][a_id][ph]
                self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][a_id][ph])
                self.gm.addConstr(a_cn_var - big_m * (1 - self.variables["P"][a_id][ph]), g.GRB.LESS_EQUAL, self.variables["Prod"]["PY"][a_id][ph])

    def define_constraints_for_phasing(self):
        self.define_constraints_for_ref_phasing()
        self.define_constraints_for_novel_phasing()

    def define_constraints_for_novel_phasing(self):
        # for novel adjacency need to force exactly one phasing to exist and for self-loop adjacencies to block (force 0) the AB/BA phasing
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            a_id = adjacency.stable_id_non_phased
            self.gm.addConstr(g.quicksum([self.variables["P"][a_id][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]), g.GRB.EQUAL, 1,
                              name="ph_nov_uniq_dip_{{{aid}}}".format(aid=str(a_id)))
            if adjacency.is_self_loop_hapl:
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables["P"][a_id][ph], g.GRB.EQUAL, 0, name="ph_nov_sl_{{{aid},{phas}}}".format(aid=str(a_id), phas=str(ph)))
        # reciprocal phasing forced on the same haplotype
        for u, v, data in self.iag.ref_adjacency_edges(data=True):
            u_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=u))
            v_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=v))
            if len(u_na_edges_w_data) == 0 or len(v_na_edges_w_data) == 0:
                continue
            if len(u_na_edges_w_data) > 1 or len(v_na_edges_w_data) > 1:
                raise Exception()
            u_na = u_na_edges_w_data[0][2]["object"]
            u_na_id = u_na.stable_id_non_phased
            v_na = v_na_edges_w_data[0][2]["object"]
            v_na_id = v_na.stable_id_non_phased
            for h in [Haplotype.A, Haplotype.B]:
                self.gm.addConstr(self.variables["P"][u_na_id][get_aabb_for_ra(haplotype=h)] +
                                  self.variables["P"][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=h)] -
                                  self.variables["P"][v_na_id][get_aabb_for_ra(haplotype=h)] -
                                  self.variables["P"][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=h)],
                                  g.GRB.EQUAL,
                                  0,
                                  name="recip_nov_{{{u},{v},{u_aid},{v_aid},{hap}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id), hap=str(h)))

    def define_constraints_for_ref_phasing(self):
        # for reference adjacecnies only need to block (force 0) AB/BA phasing
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies):
            a_id = adjacency.stable_id_non_phased
            for ph in [Phasing.AB, Phasing.BA]:
                self.gm.addConstr(self.variables["P"][a_id][ph], g.GRB.EQUAL, 0, name="ph_{{r,{aid},{phas}}}".format(aid=str(a_id), phas=str(ph)))

    def define_constraints_on_nodes(self):
        for node, data in self.iag.nodes(data=True):
            for haplotype in [Haplotype.A, Haplotype.B]:
                lin_expr = self.get_lin_expr_for_node_and_haplotype(node=node, haplotype=haplotype)
                if node in self.hapl_telomeres:
                    self.gm.addConstr(lin_expr, g.GRB.GREATER_EQUAL, 1, name="cne_{{{v},{hap}}}".format(v=str(node), hap=str(haplotype)))
                else:
                    self.gm.addConstr(lin_expr, g.GRB.EQUAL, 0, name="cnb_{{{v},{hap}}}".format(v=str(node), hap=str(haplotype)))

    def get_lin_expr_for_node_and_haplotype(self, node, haplotype):
        result = g.LinExpr()
        s_u, s_v, data = self.iag.get_segment_edge(node=node, data=True)
        s = data["object"]
        result.add(self.get_lin_expr_for_haplotype_scn(segment=s, haplotype=haplotype))
        ref_adj_edges_w_data = self.iag.ref_adjacency_edges(data=True, nbunch=node)
        ref_adjs = [data["object"] for _, __, data in ref_adj_edges_w_data]
        assert len(ref_adjs) <= 1
        for adjacency in ref_adjs:
            result.add(g.LinExpr(self.variables["Prod"]["PY"][adjacency.stable_id_non_phased][self.get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
        nov_adj_edges_w_data = self.iag.nov_adjacency_edges(data=True, nbunch=node)
        nov_adjs = [data["object"] for _, __, data in nov_adj_edges_w_data]
        for adjacency in nov_adjs:
            self_loop = adjacency.is_self_loop_hapl
            if self_loop:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][adjacency.stable_id_non_phased][self.get_aabb_for_ra(haplotype=haplotype)]), mult=-2)
            else:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][adjacency.stable_id_non_phased][self.get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
                result.add(g.LinExpr(self.variables["Prod"]["PY"][adjacency.stable_id_non_phased][self.get_abba_for_na_and_position(novel_adjacency=adjacency,
                                                                                                                                    position=node, haplotype=haplotype)]), mult=-1)
        return result

    def get_lin_expr_for_haplotype_scn(self, segment, haplotype):
        if haplotype == Haplotype.A:
            result = g.LinExpr(self.variables["B"][segment.stable_id_non_hap] * self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A))
            result.add(g.LinExpr((1 - self.variables["B"][segment.stable_id_non_hap]) * self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B)))
            return result
        elif haplotype == Haplotype.B:
            result = g.LinExpr(self.variables["B"][segment.stable_id_non_hap] * self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B))
            result.add(g.LinExpr((1 - self.variables["B"][segment.stable_id_non_hap]) * self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A)))
            return result
        else:
            raise Exception()

    def get_abba_for_na_and_position(self, novel_adjacency, position, haplotype):
        left = novel_adjacency.position1.stable_id_non_hap if novel_adjacency.is_sorted_non_phased else novel_adjacency.position2.stable_id_non_hap
        if position.stable_id_non_hap == left:
            if haplotype == Haplotype.A:
                return Phasing.AB
            elif haplotype == Haplotype.B:
                return Phasing.BA
            else:
                raise Exception()
        else:
            if haplotype == Haplotype.A:
                return Phasing.BA
            elif haplotype == Haplotype.B:
                return Phasing.AB
            else:
                raise Exception()

    def get_aabb_for_ra(self, haplotype):
        if haplotype == Haplotype.A:
            return Phasing.AA
        elif haplotype == Haplotype.B:
            return Phasing.BB
        else:
            raise Exception()

    def get_scnp_from_model(self):
        result = SegmentCopyNumberProfile()
        for s in self.hapl_segments:
            ha, hb = (Haplotype.A, Haplotype.B) if self.variables["B"][s.stable_id_non_hap].X == 1 else (Haplotype.B, Haplotype.A)
            cn_a = self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s, haplotype=ha)
            cn_b = self.scnp.get_hap_aware_cn_by_seg_and_hap(segment=s, haplotype=hb)
            result.set_cn_record_for_segment(segment=s, cn=cn_a, haplotype=Haplotype.A)
            result.set_cn_record_for_segment(segment=s, cn=cn_b, haplotype=Haplotype.B)
        return result

    def get_acnp_from_model(self):
        result = AdjacencyCopyNumberProfile()
        for adj in self.hapl_adjacencies:
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                cn = self.variables["Prod"]["PY"][adj.stable_id_non_phased][ph].X
                result.set_cnt_record_for_adjacency(adjacency=adj, cn=cn, phasing=ph)
        return result


class SatModelMultiClone(object):
    def __init__(self,
                 hapl_segments,
                 hapl_adjacencies,
                 clone_specific_scnp,
                 hapl_telomeres):
        self.hapl_segments = hapl_segments
        self.hapl_adjacencies = hapl_adjacencies
        self.clone_specific_scnp = clone_specific_scnp
        self.hapl_telomeres = hapl_telomeres
        self.clone_ids = sorted(clone_specific_scnp.keys())
        self.iag = IntervalAdjacencyGraph(segments=hapl_segments, adjacencies=self.hapl_adjacencies)
        self.iag.build_graph()
        self.clone_specific_big_Ms_by_a_ids = self.compute_big_ms()
        self.variables = {
            "B": {s.stable_id_non_hap: None for s in self.hapl_segments},
            "YR": {
                clone_id: {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies)
                } for clone_id in self.clone_ids},
            "P": {
                clone_id: {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in self.hapl_adjacencies
                } for clone_id in self.clone_ids},
            "YN": {
                clone_id: {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies)
                } for clone_id in self.clone_ids},
            "Prod": {
                "PY": {
                    clone_id: {
                        a.stable_id_non_phased: {Phasing.AA: None,
                                                 Phasing.AB: None,
                                                 Phasing.BA: None,
                                                 Phasing.BB: None}
                        for a in self.hapl_adjacencies
                    } for clone_id in self.clone_ids},
                "Pp": {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in self.hapl_adjacencies},
            }
        }
        self.model = g.Model("DASSP_sat_multi_clone")
        self.gm = self.model

    def compute_big_ms(self):
        result = defaultdict(dict)
        for clone_id in self.clone_ids:
            for adjacency in self.hapl_adjacencies:
                s1 = self.iag.get_segment_edge(node=adjacency.position1, data=True)[2]["object"]
                s2 = self.iag.get_segment_edge(node=adjacency.position2, data=True)[2]["object"]
                s1_all_segment_copy_numbers = [self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s1, haplotype=Haplotype.A),
                                               self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s1, haplotype=Haplotype.B)]
                s2_all_segment_copy_numbers = [self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s2, haplotype=Haplotype.A),
                                               self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s2, haplotype=Haplotype.B)]
                big_m = max(s1_all_segment_copy_numbers + s2_all_segment_copy_numbers)
                result[clone_id][adjacency.stable_id_non_phased] = big_m
        return result

    def build_gurobi_model(self):
        self.define_variables()
        self.define_constraints()

    def define_variables(self):
        for segment_id in list(self.variables["B"].keys()):
            self.variables["B"][segment_id] = self.gm.addVar(vtype=g.GRB.BINARY, name="b_{{{j}}}".format(j=str(segment_id)))

        for clone_id in self.clone_ids:
            for adjacency_id in list(self.variables["YR"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["YR"][clone_id][adjacency_id][ph] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="y_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                                          aid=str(adjacency_id),
                                                                                                                                                          phas=str(ph)))

        for clone_id in self.clone_ids:
            for adjacency_id in list(self.variables["YN"][clone_id].keys()):
                self.variables["YN"][clone_id][adjacency_id] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="y_{{n,{cid},{aid}}}".format(cid=str(clone_id),
                                                                                                                                           aid=str(adjacency_id)))

        for clone_id in self.clone_ids:
            for adjacency_id in list(self.variables["P"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["P"][clone_id][adjacency_id][ph] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                                aid=str(adjacency_id),
                                                                                                                                                phas=str(ph)))

        for clone_id in self.clone_ids:
            for adjacency_id in list(self.variables["Prod"]["PY"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["Prod"]["PY"][clone_id][adjacency_id][ph] = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER,
                                                                                              name="py_{{{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                      aid=str(adjacency_id),
                                                                                                                                      phas=str(ph)))
        for adjacency_id in list(self.variables["Prod"]["Pp"].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.variables["Prod"]["Pp"][adjacency_id][ph] = self.gm.addVar(vtype=g.GRB.BINARY, name="pp_{{{aid},{phas}}}".format(aid=str(adjacency_id),
                                                                                                                                      phas=str(ph)))

    def define_constraints(self):
        self.define_big_m_constraints_for_adjacency_copy_numbers()
        self.define_constraints_for_phasing()
        self.define_constraints_on_nodes()

    def solve_model(self):
        self.gm.optimize()

    def define_big_m_constraints_for_adjacency_copy_numbers(self):
        for clone_id in self.clone_ids:
            for adjacency in self.hapl_adjacencies:
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    a_id = adjacency.stable_id_non_phased
                    big_m = max(1, self.clone_specific_big_Ms_by_a_ids[clone_id][a_id])
                    self.gm.addConstr(big_m * self.variables["P"][clone_id][a_id][ph], g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][clone_id][a_id][ph])
                    self.gm.addConstr(self.variables["Prod"]["PY"][clone_id][a_id][ph], g.GRB.GREATER_EQUAL, 0.0)
                    if adjacency.adjacency_type == AdjacencyType.NOVEL:
                        a_cn_var = self.variables["YN"][clone_id][a_id]
                    else:
                        a_cn_var = self.variables["YR"][clone_id][a_id][ph]
                    self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][clone_id][a_id][ph])
                    self.gm.addConstr(a_cn_var - big_m * (1 - self.variables["P"][clone_id][a_id][ph]), g.GRB.LESS_EQUAL, self.variables["Prod"]["PY"][clone_id][a_id][ph])

    def define_constraints_for_phasing(self):
        self.define_constraints_for_ref_phasing()
        self.define_constraints_for_novel_adjacencies()

    def define_constraints_for_ref_phasing(self):
        for clone_id in self.clone_ids:
            for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies):
                a_id = adjacency.stable_id_non_phased
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables["P"][clone_id][a_id][ph], g.GRB.EQUAL, 0, name="ph_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                         aid=str(a_id),
                                                                                                                                         phas=str(ph)))

    def define_constraints_for_novel_adjacencies(self):
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            a_id = adjacency.stable_id_non_phased
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                # p_e is an or over clone specific p_{i,e}
                self.gm.addGenConstrOr(self.variables["Prod"]["Pp"][a_id][ph],
                                       [self.variables["P"][clone_id][a_id][ph] for clone_id in self.clone_ids],
                                       name="pp_{{{aid},{phas}}}".format(aid=str(a_id), phas=str(ph)))
            # exactly one p_e must be equal to 1
            self.gm.addConstr(g.quicksum([self.variables["Prod"]["Pp"][a_id][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]),
                              g.GRB.EQUAL, 1, name="ph_nov_uniq_dip_{{{aid}}}".format(aid=str(a_id)))
            if adjacency.is_self_loop_hapl:
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables["Prod"]["Pp"][a_id][ph], g.GRB.EQUAL, 0)

        for u, v, data in self.iag.ref_adjacency_edges(data=True):
            u_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=u))
            v_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=v))
            if len(u_na_edges_w_data) == 0 or len(v_na_edges_w_data) == 0:
                continue
            if len(u_na_edges_w_data) > 1 or len(v_na_edges_w_data) > 1:
                raise Exception()
            u_na = u_na_edges_w_data[0][2]["object"]
            u_na_id = u_na.stable_id_non_phased
            v_na = v_na_edges_w_data[0][2]["object"]
            v_na_id = v_na.stable_id_non_phased
            for h in [Haplotype.A, Haplotype.B]:
                self.gm.addConstr(self.variables["Prod"]["Pp"][u_na_id][get_aabb_for_ra(haplotype=h)] +
                                  self.variables["Prod"]["Pp"][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=h)] -
                                  self.variables["Prod"]["Pp"][v_na_id][get_aabb_for_ra(haplotype=h)] -
                                  self.variables["Prod"]["Pp"][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=h)],
                                  g.GRB.EQUAL,
                                  0,
                                  name="recip_nov_{{{u},{v},{u_aid},{v_aid},{hap}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id), hap=str(h)))

    def define_constraints_on_nodes(self):
        for clone_id in self.clone_ids:
            for node, data in self.iag.nodes(data=True):
                for haplotype in [Haplotype.A, Haplotype.B]:
                    lin_expr = self.get_lin_expr_for_node_haplotype(clone_id=clone_id, node=node, haplotype=haplotype)
                    if node in self.hapl_telomeres:
                        self.gm.addConstr(lin_expr, g.GRB.GREATER_EQUAL, 1, name="cne_{{{v},{hap}}}".format(v=str(node), hap=str(haplotype)))
                    else:
                        self.gm.addConstr(lin_expr, g.GRB.EQUAL, 0, name="cnb_{{{v},{hap}}}".format(v=str(node), hap=str(haplotype)))

    def get_lin_expr_for_node_haplotype(self, clone_id, node, haplotype):
        result = g.LinExpr()
        s_u, s_v, data = self.iag.get_segment_edge(node=node, data=True)
        s = data["object"]
        result.add(self.get_lin_expr_for_haplotype_scn(clone_id=clone_id, segment=s, haplotype=haplotype))
        ref_adj_edges_w_data = self.iag.ref_adjacency_edges(data=True, nbunch=node)
        ref_adjs = [data["object"] for _, __, data in ref_adj_edges_w_data]
        assert len(ref_adjs) <= 1
        for adjacency in ref_adjs:
            result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
        nov_adj_edges_w_data = self.iag.nov_adjacency_edges(data=True, nbunch=node)
        nov_adjs = [data["object"] for _, __, data in nov_adj_edges_w_data]
        for adjacency in nov_adjs:
            self_loop = adjacency.is_self_loop_hapl
            if self_loop:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-2)
            else:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_abba_for_na_and_position(novel_adjacency=adjacency,
                                                                                                                                         position=node,
                                                                                                                                         haplotype=haplotype)]), mult=-1)
        return result

    def get_lin_expr_for_haplotype_scn(self, clone_id, segment, haplotype):
        if haplotype == Haplotype.A:
            result = g.LinExpr(self.variables["B"][segment.stable_id_non_hap] * self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment,
                                                                                                                                                   haplotype=Haplotype.A))
            result.add(g.LinExpr((1 - self.variables["B"][segment.stable_id_non_hap]) * self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment,
                                                                                                                                                           haplotype=Haplotype.B)))
            return result
        elif haplotype == Haplotype.B:
            result = g.LinExpr(self.variables["B"][segment.stable_id_non_hap] * self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment,
                                                                                                                                                   haplotype=Haplotype.B))
            result.add((1 - self.variables["B"][segment.stable_id_non_hap]) * self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment,
                                                                                                                                                 haplotype=Haplotype.A))
            return result
        else:
            raise Exception()

    def get_clone_specific_scnp_from_model(self):
        result = {}
        for clone_id in self.clone_ids:
            scnp = SegmentCopyNumberProfile()
            result[clone_id] = scnp
            for s in self.hapl_segments:
                ha, hb = (Haplotype.A, Haplotype.B) if int(round(self.variables["B"][s.stable_id_non_hap].X)) == 1 else (Haplotype.B, Haplotype.A)
                cn_a = self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s, haplotype=ha)
                cn_b = self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=s, haplotype=hb)
                scnp.set_cn_record_for_segment(segment=s, cn=cn_a, haplotype=Haplotype.A)
                scnp.set_cn_record_for_segment(segment=s, cn=cn_b, haplotype=Haplotype.B)
        return result

    def get_clone_specific_acnp_from_models(self):
        result = {}
        for clone_id in self.clone_ids:
            acnp = AdjacencyCopyNumberProfile()
            result[clone_id] = acnp
            for adj in self.hapl_adjacencies:
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    cn = int(round(self.variables["Prod"]["PY"][clone_id][adj.stable_id_non_phased][ph].X))
                    acnp.set_cnt_record_for_adjacency(adjacency=adj, cn=cn, phasing=ph)
        return result


class OptModelMultiClone(object):
    def __init__(self,
                 hapl_segments,
                 hapl_adjacencies,
                 clone_specific_scnp,
                 hapl_telomeres,
                 scn_boundaries,
                 hapl_adjacencies_groups=None,
                 hapl_segments_to_fragments=None,
                 tree=None,
                 hapl_nov_adjacencies_fp=0.0):
        self.scn_boundaries = scn_boundaries
        self.tree = tree,
        self.hapl_segments = hapl_segments
        self.hapl_adjacencies = hapl_adjacencies
        self.hapl_nov_adjacencies_fp = hapl_nov_adjacencies_fp
        self.clone_specific_scnp = clone_specific_scnp
        self.hapl_telomeres = hapl_telomeres
        hapl_adjacencies_groups = hapl_adjacencies_groups if hapl_adjacencies_groups is not None else {}
        self.hapl_adjacencies_groups = hapl_adjacencies_groups
        hapl_segments_to_fragments = check_and_fill_segments_to_fragments(segments=hapl_segments, segments_to_fragments=hapl_segments_to_fragments)
        self.fragments_ids = sorted(set(hapl_segments_to_fragments.values()))
        self.hapl_segments_to_fragments = hapl_segments_to_fragments
        self.clone_ids = sorted(self.clone_specific_scnp.keys())
        self.iag = IntervalAdjacencyGraph(segments=self.hapl_segments, adjacencies=self.hapl_adjacencies)
        self.iag.build_graph()
        self.variables = self.populate_variables_dict()
        self.model = g.Model("DASSP")
        self.gm = self.model

    def populate_variables_dict(self):
        result = {}
        # defining b_f variables
        result["B"] = {fid: None for fid in sorted(self.fragments_ids)}
        # defining c_{i,j,A}, c_{i,j,B} variables
        result["C"] = {
            clone_id: {
                s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in self.hapl_segments
            } for clone_id in self.clone_ids
        }
        result["YR"] = {
            clone_id: {
                a.stable_id_non_phased: {Phasing.AA: None,
                                         Phasing.AB: None,
                                         Phasing.BA: None,
                                         Phasing.BB: None}
                for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies)
            } for clone_id in self.clone_ids
        }
        # defining N'_{i,a} for every unlabeled novel adjacency a
        result["YN"] = {
            clone_id: {
                a.stable_id_non_phased: None for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies)
            } for clone_id in self.clone_ids
        }
        # defining p_{i,e} variables
        result["P"] = {
            clone_id: {
                a.stable_id_non_phased: {Phasing.AA: None,
                                         Phasing.AB: None,
                                         Phasing.BA: None,
                                         Phasing.BB: None}
                for a in self.hapl_adjacencies
            } for clone_id in self.clone_ids
        }
        # defining p_{i,u} variables
        result["U"] = {
            clone_id: {
                key: None for key in sorted(self.hapl_adjacencies_groups.keys())
            } for clone_id in self.clone_ids
        }
        # defining \Delta_{i,j,A}, \Delta_{i,j,B} variables
        result["Delta"] = {
            clone_id: {
                s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in self.hapl_segments
            } for clone_id in self.clone_ids
        }
        result["Prod"] = {
            # defining N_{i,e} variables
            "PY": {
                clone_id: {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in self.hapl_adjacencies
                } for clone_id in self.clone_ids
            },
            # defining the p_{a} variables
            "Pp": {
                a.stable_id_non_phased: {Phasing.AA: None,
                                         Phasing.AB: None,
                                         Phasing.BA: None,
                                         Phasing.BB: None}
                for a in self.hapl_adjacencies
            }
        }
        return result

    def build_gurobi_model(self):
        self.define_variables()
        self.define_constraints()
        self.define_objective()

    def solve_model(self):
        self.gm.optimize()

    def define_variables(self):
        for fid in list(self.variables["B"].keys()):
            self.variables["B"][fid] = self.gm.addVar(vtype=g.GRB.BINARY, name="b_{{fid}}".format(fid=str(fid)))

        for clone_id in self.clone_ids:
            for aid in list(self.variables["YR"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["YR"][clone_id][aid][ph] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="N'_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                                  aid=str(aid),
                                                                                                                                                  phas=str(ph)))

        for clone_id in self.clone_ids:
            for sid in list(self.variables["C"][clone_id].keys()):
                for h in [Haplotype.A, Haplotype.B]:
                    lower = self.scn_boundaries[clone_id][sid][h][0]
                    upper = self.scn_boundaries[clone_id][sid][h][1]
                    self.variables["C"][clone_id][sid][h] = self.gm.addVar(lb=lower, ub=upper, vtype=g.GRB.INTEGER, name="c_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                                                          sid=str(sid),
                                                                                                                                                          hap=str(h)))

        for clone_id in self.clone_ids:
            for aid in list(self.variables["YN"][clone_id].keys()):
                self.variables["YN"][clone_id][aid] = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="N'_{{n,{cid},{aid}}}".format(cid=str(clone_id),
                                                                                                                                   aid=str(aid)))

        for clone_id in self.clone_ids:
            for aid in list(self.variables["P"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["P"][clone_id][aid][ph] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{a,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                         aid=str(aid),
                                                                                                                                         phas=str(ph)))
        for clone_id in self.clone_ids:
            for gid in list(self.variables["U"][clone_id].keys()):
                self.variables["U"][clone_id][gid] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{u,{cid},{gid}}}".format(cid=str(clone_id),
                                                                                                                          gid=str(gid)))
        for clone_id in self.clone_ids:
            for sid in list(self.variables["Delta"][clone_id].keys()):
                for h in [Haplotype.A, Haplotype.B]:
                    self.variables["Delta"][clone_id][sid][h] = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER, name="delta_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                                                    sid=str(sid),
                                                                                                                                                    hap=str(h)))

        for clone_id in self.clone_ids:
            for aid in list(self.variables["Prod"]["PY"][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.variables["Prod"]["PY"][clone_id][aid][ph] = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER,
                                                                                     name="N_{{{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                            aid=str(aid),
                                                                                                                            phas=str(ph)))

        for aid in list(self.variables["Prod"]["Pp"].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.variables["Prod"]["Pp"][aid][ph] = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{a,{aid},{phas}}}".format(aid=str(aid),
                                                                                                                              phas=str(ph)))

    def define_constraints(self):
        self.define_upper_bound_constraints_on_segment_copy_numbers()
        self.define_big_m_constraints_for_adjacency_copy_numbers()
        self.define_constraints_for_phasing()
        self.define_constraints_for_nov_adjacency_groups()
        self.define_constraints_on_nodes()
        self.define_constraints_on_deltas()

    def define_upper_bound_constraints_on_segment_copy_numbers(self):
        pass
        # for clone_id in self.clone_ids:
        #     for sid in list(self.variables["C"][clone_id].keys()):
        #         for h in [Haplotype.A, Haplotype.B]:
        #             if isinstance(self.scn_boundaries, int):
        #                 c_max = self.scn_boundaries
        #             else:
        #                 c_max = self.scn_boundaries[clone_id][sid][h][1]
        #             self.gm.addConstr(self.variables["C"][clone_id][sid][h], g.GRB.LESS_EQUAL, c_max)

    def define_big_m_constraints_for_adjacency_copy_numbers(self):
        for clone_id in self.clone_ids:
            for adjacency in self.hapl_adjacencies:
                u, v = self.iag.get_edge_vertices_pair_from_adjacency(adjacency=adjacency)
                s1u, s1v, s1data = self.iag.get_segment_edge(node=u, data=True)
                s1 = s1data["object"]
                s1id = s1.stable_id_non_hap
                s2u, s2v, s2data = self.iag.get_segment_edge(node=v, data=True)
                s2 = s2data["object"]
                s2id = s2.stable_id_non_hap
                s1_cn_boundaries = [self.scn_boundaries[clone_id][s1id][Haplotype.A][1],
                                    self.scn_boundaries[clone_id][s1id][Haplotype.B][1]]
                s2_cn_boundaries = [self.scn_boundaries[clone_id][s2id][Haplotype.A][1],
                                    self.scn_boundaries[clone_id][s2id][Haplotype.B][1]]
                c_max = max(s1_cn_boundaries + s2_cn_boundaries)
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    aid = adjacency.stable_id_non_phased
                    self.gm.addConstr(c_max * self.variables["P"][clone_id][aid][ph], g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][clone_id][aid][ph])
                    self.gm.addConstr(self.variables["Prod"]["PY"][clone_id][aid][ph], g.GRB.GREATER_EQUAL, 0.0)
                    if adjacency.adjacency_type == AdjacencyType.NOVEL:
                        a_cn_var = self.variables["YN"][clone_id][aid]
                    else:
                        a_cn_var = self.variables["YR"][clone_id][aid][ph]
                    self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables["Prod"]["PY"][clone_id][aid][ph])
                    self.gm.addConstr(a_cn_var - c_max * (1 - self.variables["P"][clone_id][aid][ph]), g.GRB.LESS_EQUAL,
                                      self.variables["Prod"]["PY"][clone_id][aid][ph])

    def define_constraints_for_phasing(self):
        self.define_constraints_for_ref_phasing()
        self.define_constraints_for_novel_phasing()

    def define_constraints_for_ref_phasing(self):
        for clone_id in self.clone_ids:
            for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies):
                aid = adjacency.stable_id_non_phased
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables["P"][clone_id][aid][ph], g.GRB.EQUAL, 0, name="ph_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                        aid=str(aid),
                                                                                                                                        phas=str(ph)))

    def define_constraints_for_novel_phasing(self):
        fp_lin_expr = g.LinExpr()
        hapl_nov_adjs = list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies))
        hapl_nov_adjs_cnt = len(hapl_nov_adjs)
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            aid = adjacency.stable_id_non_phased
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.gm.addGenConstrOr(self.variables["Prod"]["Pp"][aid][ph],
                                       [self.variables["P"][clone_id][aid][ph] for clone_id in self.clone_ids],
                                       name="pp_{{{aid},{phas}}}".format(aid=str(aid), phas=str(ph)))
            # force presence of a diploid counterpart adjacency for an observed unlabeled one
            self.gm.addConstr(g.quicksum([self.variables["Prod"]["Pp"][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]),
                              g.GRB.LESS_EQUAL, 1, name="pp-nov-uniq-dip_{{{aid}}}".format(aid=str(aid)))
            fp_lin_expr.add(g.quicksum([self.variables["Prod"]["Pp"][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]), mult=(1.0 / hapl_nov_adjs_cnt))
            if adjacency.is_self_loop_hapl:
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables["Prod"]["Pp"][aid][ph], g.GRB.EQUAL, 0)

        for u, v, data in self.iag.ref_adjacency_edges(data=True):
            u_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=u))
            v_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=v))
            if len(u_na_edges_w_data) == 0 or len(v_na_edges_w_data) == 0:
                continue
            if len(u_na_edges_w_data) > 1 or len(v_na_edges_w_data) > 1:
                continue
            u_na = u_na_edges_w_data[0][2]["object"]
            u_na_id = u_na.stable_id_non_phased
            v_na = v_na_edges_w_data[0][2]["object"]
            v_na_id = v_na.stable_id_non_phased
            expr = g.LinExpr()
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][u_na_id][get_aabb_for_ra(haplotype=Haplotype.A)]), mult=1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.A)]), mult=1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][v_na_id][get_aabb_for_ra(haplotype=Haplotype.A)]), mult=-1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.A)]), mult=-1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][v_na_id][get_aabb_for_ra(haplotype=Haplotype.B)]), mult=1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.B)]), mult=1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][u_na_id][get_aabb_for_ra(haplotype=Haplotype.B)]), mult=-1)
            expr.add(g.LinExpr(self.variables["Prod"]["Pp"][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.B)]), mult=-1)
            self.gm.addConstr(expr, g.GRB.GREATER_EQUAL, -1, name="recip_nov_l_{{{u},{v},{u_aid},{v_aid}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))
            self.gm.addConstr(expr, g.GRB.LESS_EQUAL, 1, name="recip_nov_u_{{{u},{v},{u_aid},{v_aid}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))
            # for h in [Haplotype.A, Haplotype.B]:
            #     self.gm.addConstr(self.variables["Prod"]["Pp"][u_na_id][get_aabb_for_ra(haplotype=h)] +
            #                       self.variables["Prod"]["Pp"][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=h)] -
            #                       self.variables["Prod"]["Pp"][v_na_id][get_aabb_for_ra(haplotype=h)] -
            #                       self.variables["Prod"]["Pp"][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=h)],
            #                       g.GRB.EQUAL,
            #                       0,
            #                       name="recip_nov_{{{u},{v},{u_aid},{v_aid},{hap}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id), hap=str(h)))

        self.gm.addConstr(fp_lin_expr, g.GRB.GREATER_EQUAL, g.LinExpr(1 - self.hapl_nov_adjacencies_fp), name="nov-adj-fp")

    def define_constraints_for_nov_adjacency_groups(self):
        for gid, group in self.hapl_adjacencies_groups.items():
            for clone_id in self.clone_ids:
                group_var = self.variables["U"][clone_id][gid]
                group_size = len(group.adjacencies)
                lin_expr = g.LinExpr()
                for adjacency in group.adjacencies:
                    aid = adjacency.stable_id_non_phased
                    self.gm.addConstr(group_var, g.GRB.LESS_EQUAL,
                                      g.quicksum([self.variables["P"][clone_id][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]),
                                      name="group-AND-upper_{{{cid},{gid},{aid}}}".format(cid=str(clone_id),
                                                                                          gid=str(gid),
                                                                                          aid=str(aid)))
                    for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                        lin_expr.add(self.variables["P"][clone_id][aid][ph])
                lin_expr.add(g.LinExpr(group_size), mult=-1)
                lin_expr.add(g.LinExpr(1))
                self.gm.addConstr(group_var, g.GRB.GREATER_EQUAL, lin_expr, name="group-AND-lower_{{{cid},{gid}}}".format(cid=str(clone_id),
                                                                                                                          gid=str(gid)))
            self.gm.addConstr(g.quicksum([self.variables["U"][clone_id][gid] for clone_id in self.clone_ids]),
                              g.GRB.GREATER_EQUAL, 1, name="group-across_{{{gid}}}".format(gid=str(gid)))

    def define_constraints_on_nodes(self):
        for clone_id in self.clone_ids:
            for node, data in self.iag.nodes(data=True):
                for haplotype in [Haplotype.A, Haplotype.B]:
                    lin_exp = self.get_lin_expr_for_node_haplotype(clone_id=clone_id, node=node, haplotype=haplotype)
                    if node in self.hapl_telomeres:
                        self.gm.addConstr(lin_exp, g.GRB.GREATER_EQUAL, 0, name="cne_{{{cid},{v},{hap}}}".format(cid=str(clone_id),
                                                                                                                 v=str(node), hap=str(haplotype)))
                    else:
                        self.gm.addConstr(lin_exp, g.GRB.EQUAL, 0, name="cnb_{{{cid},{v},{hap}}}".format(cid=str(clone_id), v=str(node), hap=str(haplotype)))

    def get_lin_expr_for_node_haplotype(self, clone_id, node, haplotype):
        result = g.LinExpr()
        s_u, s_v, data = self.iag.get_segment_edge(node=node, data=True)
        s = data["object"]
        sid = s.stable_id_non_hap
        s_var = self.variables["C"][clone_id][sid][haplotype]
        result.add(s_var)
        ref_adj_edges_w_data = self.iag.ref_adjacency_edges(data=True, nbunch=node)
        ref_adjs = [data["object"] for _, __, data in ref_adj_edges_w_data]
        assert len(ref_adjs) <= 1
        for adjacency in ref_adjs:
            result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
        nov_adj_edges_w_data = self.iag.nov_adjacency_edges(data=True, nbunch=node)
        nov_adjs = [data["object"] for _, __, data in nov_adj_edges_w_data]
        for adjacency in nov_adjs:
            self_loop = adjacency.is_self_loop_hapl
            if self_loop:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-2)
            else:
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
                result.add(g.LinExpr(self.variables["Prod"]["PY"][clone_id][adjacency.stable_id_non_phased][get_abba_for_na_and_position(novel_adjacency=adjacency,
                                                                                                                                         position=node,
                                                                                                                                         haplotype=haplotype)]), mult=-1)
        return result

    def define_constraints_on_deltas(self):
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                fid = self.hapl_segments_to_fragments[sid]
                f_var = self.variables["B"][fid]
                s_cn_a_var = self.variables["C"][clone_id][sid][Haplotype.A]
                s_cn_b_var = self.variables["C"][clone_id][sid][Haplotype.B]
                s_cn_a = self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A)
                s_cn_b = self.clone_specific_scnp[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B)
                delta_var_a = self.variables["Delta"][clone_id][sid][Haplotype.A]
                delta_var_b = self.variables["Delta"][clone_id][sid][Haplotype.B]
                s_a_lin_expr = g.LinExpr(s_cn_a_var - f_var * s_cn_a - (1 - f_var) * s_cn_b)
                self.gm.addConstr(delta_var_a, g.GRB.GREATER_EQUAL, s_a_lin_expr, name="delta-u_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                              sid=str(sid),
                                                                                                                              hap=str(Haplotype.A)))
                self.gm.addConstr(delta_var_a, g.GRB.GREATER_EQUAL, -1 * s_a_lin_expr, name="delta-l_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                                   sid=str(sid),
                                                                                                                                   hap=str(Haplotype.A)))
                s_b_lin_exp = g.LinExpr(s_cn_b_var - (1 - f_var) * s_cn_a - f_var * s_cn_b)
                self.gm.addConstr(delta_var_b, g.GRB.GREATER_EQUAL, s_b_lin_exp, name="delta-u_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                             sid=str(sid),
                                                                                                                             hap=str(Haplotype.B)))
                self.gm.addConstr(delta_var_b, g.GRB.GREATER_EQUAL, -1 * s_b_lin_exp, name="delta-l_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                                                                                                                                  sid=str(sid),
                                                                                                                                  hap=str(Haplotype.B)))
                # self.gm.addGenConstrAbs(delta_var_a, s_cn_a_var - f_var * s_cn_a - (1 - f_var) * s_cn_b, name="delta_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                #                                                                                                                                    sid=str(sid),
                #                                                                                                                                    hap=str(Haplotype.A)))
                # self.gm.addGenConstrAbs(delta_var_b, s_cn_b_var - (1 - f_var) * s_cn_a - f_var * s_cn_b, name="delta_{{{cid},{sid},{hap}}}".format(cid=str(clone_id),
                #                                                                                                                                    sid=str(sid),
                #                                                                                                                                    hap=str(Haplotype.B)))

    def define_objective(self):
        lin_exp = g.LinExpr()
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                length = segment.length_1000
                delta_var_a = self.variables["Delta"][clone_id][sid][Haplotype.A]
                delta_var_b = self.variables["Delta"][clone_id][sid][Haplotype.B]
                lin_exp.add(g.LinExpr(delta_var_a + delta_var_b), mult=length)
        self.gm.setObjective(lin_exp, g.GRB.MINIMIZE)

    def get_clone_specific_scnp_from_model(self):
        result = {}
        for clone_id in self.clone_ids:
            scnp = SegmentCopyNumberProfile()
            result[clone_id] = scnp
            for s in self.hapl_segments:
                sid = s.stable_id_non_hap
                cn_a = int(round(self.variables["C"][clone_id][sid][Haplotype.A].X))
                cn_b = int(round(self.variables["C"][clone_id][sid][Haplotype.B].X))
                scnp.set_cn_record_for_segment(segment=s, cn=cn_a, haplotype=Haplotype.A)
                scnp.set_cn_record_for_segment(segment=s, cn=cn_b, haplotype=Haplotype.B)
        return result

    def get_clone_specific_acnp_from_model(self):
        result = {}
        for clone_id in self.clone_ids:
            acnp = AdjacencyCopyNumberProfile()
            result[clone_id] = acnp
            for adj in self.hapl_adjacencies:
                aid = adj.stable_id_non_phased
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    cn = int(round(self.variables["Prod"]["PY"][clone_id][aid][ph].X))
                    acnp.set_cnt_record_for_adjacency(adjacency=adj, cn=cn, phasing=ph)
        return result
