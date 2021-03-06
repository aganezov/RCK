import gurobi as g

from rck.core.graph import IntervalAdjacencyGraph
from rck.core.io import FALSE_POSITIVE, EXTERNAL_NA_ID, AG_LABELING
from rck.core.structures import Phasing, AdjacencyType, Haplotype, get_aabb_for_ra, get_abba_for_na_and_position
from rck.core.structures import SegmentCopyNumberProfile, AdjacencyCopyNumberProfile, check_and_fill_segments_to_fragments, AdjacencyGroupType, CNBoundaries
from rck.utils.scn.process import get_haploid_scnt

FRAGMENT_ALLELE = "fragment_flipping"
SEGMENT_COPY_NUMBER = "segment_copy_number"
YR = "YR"
YN = "YN"
P = "P"
ADJ_GROUPS_M = "U-m"
ADJ_GROUPS_L = "U-l"
DELTA = "Delta"
PROD = "Prod"
PY = "PY"
PP = "Pp"

DEFAULT_GROUP_M_FP = "DEFAULT_GROUP_M_FP"
DEFAULT_GROUP_N_FP = "DEFAULT_GROUP_N_FP"
SEGMENT_LENGTH_ATTRIBUTE = "SEGMENT_LENGTH_ATTRIBUTE"


def merge_variables_from_presolve(*variables_dicts, result=None):
    if result is None:
        result = {}
    for variables_dict in variables_dicts:
        if FRAGMENT_ALLELE in variables_dicts:
            if FRAGMENT_ALLELE not in result:
                result[FRAGMENT_ALLELE] = {}
            for key, value in variables_dict[FRAGMENT_ALLELE].items():
                result[key] = value.X
        if SEGMENT_COPY_NUMBER in variables_dict:
            if SEGMENT_COPY_NUMBER not in result:
                result[SEGMENT_COPY_NUMBER] = {}
            result_dict = result[SEGMENT_COPY_NUMBER]
            for clone_id, value in variables_dict[SEGMENT_COPY_NUMBER].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for sid, clone_specific_dict in value.items():
                    result_dict[clone_id][sid] = {
                        Haplotype.A: clone_specific_dict[Haplotype.A].X,
                        Haplotype.B: clone_specific_dict[Haplotype.B].X
                    }
        if YR in variables_dict:
            if YR not in result:
                result[YR] = {}
            result_dict = result[YR]
            for clone_id, value in variables_dict[YR].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for aid, clone_specific_dict in value.items():
                    result_dict[clone_id][aid] = {
                        Phasing.AA: clone_specific_dict[Phasing.AA].X,
                        Phasing.AB: clone_specific_dict[Phasing.AB].X,
                        Phasing.BA: clone_specific_dict[Phasing.BA].X,
                        Phasing.BB: clone_specific_dict[Phasing.BB].X
                    }
        if YN in variables_dict:
            if YN not in result:
                result[YN] = {}
            result_dict = result[YN]
            for clone_id, value in variables_dict[YN].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for aid, cs_value in value.items():
                    result_dict[clone_id][aid] = cs_value.X
        if P in variables_dict:
            if P not in result:
                result[P] = {}
            result_dict = result[P]
            for clone_id, value in variables_dict[P].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for aid, clone_specific_dict in value.items():
                    result_dict[clone_id][aid] = {
                        Phasing.AA: clone_specific_dict[Phasing.AA].X,
                        Phasing.AB: clone_specific_dict[Phasing.AB].X,
                        Phasing.BA: clone_specific_dict[Phasing.BA].X,
                        Phasing.BB: clone_specific_dict[Phasing.BB].X
                    }
        if ADJ_GROUPS_M in variables_dict:
            if ADJ_GROUPS_M not in result:
                result[ADJ_GROUPS_M] = {}
            result_dict = result[ADJ_GROUPS_M]
            for clone_id, value in variables_dict[ADJ_GROUPS_M].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for gid, cs_value in value.items():
                    result_dict[clone_id][gid] = cs_value.X
        if ADJ_GROUPS_L in variables_dict:
            if ADJ_GROUPS_L not in result:
                result[ADJ_GROUPS_L] = {}
            result_dict = result[ADJ_GROUPS_L]
            for gid, value in variables_dict[ADJ_GROUPS_L].items():
                result_dict[gid] = {
                    Haplotype.A: value[Haplotype.A].X,
                    Haplotype.B: value[Haplotype.B].X
                }
        if DELTA in variables_dict:
            if DELTA not in result:
                result[DELTA] = {}
            result_dict = result[DELTA]
            for clone_id, value in variables_dict[DELTA].items():
                if clone_id not in result_dict:
                    result_dict[clone_id] = {}
                for sid, cs_value in value.items():
                    result_dict[clone_id][sid] = {
                        Haplotype.A: cs_value[Haplotype.A].X,
                        Haplotype.B: cs_value[Haplotype.B].X
                    }
        if PROD in variables_dict:
            if PROD not in result:
                result[PROD] = {}
            result_dict = result[PROD]
            if PY in variables_dict[PROD]:
                if PY not in result_dict:
                    result_dict[PY] = {}
                result_dict = result_dict[PY]
                for clone_id, value in variables_dict[PROD][PY].items():
                    if clone_id not in result_dict:
                        result_dict[clone_id] = {}
                    for aid, cs_value in value.items():
                        result_dict[clone_id][aid] = {
                            Phasing.AA: cs_value[Phasing.AA].X,
                            Phasing.AB: cs_value[Phasing.AB].X,
                            Phasing.BA: cs_value[Phasing.BA].X,
                            Phasing.BB: cs_value[Phasing.BB].X,
                        }
            result_dict = result[PROD]
            if PP in variables_dict[PROD]:
                if PP not in result_dict:
                    result_dict[PP] = {}
                result_dict = result_dict[PP]
                for aid, value in variables_dict[PROD][PP].items():
                    result_dict[aid] = {
                        Phasing.AA: value[Phasing.AA].X,
                        Phasing.AB: value[Phasing.AB].X,
                        Phasing.BA: value[Phasing.BA].X,
                        Phasing.BB: value[Phasing.BB].X
                    }
    return result


class OptModelMultiClone(object):
    def __init__(self,
                 hapl_segments,
                 hapl_adjacencies,
                 scnt,
                 hapl_telomeres,
                 scnb,
                 hapl_adjacencies_groups=None,
                 hapl_segments_to_fragments=None,
                 hapl_nov_adjacencies_fp=0.0,
                 starting_vars=None,
                 solve_as_haploid=False,
                 extra=None):
        self.scnb = scnb
        self.hapl_segments = hapl_segments
        self.hapl_adjacencies = hapl_adjacencies
        self.hapl_nov_adjacencies_fp = hapl_nov_adjacencies_fp
        self.solve_as_haploid = solve_as_haploid
        self.scnt = scnt
        if self.solve_as_haploid:
            self.scnt = get_haploid_scnt(segments=self.hapl_segments, scnt=self.scnt)
        self.hapl_telomeres = hapl_telomeres
        self.hapl_adjacencies_groups = hapl_adjacencies_groups if hapl_adjacencies_groups is not None else []
        hapl_segments_to_fragments = check_and_fill_segments_to_fragments(segments=hapl_segments, segments_to_fragments=hapl_segments_to_fragments)
        self.fragments_ids = sorted(set(hapl_segments_to_fragments.values()))
        self.hapl_segments_to_fragments = hapl_segments_to_fragments
        self.clone_ids = sorted(self.scnt.keys())
        self.extra = extra if extra is not None else {}
        self.starting_vars = starting_vars
        self.iag = IntervalAdjacencyGraph(segments=self.hapl_segments, adjacencies=self.hapl_adjacencies)
        self.iag.build_graph()
        self.variables = self.populate_variables_dict()
        self.model = g.Model("RCK-mc-mln")  # multi-clone, molecule, labeling, and general groups; change when other features (e.g., multi-sample, labeling constraints, trees, etc)
        self.gm = self.model
        self.reciprocal_locations = []

    def populate_variables_dict(self):
        result = {}

        # defining b_f variables, encoding allele "flipping"
        result[FRAGMENT_ALLELE] = {fid: None for fid in sorted(self.fragments_ids)}

        # defining c_{i,j,A}, c_{i,j,B} variables encoding segment copy numbers
        result[SEGMENT_COPY_NUMBER] = {
            clone_id: {
                s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in self.hapl_segments
            } for clone_id in self.clone_ids
        }

        # defining Y_r for every reference adjacency
        result[YR] = {
            clone_id: {
                a.stable_id_non_phased: {Phasing.AA: None,
                                         Phasing.AB: None,
                                         Phasing.BA: None,
                                         Phasing.BB: None}
                for a in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies)
            } for clone_id in self.clone_ids
        }

        # defining N'_{i,a} for every unlabeled novel adjacency a, as its unique copy number

        result[YN] = {
            clone_id: {
                a.stable_id_non_phased: None for a in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies)
            } for clone_id in self.clone_ids
        }
        # defining p_{i,e} variables (phasing for every adjacency)
        result[P] = {
            clone_id: {
                a.stable_id_non_phased: {Phasing.AA: None,
                                         Phasing.AB: None,
                                         Phasing.BA: None,
                                         Phasing.BB: None}
                for a in self.hapl_adjacencies
            } for clone_id in self.clone_ids
        }

        # defining p_{i,u} variables,
        # molecule adjacency groups
        result[ADJ_GROUPS_M] = {
            clone_id: {
                ag.gid: None for ag in sorted(filter(lambda ag: ag.group_type == AdjacencyGroupType.MOLECULE, self.hapl_adjacencies_groups), key=lambda ag: ag.gid)
            } for clone_id in self.clone_ids
        }

        # defining p_{A|B,u} variable
        # labeling adjacency groups
        result[ADJ_GROUPS_L] = {
            ag.gid: {Haplotype.A: None, Haplotype.B: None} for ag in
            sorted(filter(lambda ag: ag.group_type == AdjacencyGroupType.LABELING, self.hapl_adjacencies_groups), key=lambda ag: ag.gid)
        }

        # defining \Delta_{i,j,A}, \Delta_{i,j,B} variables (distance from inferred segment copy numbers, and the starting ones)
        result[DELTA] = {
            clone_id: {
                s.stable_id_non_hap: {Haplotype.A: None, Haplotype.B: None} for s in self.hapl_segments
            } for clone_id in self.clone_ids
        }

        # helper variables used to encode product of other binary/integer variables
        result[PROD] = {
            # defining N_{i,e} variables (actual copy number of every phased realization of an adjacency)
            PY: {
                clone_id: {
                    a.stable_id_non_phased: {Phasing.AA: None,
                                             Phasing.AB: None,
                                             Phasing.BA: None,
                                             Phasing.BB: None}
                    for a in self.hapl_adjacencies
                } for clone_id in self.clone_ids
            },

            # defining the p_{a} variables (presence/absence of adjacency (any of its phasing realizations) across clones)
            PP: {
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
        ###
        # Binary variables that define whether inferred haplotype copy number pairing for a given segment/fragments in synch (1) on the allele-specific input order, or flips it (0)
        ###
        for fid in list(self.variables[FRAGMENT_ALLELE].keys()):
            var = self.gm.addVar(vtype=g.GRB.BINARY, name="b_{{{fid}}}".format(fid=str(fid)))
            if self.starting_vars is not None and FRAGMENT_ALLELE in self.starting_vars and fid in self.starting_vars[FRAGMENT_ALLELE]:
                var.start = self.starting_vars[FRAGMENT_ALLELE][fid]
            self.variables[FRAGMENT_ALLELE][fid] = var
        ######

        ###
        # Integer copy number of a phased reference adjacency (AB, BA variations will be set to 0 in constraints definitions, via respective binary indicators)
        #   lower bound is set at 1, as the absence (i.e., copy number 0) is achieved via a binary indicator variable.
        ###
        for clone_id in self.clone_ids:
            for aid in list(self.variables[YR][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    var = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="N'_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id), aid=str(aid), phas=str(ph)))
                    if self.starting_vars is not None and \
                            YR in self.starting_vars and clone_id in self.starting_vars[YR] and aid in self.starting_vars[YR][clone_id] and \
                            ph in self.starting_vars[YR][clone_id][aid]:
                        var.start = self.starting_vars[YR][clone_id][aid][ph]
                    self.variables[YR][clone_id][aid][ph] = var
        ######

        ###
        # Integer segment copy numbers, bounded by lower and upper bounds, provided as part of the input
        ###
        start_cnt = 0
        for clone_id in self.clone_ids:
            for sid in list(self.variables[SEGMENT_COPY_NUMBER][clone_id].keys()):
                min_lower = min([self.scnb[clone_id].get_cnb(sid=sid, hap=h, boundary_type=CNBoundaries.LOWER) for h in [Haplotype.A, Haplotype.B]])
                max_upper = max([self.scnb[clone_id].get_cnb(sid=sid, hap=h, boundary_type=CNBoundaries.UPPER) for h in [Haplotype.A, Haplotype.B]])
                for h in [Haplotype.A, Haplotype.B]:
                    var = self.gm.addVar(lb=min_lower, ub=max_upper, vtype=g.GRB.INTEGER, name="c_{{{cid},{sid},{hap}}}".format(cid=str(clone_id), sid=str(sid), hap=str(h)))
                    if self.starting_vars is not None and \
                            SEGMENT_COPY_NUMBER in self.starting_vars and clone_id in self.starting_vars[SEGMENT_COPY_NUMBER] and sid in self.starting_vars[SEGMENT_COPY_NUMBER][
                        clone_id] and \
                            h in self.starting_vars[SEGMENT_COPY_NUMBER][clone_id][sid]:
                        var.start = self.starting_vars[SEGMENT_COPY_NUMBER][clone_id][sid][h]
                        start_cnt += 1
                    self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][h] = var
        # print("segment start cnt = {cnt}".format(cnt=start_cnt))
        ######

        ###
        # Integer copy number of a novel adjacency (single phased realization). Every unphased novel adjacency can have at most one phased realization, thus 1 variable.
        #   lower bound is set at 1, as absence (i.e., copy number 0) is achieved via a binary indicator variable.
        ###
        start_cnt = 0
        for clone_id in self.clone_ids:
            for aid in list(self.variables[YN][clone_id].keys()):
                var = self.gm.addVar(lb=1, vtype=g.GRB.INTEGER, name="N'_{{n,{cid},{aid}}}".format(cid=str(clone_id), aid=str(aid)))
                if self.starting_vars is not None and \
                        YN in self.starting_vars and clone_id in self.starting_vars[YN] and aid in self.starting_vars[YN][clone_id]:
                    var.start = self.starting_vars[YN][clone_id][aid]
                    start_cnt += 1
                self.variables[YN][clone_id][aid] = var
        # print("yn start cnt = {cnt}".format(cnt=start_cnt))
        ######

        ###
        # Binary indicator variable that determines presence/absence of phased adjacency realizations
        ###
        for clone_id in self.clone_ids:
            for aid in list(self.variables[P][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    var = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{a,{cid},{aid},{phas}}}".format(cid=str(clone_id), aid=str(aid), phas=str(ph)))
                    if self.starting_vars is not None and \
                            P in self.starting_vars and clone_id in self.starting_vars[P] and aid in self.starting_vars[P][clone_id] and \
                            ph in self.starting_vars[P][clone_id][aid]:
                        var.start = self.starting_vars[P][clone_id][aid][ph]
                    self.variables[P][clone_id][aid][ph] = var
        ######

        ###
        # Binary indicator presence variable that determines whether a group of adjacencies is present in a given clone, or not
        ###
        start_cnt = 0
        for clone_id in self.clone_ids:
            for gid in list(self.variables[ADJ_GROUPS_M][clone_id].keys()):
                var = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{u,{cid},{gid}}}".format(cid=str(clone_id), gid=str(gid)))
                if self.starting_vars is not None and \
                        ADJ_GROUPS_M in self.starting_vars and clone_id in self.starting_vars[ADJ_GROUPS_M] and gid in self.starting_vars[ADJ_GROUPS_M][clone_id]:
                    var.start = self.starting_vars[ADJ_GROUPS_M][clone_id][gid]
                    start_cnt += 1
                self.variables[ADJ_GROUPS_M][clone_id][gid] = var
        # print("adj group m start cnt = {cnt}".format(cnt=start_cnt))
        ######

        ###
        # Binary indicator presence variable that determines whether any of the adjacencies in the given group are present on respective haplotype
        ###
        start_cnt = 0
        for gid in list(self.variables[ADJ_GROUPS_L].keys()):
            for haplotype in [Haplotype.A, Haplotype.B]:
                var = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{u,{gid},{hap}}}".format(gid=str(gid), hap=str(haplotype.value)))
                if self.starting_vars is not None and \
                        ADJ_GROUPS_L in self.starting_vars and gid in self.starting_vars[ADJ_GROUPS_L] and haplotype in self.starting_vars[ADJ_GROUPS_L][gid]:
                    var.start = self.starting_vars[ADJ_GROUPS_L][gid][haplotype]
                    start_cnt += 1
                self.variables[ADJ_GROUPS_L][gid][haplotype] = var
        # print("adj group l start cnt = {cnt}".format(cnt=start_cnt))
        ###
        # Internal variable that encodes the difference between the inferred segment copy numbers and the original values
        ###
        for clone_id in self.clone_ids:
            for sid in list(self.variables[DELTA][clone_id].keys()):
                for h in [Haplotype.A, Haplotype.B]:
                    var = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER, name="delta_{{{cid},{sid},{hap}}}".format(cid=str(clone_id), sid=str(sid), hap=str(h)))
                    if self.starting_vars is not None and \
                            DELTA in self.starting_vars and clone_id in self.starting_vars[DELTA] and sid in self.starting_vars[DELTA][clone_id] and \
                            h in self.starting_vars[DELTA][clone_id][sid]:
                        var.start = self.starting_vars[DELTA][clone_id][sid][h]
                    self.variables[DELTA][clone_id][sid][h] = var
        ######

        ###
        # Internal integer variable encoding internal copy number for all phased realizations of the adjacencies. Used in conjunction with respective binary variable
        #   and bigM notation technique to encode the true copy number (i.e., the result)
        ###
        start_cnt = 0
        for clone_id in self.clone_ids:
            for aid in list(self.variables[PROD][PY][clone_id].keys()):
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    var = self.gm.addVar(lb=0, vtype=g.GRB.INTEGER, name="N_{{{cid},{aid},{phas}}}".format(cid=str(clone_id), aid=str(aid), phas=str(ph)))
                    if self.starting_vars is not None and \
                            PROD in self.starting_vars and PY in self.starting_vars[PROD] and clone_id in self.starting_vars[PROD][PY] and \
                            aid in self.starting_vars[PROD][PY][clone_id] and ph in self.starting_vars[PROD][PY][clone_id][aid]:
                        var.start = self.starting_vars[PROD][PY][clone_id][aid][ph]
                        start_cnt += 1
                    self.variables[PROD][PY][clone_id][aid][ph] = var
        # print("prod py start cnt = {cnt}".format(cnt=start_cnt))
        ######

        ###
        # Binary variables used for enforcing phasing constraints across clones
        ###
        start_cnt = 0
        for aid in list(self.variables[PROD][PP].keys()):
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                var = self.gm.addVar(vtype=g.GRB.BINARY, name="p_{{a,{aid},{phas}}}".format(aid=str(aid), phas=str(ph)))
                if self.starting_vars is not None and \
                        PROD in self.starting_vars and PP in self.starting_vars[PROD] and aid in self.starting_vars[PROD][PP] and ph in self.starting_vars[PROD][PP][aid]:
                    var.start = self.starting_vars[PROD][PP][aid][ph]
                    start_cnt += 1
                self.variables[PROD][PP][aid][ph] = var
        # print("prod pp start cnt = {cnt}".format(cnt=start_cnt))

    def define_constraints(self):
        self.define_segment_copy_number_boundary_constraints()
        self.define_big_m_constraints_for_adjacencies_copy_numbers()
        self.define_constraints_nov_adjacency_overall_presence()
        self.define_constraints_for_labeling()
        self.define_constraints_for_adjacency_groups()
        self.define_constraints_on_nodes()
        self.define_constraints_on_deltas()
        if self.solve_as_haploid:
            self.define_allele_flipping_haploid_constraints()

    def define_allele_flipping_haploid_constraints(self):
        """
        When solving a haploid version of the problem, we don't need allele flipping possibility enabled.
        We thus fix the respecting binary variables at 1 (i.e., cn on haplotype A equal to allele A in the input),
            which coupled with the "haploidazation" of the input SCNT in the init, as well as forced 0-0 bounds on haplotype B to have the CN on segment B always at 0,
            and not contributing to the problem solution at all.

        """
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                fid = self.hapl_segments_to_fragments[sid]
                f_var = self.variables[FRAGMENT_ALLELE][fid]
                self.gm.addConstr(f_var, g.GRB.EQUAL, 1, name="allele-flipping-fix-haploid-{{{cid},{fid}}}".format(cid=clone_id, fid=str(fid)))

    def define_segment_copy_number_boundary_constraints(self):
        """
        Lower upper copy number boundaries on the inferred segment copy numbers must be in sync with the flipping nature of the input allele-specific segment copy number boundaries

        Example input allele-specific copy number values (e.g., cn={'A': 1, 'B': 10))
        If the input has copy number boundaries of 0-2  for A and 9-12 for B the inferred segment copy number variables (both for A and B) are set to have lower = 0 (min) and upper = 12 (max)
        But when inferring we need to have copy numbers if
        (i) no-flipping alleles, then (i.e., fragment binary indicator = 1), then we'll need to have 0 <= A <= 2 and 9 <= B <= 12
        (ii) flipping alleles, then (i.e., fragments binary indicator = 0), then we'll need to have 9 <= A <= 12 and 0 <= B <= 2
        We achive this by constructing a constraint that utilize the flipping binary indicator variable to nullify the incorrect cn boundary value.
        """
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                fid = self.hapl_segments_to_fragments[sid]
                f_var = self.variables[FRAGMENT_ALLELE][fid]
                s_cn_a_var = self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.A]
                s_cn_b_var = self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.B]
                ####
                # allele-tied copy number boundaries
                ####
                lower_a = self.scnb[clone_id].get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.LOWER)
                lower_b = self.scnb[clone_id].get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.LOWER)
                upper_a = self.scnb[clone_id].get_cnb(sid=sid, hap=Haplotype.A, boundary_type=CNBoundaries.UPPER)
                upper_b = self.scnb[clone_id].get_cnb(sid=sid, hap=Haplotype.B, boundary_type=CNBoundaries.UPPER)
                ####
                # IF we are solving a haploid version of the problem, we must force the segment copy numbers on haplotype B to be equal to 0
                ####
                if self.solve_as_haploid:
                    lower_b, upper_b = 0, 0
                ####
                self.gm.addConstr(s_cn_a_var, g.GRB.LESS_EQUAL, f_var * upper_a + (1 - f_var) * upper_b, name="scnb-A-upper-{{{cid},{sid}}}".format(cid=clone_id, sid=sid))
                self.gm.addConstr(s_cn_a_var, g.GRB.GREATER_EQUAL, f_var * lower_a + (1 - f_var) * lower_b, name="scnb-A-lower-{{{cid},{sid}}}".format(cid=clone_id, sid=sid))
                self.gm.addConstr(s_cn_b_var, g.GRB.LESS_EQUAL, (1 - f_var) * upper_a + f_var * upper_b, name="scnb-B-upper-{{{cid},{sid}}}".format(cid=clone_id, sid=sid))
                self.gm.addConstr(s_cn_b_var, g.GRB.GREATER_EQUAL, (1 - f_var) * lower_a + f_var * lower_b, name="scnb-B-lower-{{{cid},{sid}}}".format(cid=clone_id, sid=sid))

    def define_big_m_constraints_for_adjacencies_copy_numbers(self):
        for clone_id in self.clone_ids:
            for adjacency in self.hapl_adjacencies:
                u, v = self.iag.get_edge_vertices_pair_from_adjacency(adjacency=adjacency)
                s1u, s1v, s1data = self.iag.get_segment_edge(node=u, data=True)
                s1 = s1data["object"]
                s1id = s1.stable_id_non_hap
                s2u, s2v, s2data = self.iag.get_segment_edge(node=v, data=True)
                s2 = s2data["object"]
                s2id = s2.stable_id_non_hap
                s1_cn_boundaries = [self.scnb[clone_id].get_cnb(sid=s1id, hap=Haplotype.A, boundary_type=CNBoundaries.UPPER),
                                    self.scnb[clone_id].get_cnb(sid=s1id, hap=Haplotype.B, boundary_type=CNBoundaries.UPPER)]
                s2_cn_boundaries = [self.scnb[clone_id].get_cnb(sid=s2id, hap=Haplotype.A, boundary_type=CNBoundaries.UPPER),
                                    self.scnb[clone_id].get_cnb(sid=s2id, hap=Haplotype.B, boundary_type=CNBoundaries.UPPER)]
                c_max = max(s1_cn_boundaries + s2_cn_boundaries)
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    aid = adjacency.stable_id_non_phased
                    self.gm.addConstr(c_max * self.variables[P][clone_id][aid][ph], g.GRB.GREATER_EQUAL, self.variables[PROD][PY][clone_id][aid][ph])
                    self.gm.addConstr(self.variables[PROD][PY][clone_id][aid][ph], g.GRB.GREATER_EQUAL, 0.0)
                    if adjacency.adjacency_type == AdjacencyType.NOVEL:
                        a_cn_var = self.variables[YN][clone_id][aid]
                    else:
                        a_cn_var = self.variables[YR][clone_id][aid][ph]
                    self.gm.addConstr(a_cn_var, g.GRB.GREATER_EQUAL, self.variables[PROD][PY][clone_id][aid][ph])
                    self.gm.addConstr(a_cn_var - c_max * (1 - self.variables[P][clone_id][aid][ph]), g.GRB.LESS_EQUAL,
                                      self.variables[PROD][PY][clone_id][aid][ph])

    def define_constraints_for_labeling(self):
        self.define_constraints_for_ref_labeling()
        self.define_constraints_for_novel_labeling()
        self.define_constraints_on_each_location()

    def define_constraints_for_ref_labeling(self):
        """
        Using binary indicator variable p_{i,e} stored in `variables[P][clone_id][adjacency_id][phasing]` to ensure that reference adjacency (unlabeled version) has at most two phasing
            (AA and BB, as they are in the reference);

            phasing on Y chromosome (and X in case of males, respectively) are enforced with balancing on segments extremities and 0-0 boundaries on homologous copies
        """
        for clone_id in self.clone_ids:
            for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.REFERENCE, self.hapl_adjacencies):
                aid = adjacency.stable_id_non_phased
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables[P][clone_id][aid][ph], g.GRB.EQUAL, 0, name="ph_{{r,{cid},{aid},{phas}}}".format(cid=str(clone_id),
                                                                                                                                      aid=str(aid),
                                                                                                                                      phas=str(ph)))
                ###
                # If we are solving the problem in the haploid setting, we can not have BB reference adjacency present (i.e., CN must be equal to 0, we are achieving it
                #   via fixing binary indicator at 0)
                ###
                if self.solve_as_haploid:
                    self.gm.addConstr(self.variables[P][clone_id][aid][Phasing.BB], g.GRB.EQUAL, 0, name="ph_{{r,{cid},{aid},hapl-{phas}}}".format(cid=str(clone_id),
                                                                                                                                                   aid=str(aid),
                                                                                                                                                   phas=str(Phasing.BB)))

    def define_constraints_nov_adjacency_overall_presence(self):
        fp_lin_expr = g.LinExpr()
        hapl_nov_adjs_cnt = len(list(filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies)))
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            aid = adjacency.stable_id_non_phased
            fp_lin_expr.add(g.quicksum([self.variables[PROD][PP][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]), mult=(1.0 / hapl_nov_adjs_cnt))
        self.gm.addConstr(fp_lin_expr, g.GRB.GREATER_EQUAL, g.LinExpr(1 - self.hapl_nov_adjacencies_fp), name="nov-adj-fp")

    def define_constraints_for_novel_labeling(self):
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            aid = adjacency.stable_id_non_phased
            ###
            # tying clone specific phasing constraints to the overarching one
            #   Note: using Gurobi built-in OR constraint, rather than explicitly writing it down
            ###
            for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                self.gm.addGenConstrOr(self.variables[PROD][PP][aid][ph],
                                       [self.variables[P][clone_id][aid][ph] for clone_id in self.clone_ids],
                                       name="pp_{{{aid},{phas}}}".format(aid=str(aid), phas=str(ph)))
            ###
            # force presence of a diploid counterpart adjacency for an observed unlabeled one
            ###
            self.gm.addConstr(g.quicksum([self.variables[PROD][PP][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]),
                              g.GRB.LESS_EQUAL, 1, name="pp-nov-uniq-dip_{{{aid}}}".format(aid=str(aid)))

            ###
            # If the adjacency is a self-loop, then both AB, and BA phasing are prohibited
            #   Note, that the the constraint is achieved via a NON-clone specific binary variable, as the clone-specific ones are tied to the overarching one
            ###
            if adjacency.is_self_loop_hapl:
                for ph in [Phasing.AB, Phasing.BA]:
                    self.gm.addConstr(self.variables[PROD][PP][aid][ph], g.GRB.EQUAL, 0)
            ###
            # If we are solving the problem in a haploid setting, then we must force every labeling choice, except for AA one, to not be present
            # (i.e., indicator forced to equal to 0)
            ###
            if self.solve_as_haploid:
                for ph in [Phasing.AB, Phasing.BA, Phasing.BB]:
                    self.gm.addConstr(self.variables[PROD][PP][aid][ph], g.GRB.EQUAL, 0)

    def define_constraints_on_each_location(self):
        self.reciprocal_locations = []
        for u, v, data in self.iag.ref_adjacency_edges(data=True):
            u_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=u))
            v_na_edges_w_data = list(self.iag.nov_adjacency_edges(data=True, nbunch=v))
            if len(u_na_edges_w_data) == 0 or len(v_na_edges_w_data) == 0:
                continue
            self.reciprocal_locations.append((u, v))
            for u_na_edge_w_data in u_na_edges_w_data:
                for v_na_edge_w_data in v_na_edges_w_data:
                    u_na = u_na_edge_w_data[2]["object"]
                    v_na = v_na_edge_w_data[2]["object"]
                    u_na_id = u_na.stable_id_non_phased
                    v_na_id = v_na.stable_id_non_phased
                    expr1 = g.LinExpr()
                    expr2 = g.LinExpr()
                    expr1.add(self.variables[PROD][PP][u_na_id][get_aabb_for_ra(haplotype=Haplotype.A)], mult=1)
                    expr1.add(self.variables[PROD][PP][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.A)], mult=1)
                    expr1.add(self.variables[PROD][PP][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.B)], mult=1)
                    expr1.add(self.variables[PROD][PP][v_na_id][get_aabb_for_ra(haplotype=Haplotype.B)], mult=1)
                    self.gm.addConstr(expr1, g.GRB.LESS_EQUAL, 1, name="reciprocal_location_1_{{{u},{v},{u_aid},{v_aid}}}"
                                                                       "".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))
                    expr2.add(self.variables[PROD][PP][u_na_id][get_aabb_for_ra(haplotype=Haplotype.B)], mult=1)
                    expr2.add(self.variables[PROD][PP][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.B)], mult=1)
                    expr2.add(self.variables[PROD][PP][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.A)], mult=1)
                    expr2.add(self.variables[PROD][PP][v_na_id][get_aabb_for_ra(haplotype=Haplotype.A)], mult=1)
                    self.gm.addConstr(expr2, g.GRB.LESS_EQUAL, 1, name="reciprocal_location_2_{{{u},{v},{u_aid},{v_aid}}}"
                                                                       "".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))
            # if len(u_na_edges_w_data) > 1 or len(v_na_edges_w_data) > 1:
            #     continue
            # u_na = u_na_edges_w_data[0][2]["object"]
            # u_na_id = u_na.stable_id_non_phased
            # v_na = v_na_edges_w_data[0][2]["object"]
            # v_na_id = v_na.stable_id_non_phased
            # expr = g.LinExpr()
            # expr.add(g.LinExpr(self.variables[PROD][PP][u_na_id][get_aabb_for_ra(haplotype=Haplotype.A)]), mult=1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.A)]), mult=1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][v_na_id][get_aabb_for_ra(haplotype=Haplotype.A)]), mult=-1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.A)]), mult=-1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][v_na_id][get_aabb_for_ra(haplotype=Haplotype.B)]), mult=1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][v_na_id][get_abba_for_na_and_position(novel_adjacency=v_na, position=v, haplotype=Haplotype.B)]), mult=1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][u_na_id][get_aabb_for_ra(haplotype=Haplotype.B)]), mult=-1)
            # expr.add(g.LinExpr(self.variables[PROD][PP][u_na_id][get_abba_for_na_and_position(novel_adjacency=u_na, position=u, haplotype=Haplotype.B)]), mult=-1)
            # self.gm.addConstr(expr, g.GRB.GREATER_EQUAL, -1, name="recip_nov_l_{{{u},{v},{u_aid},{v_aid}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))
            # self.gm.addConstr(expr, g.GRB.LESS_EQUAL, 1, name="recip_nov_u_{{{u},{v},{u_aid},{v_aid}}}".format(u=str(u), v=str(v), u_aid=str(u_na_id), v_aid=str(v_na_id)))

    def define_constraints_for_adjacency_groups(self):
        self.define_constraints_for_adjacency_groups_molecule()
        self.define_constraints_for_adjacency_groups_general()
        self.define_constraints_for_adjacency_groups_labeling()

    def define_constraints_for_adjacency_groups_molecule(self):
        for adj_group in filter(lambda ag: ag.group_type == AdjacencyGroupType.MOLECULE, self.hapl_adjacencies_groups):
            fp = adj_group.extra.get(FALSE_POSITIVE, self.extra.get(DEFAULT_GROUP_M_FP, 0.1))
            for clone_id in self.clone_ids:
                group_var = self.variables[ADJ_GROUPS_M][clone_id][adj_group.gid]
                group_size = len(adj_group.adjacencies_ids)
                lin_expr = g.LinExpr()
                for adjacency in adj_group.adjacencies:
                    aid = adjacency.stable_id_non_phased
                    for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                        lin_expr.add(self.variables[P][clone_id][aid][ph])
                ###
                # in every clone, we force the total number of "realized" novel adjacencies in the group to be greater or equal to (1 - fp) * |u| * p_{i,u}
                #   when p_{i,u} is 0, this inequality always holds, while when it is 1, the total number is forced to be greater or equal to the (1 - fp) * |u| of the group
                ###
                self.gm.addConstr(g.LinExpr(group_var * group_size * (1 - fp)), g.GRB.LESS_EQUAL, lin_expr, name="group-molecule_{{{cid},{gid}}}".format(cid=str(clone_id),
                                                                                                                                                         gid=str(adj_group.gid)))
            ###
            # we force in at least 1 clone the p_{i, u} variable to be equal to 1, thus forcing at least one clone-specific inequality to force the group (fraction) realization
            ###
            self.gm.addConstr(g.quicksum([self.variables[ADJ_GROUPS_M][clone_id][adj_group.gid] for clone_id in self.clone_ids]),
                              g.GRB.GREATER_EQUAL, 1, name="group-molecule-across_{{{gid}}}".format(gid=str(adj_group.gid)))

    def define_constraints_for_adjacency_groups_general(self):
        for adj_group in filter(lambda ag: ag.group_type == AdjacencyGroupType.GENERAL, self.hapl_adjacencies_groups):
            fp = adj_group.extra.get(FALSE_POSITIVE, self.extra.get(DEFAULT_GROUP_N_FP, 0.1))
            group_size = len(adj_group.adjacencies_ids)
            fp_lin_expr = g.LinExpr()
            for adjacency in adj_group.adjacencies:
                aid = adjacency.stable_id_non_phased
                fp_lin_expr.add(g.quicksum([self.variables[PROD][PP][aid][ph] for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]]), mult=(1.0 / group_size))
            self.gm.addConstr(fp_lin_expr, g.GRB.GREATER_EQUAL, g.LinExpr(1 - fp), name="general_group-fp_{{{gid}}}".format(gid=adj_group.gid))

    def define_constraints_for_adjacency_groups_labeling(self):
        if self.solve_as_haploid:
            return
        self.hapl_adjacencies_by_external_ids = {adj.extra.get(EXTERNAL_NA_ID, adj.stable_id_non_phased): adj for adj in self.hapl_adjacencies}
        for adj_group in filter(lambda ag: ag.group_type == AdjacencyGroupType.LABELING, self.hapl_adjacencies_groups):
            aids = adj_group.adjacencies_ids
            indexes = adj_group.extra.get(AG_LABELING, [])
            if len(aids) != len(indexes):
                continue
            positions = []
            for aid, index in zip(aids, indexes):
                adjacency = self.hapl_adjacencies_by_external_ids[aid]
                position = adjacency.position1 if index == 0 else adjacency.position2
                positions.append(position)
            group_vars = []
            for haplotype in [Haplotype.A, Haplotype.B]:
                positions_adjacency_variables = []
                group_var = self.variables[ADJ_GROUPS_L][adj_group.gid][haplotype]
                group_vars.append(group_var)
                for position in positions:
                    for na_w_data in self.iag.nov_adjacency_edges(nbunch=position, data=True):
                        u, v, data = na_w_data
                        na = data["object"]
                        aid = na.stable_id_non_phased
                        positions_adjacency_variables.append(self.variables[PROD][PP][aid][get_aabb_for_ra(haplotype=haplotype)])
                        positions_adjacency_variables.append(self.variables[PROD][PP][aid][get_abba_for_na_and_position(novel_adjacency=na,
                                                                                                                        position=position, haplotype=haplotype)])
                self.gm.addGenConstrOr(group_var, positions_adjacency_variables, name="labeling_group-{{{gid},{hap}}}".format(gid=adj_group.gid,
                                                                                                                              hap=haplotype.value))
            self.gm.addConstr(g.quicksum(group_vars), g.GRB.LESS_EQUAL, 1, name="labeling_group-{{{gid}}}".format(gid=adj_group.gid))

    def define_constraints_on_nodes(self):
        """
        Balancing constraints with copy-number balanced enforced on nodes that are not specified as possible telomeres locations
            and copy-number excess allowed (though not enforced) on nodes that are specified as possible telomeres
        """
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
        s_var = self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][haplotype]
        result.add(s_var)
        ###
        # Here IAG is haploid, not diploid.
        ###
        ref_adj_edges_w_data = self.iag.ref_adjacency_edges(data=True, nbunch=node)
        ref_adjs = [data["object"] for _, __, data in ref_adj_edges_w_data]
        ###
        # Constraint on at most reference adjacency per node is
        #   subject to change is we were to identify segments as being present more than once in the reference. May require pan-graph version of the reference, if need be.
        ###
        assert len(ref_adjs) <= 1
        for adjacency in ref_adjs:
            result.add(g.LinExpr(self.variables[PROD][PY][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
        nov_adj_edges_w_data = self.iag.nov_adjacency_edges(data=True, nbunch=node)
        nov_adjs = [data["object"] for _, __, data in nov_adj_edges_w_data]
        for adjacency in nov_adjs:
            self_loop = adjacency.is_self_loop_hapl
            if self_loop:
                result.add(g.LinExpr(self.variables[PROD][PY][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-2)
            else:
                result.add(g.LinExpr(self.variables[PROD][PY][clone_id][adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)]), mult=-1)
                result.add(g.LinExpr(self.variables[PROD][PY][clone_id][adjacency.stable_id_non_phased][get_abba_for_na_and_position(novel_adjacency=adjacency,
                                                                                                                                     position=node,
                                                                                                                                     haplotype=haplotype)]), mult=-1)
        return result

    def define_constraints_on_deltas(self):
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                fid = self.hapl_segments_to_fragments[sid]
                f_var = self.variables[FRAGMENT_ALLELE][fid]
                s_cn_a_var = self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.A]
                s_cn_b_var = self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.B]
                ####
                # input allele-specific values
                ####
                s_cn_a = self.scnt[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.A)
                s_cn_b = self.scnt[clone_id].get_hap_aware_cn_by_seg_and_hap(segment=segment, haplotype=Haplotype.B)
                ####
                delta_var_a = self.variables[DELTA][clone_id][sid][Haplotype.A]
                delta_var_b = self.variables[DELTA][clone_id][sid][Haplotype.B]
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

    def define_objective(self):
        lin_exp = g.LinExpr()
        for clone_id in self.clone_ids:
            for segment in self.hapl_segments:
                sid = segment.stable_id_non_hap
                length = getattr(segment, self.extra.get(SEGMENT_LENGTH_ATTRIBUTE, "length_100"))
                delta_var_a = self.variables[DELTA][clone_id][sid][Haplotype.A]
                delta_var_b = self.variables[DELTA][clone_id][sid][Haplotype.B]
                lin_exp.add(g.LinExpr(delta_var_a + delta_var_b), mult=length)
        self.gm.setObjective(lin_exp, g.GRB.MINIMIZE)

    def get_scnt_from_model(self):
        # self.gm.update()
        result = {}
        for clone_id in self.clone_ids:
            scnp = SegmentCopyNumberProfile()
            result[clone_id] = scnp
            for s in self.hapl_segments:
                sid = s.stable_id_non_hap
                cn_a = int(round(self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.A].X))
                cn_b = int(round(self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][Haplotype.B].X))
                scnp.set_cn_record_for_segment(segment=s, cn=cn_a, haplotype=Haplotype.A)
                scnp.set_cn_record_for_segment(segment=s, cn=cn_b, haplotype=Haplotype.B)
        return result

    def get_acnt_from_model(self):
        # self.gm.update()
        result = {}
        for clone_id in self.clone_ids:
            acnp = AdjacencyCopyNumberProfile()
            result[clone_id] = acnp
            for adj in self.hapl_adjacencies:
                aid = adj.stable_id_non_phased
                for ph in [Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB]:
                    cn = int(round(self.variables[PROD][PY][clone_id][aid][ph].X))
                    acnp.set_cn_record_for_adjacency(adjacency=adj, cn=cn, phasing=ph)
        return result

    def alleles_sync_result(self, segment):
        sid = segment.stable_id_non_hap
        fid = self.hapl_segments_to_fragments[sid]
        return int(round(self.variables[FRAGMENT_ALLELE][fid].X))
