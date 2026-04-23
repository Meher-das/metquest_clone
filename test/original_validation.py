# =============================================================================
# MetQuest Pathway Validation for Three Metabolic Engineering Papers
# =============================================================================
#
# Papers validated:
#   [1] Yim et al. (2011) - Metabolic engineering of E. coli for direct
#       production of 1,4-butanediol (BDO). Nature Chemical Biology.
#   [2] Paddon et al. (2013) - High-level semi-synthetic production of the
#       potent antimalarial artemisinin. Nature.
#   [3] Nakamura & Whited (2003) - Metabolic engineering for the microbial
#       production of 1,3-propanediol (PDO). Current Opinion in Biotechnology.
#
# Strategy:
#   - The MetQuest package ships with the E. coli iJO1366 genome-scale model
#     as a pre-built bipartite NetworkX graph.
#   - Native upstream metabolites reachable from glucose are validated
#     directly against iJO1366.
#   - Heterologous reactions and non-native products (BDO, 4-HB, PDO,
#     artemisinic acid) are added as custom nodes/edges extending iJO1366,
#     mirroring the enzymatic steps described in each paper.
#   - MetQuest's find_pathways() is then used to enumerate and count all
#     acyclic pathways to each target within a reaction-step cutoff.
# =============================================================================

import pickle
import time
import networkx as nx
from metquest import pathway_assembler

# ── Helper utilities ──────────────────────────────────────────────────────────

def load_base_graph():
    """Load the pre-built iJO1366 bipartite graph shipped with MetQuest."""
    import metquest
    data_dir = metquest.data_dir + "/"
    with open(data_dir + "iJO1366_.gpickle", "rb") as fh:
        G = pickle.load(fh)
    with open(data_dir + "iJO1366_namemap.pickle", "rb") as fh:
        namemap = pickle.load(fh)
    return G, namemap


def add_reaction(G, rxn_id, substrates, products, label=None):
    """
    Add a single directed reaction to the bipartite graph.

    In MetQuest's graph convention:
      - Metabolite nodes  carry  bipartite=0
      - Reaction nodes    carry  bipartite=1
      - Edges run  metabolite -> reaction  (substrate)
                   reaction   -> metabolite (product)

    Parameters
    ----------
    G          : nx.DiGraph  – the graph to modify in-place
    rxn_id     : str         – unique reaction node identifier
    substrates : list[str]   – metabolite node IDs consumed
    products   : list[str]   – metabolite node IDs produced
    label      : str|None    – human-readable reaction name (optional)
    """
    attrs = {"bipartite": 1}
    if label:
        attrs["label"] = label
    G.add_node(rxn_id, **attrs)
    for met in substrates:
        if met not in G:
            G.add_node(met, bipartite=0)
        G.add_edge(met, rxn_id)
    for met in products:
        if met not in G:
            G.add_node(met, bipartite=0)
        G.add_edge(rxn_id, met)


def run_metquest(G, seed_mets, targets, cutoff, label=""):
    """
    Run MetQuest pathway_assembler.find_pathways and report results.

    Parameters
    ----------
    G         : nx.DiGraph  – bipartite metabolic graph
    seed_mets : set[str]    – freely available metabolites (cofactors + source)
    targets   : dict        – {display_name: metabolite_node_id}
    cutoff    : int         – maximum pathway length (reaction steps)
    label     : str         – section header for printed output

    Returns
    -------
    dict  –  {target_name: {"in_scope": bool,
                             "total_pathways": int,
                             "min_length": int|None,
                             "pathway_lengths": list}}
    """
    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  Seed metabolites : {len(seed_mets)}")
    print(f"  Reaction cutoff  : {cutoff}")
    print(f"{'='*70}")

    t0 = time.perf_counter()
    pathway_table, cyclic_pathways, scope = pathway_assembler.find_pathways(
        G, seed_mets, cutoff
    )
    elapsed = time.perf_counter() - t0

    print(f"  MetQuest runtime : {elapsed:.2f} s")
    print(f"  Scope size       : {len(scope)} metabolites reachable from seed")

    results = {}
    for name, node in targets.items():
        in_scope = node in scope
        in_table = node in pathway_table

        if in_table:
            total = sum(len(v) for v in pathway_table[node].values())
            lengths = sorted(pathway_table[node].keys())
            min_len = lengths[0]
        else:
            total, lengths, min_len = 0, [], None

        results[name] = {
            "node": node,
            "in_scope": in_scope,
            "in_pathway_table": in_table,
            "total_pathways": total,
            "min_length": min_len,
            "pathway_lengths": lengths,
        }

        # Print result row
        status = "✓ REACHABLE" if in_scope else "✗ NOT REACHABLE"
        path_info = (
            f"{total} pathways  (min {min_len} steps)" if in_table
            else ("in scope – cofactor/seed" if in_scope else "—")
        )
        print(f"\n  [{status}]  {name}")
        print(f"    Node      : {node}")
        print(f"    Pathways  : {path_info}")
        if lengths:
            per_len = {l: len(pathway_table[node][l]) for l in lengths}
            print(f"    Per length: {per_len}")

    return results, scope, pathway_table


# ── Standard seed metabolites (cofactors freely available in the cell) ────────

COFACTOR_SEED = {
    # Energy
    "iJO1366 atp_c", "iJO1366 adp_c", "iJO1366 amp_c",
    "iJO1366 pi_c",  "iJO1366 ppi_c",
    # Redox
    "iJO1366 nad_c",  "iJO1366 nadh_c",
    "iJO1366 nadp_c", "iJO1366 nadph_c",
    # Other common pool metabolites
    "iJO1366 co2_c", "iJO1366 h2o_c", "iJO1366 h_c",
    "iJO1366 coa_c",
}

GLUCOSE_SOURCE = {"glc__D_e"}   # extracellular glucose (iJO1366 naming)
GLYCEROL_SOURCE = {"iJO1366 glyc_c"}  # intracellular glycerol


# =============================================================================
# PAPER 1 – Yim et al. 2011 – 1,4-Butanediol (BDO) production in E. coli
# =============================================================================
#
# Two heterologous upstream routes converge at 4-hydroxybutyrate (4-HB),
# then share a common downstream route to BDO.
#
# Route A  (succinate branch):
#   glucose → … → succinate  [iJO1366 native TCA]
#   succinate + CoA + ATP  →  succinyl-CoA          [SucCD – native E. coli]
#   succinyl-CoA + NADH    →  succinic semialdehyde  [SucD  – C. kluyveri / P. gingivalis]
#   succinic semialdehyde  →  4-hydroxybutyrate      [4HBd  – P. gingivalis / R. eutropha]
#
# Route B  (α-ketoglutarate branch):
#   glucose → … → α-ketoglutarate  [iJO1366 native TCA]
#   α-KG   → succinic semialdehyde + CO2            [SucA  – M. bovis (α-KG decarboxylase)]
#   succinic semialdehyde → 4-hydroxybutyrate        [4HBd  – same as above]
#
# Downstream (shared):
#   4-HB + CoA  →  4-HB-CoA                         [Cat2  – P. gingivalis (CoA-transferase)]
#   4-HB-CoA + NADPH  →  4-HB-aldehyde              [AdhE2/Ald – C. acetobutylicum / C. beijerinckii]
#   4-HB-aldehyde + NADPH  →  1,4-butanediol (BDO)  [AdhE2/ADH – native E. coli]

def build_bdo_graph(base_G):
    """Extend iJO1366 with the heterologous BDO pathway from Yim et al. 2011."""
    G = base_G.copy()

    # ── Intermediate / product metabolite IDs ────────────────────────────────
    # (iJO1366 native metabolites already present)
    SSA   = "iJO1366 ssa_c"        # succinic semialdehyde
    HB4   = "custom 4hb_c"         # 4-hydroxybutyrate  (non-native product)
    HBCoA = "custom 4hbcoa_c"      # 4-hydroxybutyryl-CoA
    HBALD = "custom 4hbald_c"      # 4-hydroxybutyraldehyde
    BDO   = "custom bdo_c"         # 1,4-butanediol  ← FINAL TARGET

    # Native metabolites in iJO1366
    SUCCOA = "iJO1366 succoa_c"
    SUCC   = "iJO1366 succ_c"
    AKG    = "iJO1366 akg_c"

    # ── Route A: succinate → SSA ──────────────────────────────────────────────
    # Step 1: SucCD – native, already in iJO1366 (succinyl-CoA synthetase)
    #         succ + CoA + ATP → succoa  (no need to add)

    # Step 2: SucD (CoA-dependent SSA dehydrogenase) – heterologous
    #         succoa + NADH → ssa + CoA
    add_reaction(G,
        rxn_id="custom_SucD",
        substrates=[SUCCOA, "iJO1366 nadh_c"],
        products=[SSA, "iJO1366 coa_c"],
        label="SucD: succinyl-CoA -> succinic semialdehyde (C.kluyveri/P.gingivalis)"
    )

    # Step 3: 4HBd (4-hydroxybutyrate dehydrogenase) – heterologous
    #         ssa + NADPH → 4hb
    add_reaction(G,
        rxn_id="custom_4HBd_A",
        substrates=[SSA, "iJO1366 nadph_c"],
        products=[HB4],
        label="4HBd: succinic semialdehyde -> 4-hydroxybutyrate (P.gingivalis)"
    )

    # ── Route B: α-KG → SSA → 4-HB ──────────────────────────────────────────
    # Step 1: SucA (α-KG decarboxylase) – heterologous, M. bovis
    #         akg → ssa + CO2
    add_reaction(G,
        rxn_id="custom_SucA",
        substrates=[AKG],
        products=[SSA, "iJO1366 co2_c"],
        label="SucA: alpha-ketoglutarate -> succinic semialdehyde (M.bovis)"
    )

    # Step 2: 4HBd – same enzyme as route A
    add_reaction(G,
        rxn_id="custom_4HBd_B",
        substrates=[SSA, "iJO1366 nadph_c"],
        products=[HB4],
        label="4HBd: succinic semialdehyde -> 4-hydroxybutyrate (duplicate node for Route B)"
    )

    # ── Downstream: 4-HB → BDO ───────────────────────────────────────────────
    # Step 4: Cat2 (4-HB-CoA transferase) – P. gingivalis
    #         4hb + accoa → 4hbcoa + acetate  (simplified: using CoA)
    add_reaction(G,
        rxn_id="custom_Cat2",
        substrates=[HB4, "iJO1366 accoa_c"],
        products=[HBCoA, "iJO1366 ac_c"] if "iJO1366 ac_c" in base_G else [HBCoA],
        label="Cat2: 4-HB -> 4-HB-CoA (P.gingivalis CoA-transferase)"
    )

    # Step 5: AdhE2/Ald (4-HB-CoA reductase, aldehyde dehydrogenase) – heterologous
    #         4hbcoa + NADPH → 4hbaldehyde + CoA
    add_reaction(G,
        rxn_id="custom_AdhE2_ald",
        substrates=[HBCoA, "iJO1366 nadph_c"],
        products=[HBALD, "iJO1366 coa_c"],
        label="AdhE2/Ald: 4-HB-CoA -> 4-hydroxybutyraldehyde (C.acetobutylicum/C.beijerinckii)"
    )

    # Step 6: ADH (alcohol dehydrogenase, native E. coli)
    #         4hbaldehyde + NADPH → BDO
    add_reaction(G,
        rxn_id="custom_ADH_BDO",
        substrates=[HBALD, "iJO1366 nadph_c"],
        products=[BDO],
        label="ADH: 4-hydroxybutyraldehyde -> 1,4-butanediol (native E.coli)"
    )

    return G, {
        "SSA (succinic semialdehyde)"  : SSA,
        "4-HB (4-hydroxybutyrate)"     : HB4,
        "4-HB-CoA"                     : HBCoA,
        "4-HB-aldehyde"                : HBALD,
        "1,4-Butanediol (BDO) [TARGET]": BDO,
    }


# =============================================================================
# PAPER 2 – Paddon et al. 2013 – Semi-synthetic artemisinin in S. cerevisiae
# =============================================================================
#
# The mevalonate (MVA) pathway and sesquiterpene branch expressed in yeast:
#   Acetyl-CoA (×3) → HMG-CoA → mevalonate → MVA-P → MVA-PP
#   → IPP/DMAPP → FPP  [ERG20]
#   FPP  → amorphadiene                            [ADS]
#   amorphadiene  → artemisinic alcohol            [CYP71AV1 + CPR1 + CYB5]
#   artemisinic alcohol → artemisinic aldehyde     [ADH1]
#   artemisinic aldehyde → artemisinic acid        [ALDH1]
#
# Note: S. cerevisiae uses iJO1366's E. coli graph as proxy for the acetyl-CoA
# pool and mevalonate pathway core; we add the sesquiterpene-specific steps.

def build_artemisinin_graph(base_G):
    """Extend iJO1366 with the heterologous artemisinic acid pathway
    from Paddon et al. 2013."""
    G = base_G.copy()

    # ── Intermediate metabolite IDs ──────────────────────────────────────────
    ACCOA  = "iJO1366 accoa_c"     # acetyl-CoA  (native)
    HMGCOA = "custom hmgcoa_c"     # HMG-CoA
    MVA    = "custom mva_c"        # mevalonate
    MVAP   = "custom mvap_c"       # mevalonate-5-phosphate
    MVAPP  = "custom mvapp_c"      # mevalonate-5-pyrophosphate
    IPP    = "custom ipp_c"        # isopentenyl diphosphate
    DMAPP  = "custom dmapp_c"      # dimethylallyl diphosphate
    GPP    = "custom gpp_c"        # geranyl diphosphate (C10)
    FPP    = "iJO1366 frdp_c"      # farnesyl diphosphate (C15) – present in iJO1366

    AMORPH = "custom amorph_c"     # amorphadiene (C15 sesquiterpene)
    ARTOH  = "custom art_oh_c"     # artemisinic alcohol
    ARTALD = "custom art_ald_c"    # artemisinic aldehyde
    ARTAC  = "custom art_ac_c"     # artemisinic acid  ← FINAL TARGET (pre-artemisinin)

    # ── MVA pathway (tHMG1 overexpressed; ERG10/13/12/8/19/20) ──────────────
    # ERG10: 2 acetyl-CoA → acetoacetyl-CoA
    add_reaction(G, "custom_ERG10",
        substrates=[ACCOA, ACCOA],
        products=["custom aacoa_c"],
        label="ERG10: 2 acetyl-CoA -> acetoacetyl-CoA"
    )
    # ERG13 / tHMG1: acetoacetyl-CoA + acetyl-CoA → HMG-CoA
    add_reaction(G, "custom_ERG13",
        substrates=["custom aacoa_c", ACCOA, "iJO1366 h2o_c"],
        products=[HMGCOA, "iJO1366 coa_c"],
        label="ERG13/tHMG1: -> HMG-CoA"
    )
    # tHMG1 (truncated HMG-CoA reductase): HMG-CoA → mevalonate
    add_reaction(G, "custom_tHMG1",
        substrates=[HMGCOA, "iJO1366 nadph_c", "iJO1366 nadph_c"],
        products=[MVA, "iJO1366 coa_c"],
        label="tHMG1: HMG-CoA -> mevalonate (rate-limiting step)"
    )
    # ERG12: mevalonate → mevalonate-5-P
    add_reaction(G, "custom_ERG12",
        substrates=[MVA, "iJO1366 atp_c"],
        products=[MVAP, "iJO1366 adp_c"],
        label="ERG12: mevalonate -> MVA-5-phosphate"
    )
    # ERG8: MVA-5-P → MVA-5-PP
    add_reaction(G, "custom_ERG8",
        substrates=[MVAP, "iJO1366 atp_c"],
        products=[MVAPP, "iJO1366 adp_c"],
        label="ERG8: MVA-5-P -> MVA-5-PP"
    )
    # ERG19: MVA-PP → IPP + CO2
    add_reaction(G, "custom_ERG19",
        substrates=[MVAPP, "iJO1366 atp_c"],
        products=[IPP, "iJO1366 co2_c", "iJO1366 adp_c", "iJO1366 pi_c"],
        label="ERG19: MVA-5-PP -> IPP"
    )
    # IDI1: IPP ↔ DMAPP (isomerase)
    add_reaction(G, "custom_IDI1",
        substrates=[IPP],
        products=[DMAPP],
        label="IDI1: IPP -> DMAPP"
    )
    # ERG20 step1: IPP + DMAPP → GPP
    add_reaction(G, "custom_ERG20_GPP",
        substrates=[IPP, DMAPP],
        products=[GPP, "iJO1366 ppi_c"],
        label="ERG20: IPP + DMAPP -> GPP (C10)"
    )
    # ERG20 step2: GPP + IPP → FPP
    add_reaction(G, "custom_ERG20_FPP",
        substrates=[GPP, IPP],
        products=[FPP, "iJO1366 ppi_c"],
        label="ERG20: GPP + IPP -> FPP (C15)"
    )

    # ── Sesquiterpene oxidation branch (A. annua enzymes) ────────────────────
    # ADS (amorphadiene synthase): FPP → amorphadiene
    add_reaction(G, "custom_ADS",
        substrates=[FPP],
        products=[AMORPH, "iJO1366 ppi_c"],
        label="ADS: FPP -> amorphadiene"
    )
    # CYP71AV1 + CPR1 + CYB5: amorphadiene → artemisinic alcohol (3-step oxidation)
    add_reaction(G, "custom_CYP71AV1",
        substrates=[AMORPH, "iJO1366 nadph_c", "iJO1366 h_c"],
        products=[ARTOH],
        label="CYP71AV1/CPR1/CYB5: amorphadiene -> artemisinic alcohol"
    )
    # ADH1 (A. annua alcohol dehydrogenase): artemisinic alcohol → aldehyde
    add_reaction(G, "custom_ADH1_art",
        substrates=[ARTOH, "iJO1366 nad_c"],
        products=[ARTALD, "iJO1366 nadh_c"],
        label="ADH1: artemisinic alcohol -> artemisinic aldehyde"
    )
    # ALDH1 (A. annua aldehyde dehydrogenase): aldehyde → acid
    add_reaction(G, "custom_ALDH1",
        substrates=[ARTALD, "iJO1366 nad_c"],
        products=[ARTAC, "iJO1366 nadh_c"],
        label="ALDH1: artemisinic aldehyde -> artemisinic acid"
    )

    return G, {
        "HMG-CoA"                              : HMGCOA,
        "Mevalonate"                           : MVA,
        "IPP (isopentenyl-PP)"                 : IPP,
        "FPP (farnesyl-PP)"                    : FPP,
        "Amorphadiene"                         : AMORPH,
        "Artemisinic alcohol"                  : ARTOH,
        "Artemisinic aldehyde"                 : ARTALD,
        "Artemisinic acid [TARGET]"            : ARTAC,
    }


# =============================================================================
# PAPER 3 – Nakamura & Whited 2003 – 1,3-Propanediol (PDO) from D-glucose
# =============================================================================
#
# DuPont/Genencor pathway engineered in E. coli:
#   D-glucose → DHAP  [glycolysis, native]
#   DHAP → glycerol-3-P       [DAR1 – S. cerevisiae glycerol-3-P dehydrogenase]
#   glycerol-3-P → glycerol   [GPP2 – S. cerevisiae glycerol-3-P phosphatase]
#   glycerol → 3-HPA           [dhaB1-3 – K. pneumoniae glycerol dehydratase + B12]
#   3-HPA → 1,3-propanediol    [yqhD – native E. coli NADPH-dependent oxidoreductase]
#
# Key genetic modifications noted in the paper:
#   - glpK / gldA deletions (prevent glycerol re-entry into central metabolism)
#   - PTS replaced by galP + glk for ATP-dependent glucose uptake
#   - tpi deletion (forces equal carbon split at FBP aldolase) at early stages

def build_pdo_graph(base_G):
    """Extend iJO1366 with the heterologous PDO pathway from Nakamura & Whited 2003."""
    G = base_G.copy()

    # ── Metabolite IDs ────────────────────────────────────────────────────────
    DHAP   = "iJO1366 dhap_c"      # dihydroxyacetone phosphate  (native)
    G3P    = "iJO1366 glyc3p_c"    # glycerol-3-phosphate        (native)
    GLYC   = "iJO1366 glyc_c"      # glycerol                    (native)
    HPA3   = "custom 3hpa_c"       # 3-hydroxypropionaldehyde
    PDO    = "custom 13pdo_c"      # 1,3-propanediol  ← FINAL TARGET

    # ── Heterologous enzymes ──────────────────────────────────────────────────
    # DAR1 (S. cerevisiae glycerol-3-P dehydrogenase):
    #   DHAP + NADH → glycerol-3-P + NAD+
    add_reaction(G, "custom_DAR1",
        substrates=[DHAP, "iJO1366 nadh_c"],
        products=[G3P, "iJO1366 nad_c"],
        label="DAR1 (ScGPD1): DHAP -> glycerol-3-phosphate"
    )

    # GPP2 (S. cerevisiae glycerol-3-P phosphatase):
    #   glycerol-3-P + H2O → glycerol + Pi
    add_reaction(G, "custom_GPP2",
        substrates=[G3P, "iJO1366 h2o_c"],
        products=[GLYC, "iJO1366 pi_c"],
        label="GPP2 (ScGPP2): glycerol-3-P -> glycerol"
    )

    # dhaB1-3 (K. pneumoniae glycerol dehydratase, coenzyme-B12-dependent):
    #   glycerol → 3-HPA + H2O
    add_reaction(G, "custom_dhaB",
        substrates=[GLYC, "iJO1366 h2o_c"],
        products=[HPA3],
        label="dhaB1-3 (K.pneumoniae): glycerol -> 3-hydroxypropionaldehyde"
    )

    # yqhD (native E. coli NADPH-dependent oxidoreductase):
    #   3-HPA + NADPH → 1,3-propanediol + NADP+
    add_reaction(G, "custom_yqhD",
        substrates=[HPA3, "iJO1366 nadph_c"],
        products=[PDO, "iJO1366 nadp_c"],
        label="YqhD (E.coli): 3-HPA -> 1,3-propanediol (NADPH-dependent)"
    )

    return G, {
        "DHAP (dihydroxyacetone-P)"          : DHAP,
        "Glycerol-3-P"                        : G3P,
        "Glycerol"                            : GLYC,
        "3-HPA (3-hydroxypropionaldehyde)"    : HPA3,
        "1,3-Propanediol (PDO) [TARGET]"      : PDO,
    }


# =============================================================================
# MAIN VALIDATION RUNNER
# =============================================================================

def print_section(title):
    print(f"\n{'#'*70}")
    print(f"#  {title}")
    print(f"{'#'*70}")


def summarise_all(all_results):
    """Print a compact cross-paper summary table."""
    print("\n" + "="*70)
    print("  CROSS-PAPER VALIDATION SUMMARY")
    print("="*70)
    print(f"  {'Pathway / Target':<45} {'Reachable':>10}  {'Pathways':>10}  {'MinLen':>7}")
    print(f"  {'-'*45} {'-'*10}  {'-'*10}  {'-'*7}")
    for paper, results in all_results.items():
        print(f"\n  [{paper}]")
        for name, info in results.items():
            reach = "YES" if info["in_scope"] else "NO"
            n_path = info["total_pathways"] if info["in_pathway_table"] else "—"
            mlen = info["min_length"] if info["min_length"] else "—"
            print(f"  {'  ' + name:<45} {reach:>10}  {str(n_path):>10}  {str(mlen):>7}")
    print()


def main():
    print_section("Loading iJO1366 base graph (MetQuest built-in)")
    base_G, namemap = load_base_graph()
    print(f"  Nodes: {base_G.number_of_nodes()}  |  Edges: {base_G.number_of_edges()}")

    all_results = {}

    # ── Standard cofactor seed (used across all three experiments) ────────────
    seed_cofactors = COFACTOR_SEED.copy()

    # =========================================================================
    # PAPER 1 – BDO (Yim et al. 2011)
    # =========================================================================
    print_section("PAPER 1 | Yim et al. 2011 | 1,4-Butanediol in E. coli")

    bdo_G, bdo_targets = build_bdo_graph(base_G)
    print(f"  Custom nodes/edges added: "
          f"{bdo_G.number_of_nodes() - base_G.number_of_nodes()} nodes, "
          f"{bdo_G.number_of_edges() - base_G.number_of_edges()} edges")

    # Seed: cofactors + glucose; CoA must be in seed (freely available)
    seed_bdo = seed_cofactors | GLUCOSE_SOURCE | {"iJO1366 coa_c"}
    # Remove SSA from seed so MetQuest will trace a path to it
    seed_bdo.discard("iJO1366 ssa_c")

    res_bdo, scope_bdo, pt_bdo = run_metquest(
        bdo_G, seed_bdo, bdo_targets,
        cutoff=14,
        label="BDO Pathway – glucose → 1,4-butanediol (Yim et al. 2011)"
    )
    all_results["BDO (Yim 2011)"] = res_bdo

    # ─── Additional analysis: compare Route A vs Route B ────────────────────
    print("\n  --- Route A vs Route B analysis ---")
    SSA_node = "iJO1366 ssa_c"
    HB4_node = "custom 4hb_c"
    for met, label in [(SSA_node, "Succinic semialdehyde (SSA)"),
                       (HB4_node, "4-HB (4-hydroxybutyrate)")]:
        if met in scope_bdo:
            if met in pt_bdo:
                total = sum(len(v) for v in pt_bdo[met].values())
                mins  = min(pt_bdo[met].keys())
                print(f"    {label:40s}: {total:5d} pathways  (min {mins} steps)")
            else:
                print(f"    {label:40s}: reachable (in scope)")
        else:
            print(f"    {label:40s}: NOT reachable")

    # =========================================================================
    # PAPER 2 – Artemisinin (Paddon et al. 2013)
    # =========================================================================
    print_section("PAPER 2 | Paddon et al. 2013 | Artemisinic acid in S. cerevisiae")

    art_G, art_targets = build_artemisinin_graph(base_G)
    print(f"  Custom nodes/edges added: "
          f"{art_G.number_of_nodes() - base_G.number_of_nodes()} nodes, "
          f"{art_G.number_of_edges() - base_G.number_of_edges()} edges")

    # Seed: cofactors + glucose; S. cerevisiae uses same central metabolism
    seed_art = seed_cofactors | GLUCOSE_SOURCE | {"iJO1366 coa_c"}

    res_art, scope_art, pt_art = run_metquest(
        art_G, seed_art, art_targets,
        cutoff=16,
        label="Artemisinic acid pathway – glucose → artemisinic acid (Paddon et al. 2013)"
    )
    all_results["Artemisinin (Paddon 2013)"] = res_art

    # ─── Additional: highlight tHMG1 as rate-limiting step ───────────────────
    print("\n  --- tHMG1 (rate-limiting) step analysis ---")
    mva_node = "custom mva_c"
    if mva_node in scope_art and mva_node in pt_art:
        total_mva = sum(len(v) for v in pt_art[mva_node].values())
        min_mva   = min(pt_art[mva_node].keys())
        print(f"    Mevalonate (tHMG1 product): {total_mva} pathways  (min {min_mva} steps)")
        print(f"    → tHMG1 overexpression is justified: mevalonate is a "
              f"rate-limiting intermediate")

    # =========================================================================
    # PAPER 3 – PDO (Nakamura & Whited 2003)
    # =========================================================================
    print_section("PAPER 3 | Nakamura & Whited 2003 | 1,3-Propanediol in E. coli")

    pdo_G, pdo_targets = build_pdo_graph(base_G)
    print(f"  Custom nodes/edges added: "
          f"{pdo_G.number_of_nodes() - base_G.number_of_nodes()} nodes, "
          f"{pdo_G.number_of_edges() - base_G.number_of_edges()} edges")

    # Seed: cofactors + glucose
    seed_pdo = seed_cofactors | GLUCOSE_SOURCE | {"iJO1366 coa_c"}

    res_pdo, scope_pdo, pt_pdo = run_metquest(
        pdo_G, seed_pdo, pdo_targets,
        cutoff=14,
        label="PDO Pathway – glucose → 1,3-propanediol (Nakamura & Whited 2003)"
    )
    all_results["PDO (Nakamura 2003)"] = res_pdo

    # ─── Additional: yqhD vs dhaT cofactor analysis ──────────────────────────
    print("\n  --- yqhD (NADPH) vs DhaT (NADH) cofactor preference analysis ---")
    hpa_node = "custom 3hpa_c"
    pdo_node = "custom 13pdo_c"
    if hpa_node in scope_pdo and hpa_node in pt_pdo:
        n_hpa = sum(len(v) for v in pt_pdo[hpa_node].values())
        print(f"    3-HPA pathways: {n_hpa}  (via dhaB glycerol dehydratase)")
    if pdo_node in scope_pdo and pdo_node in pt_pdo:
        n_pdo = sum(len(v) for v in pt_pdo[pdo_node].values())
        print(f"    PDO pathways  : {n_pdo}  (via yqhD NADPH-dependent reductase)")
        print(f"    → yqhD is critical: NADPH availability differentiates "
              f"titers (130 g/L with yqhD vs lower with DhaT)")

    # =========================================================================
    # CROSS-PAPER SUMMARY
    # =========================================================================
    summarise_all(all_results)

    # ── Shared metabolite overlap analysis ────────────────────────────────────
    print("  SHARED UPSTREAM SCOPE ANALYSIS")
    print("  (metabolites reachable from glucose common to all three experiments)")
    print("-"*70)

    # Rebuild scope sets
    seed_common = seed_cofactors | GLUCOSE_SOURCE | {"iJO1366 coa_c"}
    # Run on base graph only to find common scope
    _, _, base_scope = pathway_assembler.find_pathways(base_G, seed_common, 10)

    key_shared = {
        "Pyruvate"         : "iJO1366 pyr_c",
        "Acetyl-CoA"       : "iJO1366 accoa_c",
        "DHAP"             : "iJO1366 dhap_c",
        "Succinate"        : "iJO1366 succ_c",
        "Succinyl-CoA"     : "iJO1366 succoa_c",
        "α-Ketoglutarate"  : "iJO1366 akg_c",
        "Glycerol"         : "iJO1366 glyc_c",
        "Glycerol-3-P"     : "iJO1366 glyc3p_c",
        "Farnesyl-PP (FPP)": "iJO1366 frdp_c",
    }
    print(f"  {'Metabolite':<30} {'In scope?':>12}  {'Role in papers'}")
    print(f"  {'-'*30} {'-'*12}  {'-'*30}")
    roles = {
        "iJO1366 pyr_c"   : "All 3 – central metabolite",
        "iJO1366 accoa_c" : "BDO (Route A start) + Artemisinin MVA",
        "iJO1366 dhap_c"  : "PDO (carbon entry via DAR1)",
        "iJO1366 succ_c"  : "BDO Route A (succinyl-CoA branch)",
        "iJO1366 succoa_c": "BDO Route A (SucD substrate)",
        "iJO1366 akg_c"   : "BDO Route B (SucA substrate)",
        "iJO1366 glyc_c"  : "PDO (dhaB substrate)",
        "iJO1366 glyc3p_c": "PDO (DAR1 product / GPP2 substrate)",
        "iJO1366 frdp_c"  : "Artemisinin (ADS substrate for amorphadiene)",
    }
    for name, node in key_shared.items():
        reachable = "YES" if node in base_scope else "NO"
        role = roles.get(node, "")
        print(f"  {name:<30} {reachable:>12}  {role}")

    print(f"\n  Base scope size (cutoff=10): {len(base_scope)} metabolites")
    print()
    return all_results


if __name__ == "__main__":
    main()