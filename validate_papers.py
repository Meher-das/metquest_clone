import os
import sys
from run import run

# =============================================================================
# Helper functions for formatting and output
# =============================================================================

def format_reaction(rxn_id):
    """
    Cleans up the parser's directional suffixes for terminal output.
    Converts 'rxn_id_fwd' -> 'rxn_id' and 'rxn_id_rev' -> 'rxn_id (reverse)'
    """
    if rxn_id.endswith('_fwd'):
        return rxn_id[:-4]
    elif rxn_id.endswith('_rev'):
        return rxn_id[:-4] + " (reverse)"
    return rxn_id

def print_validation_summary(paper_name, target, results):
    print(f"\n[{paper_name}] Validation Results:")
    pathways = results.get("pathways", {})
    target_pathways = pathways.get(target, [])
    
    if target_pathways:
        print(f"  SUCCESS: Found {len(target_pathways)} pathway(s) to {target}.")
        # Print the shortest pathway found
        shortest = min(target_pathways, key=len)
        formatted_rxns = [format_reaction(r) for r in shortest]
        
        print(f"  Shortest pathway length: {len(shortest)} reactions.")
        print(f"  Shortest pathway reactions: {', '.join(formatted_rxns)}")
    else:
        print(f"  FAILED: Could not find a pathway to {target} within the beta limit.")
        if target not in results.get("scope", set()):
            print(f"  Reason: {target} never entered the reachable scope during Phase 1 BFS.")
    print("-" * 60)


# =============================================================================
# Validation Setups for Each Paper
# =============================================================================

def validate_1_4_bdo():
    """
    Paper 1: Direct production of 1,4-butanediol (BDO) in E. coli.
    Pathway: Succinate -> Succinyl-CoA -> Succinate semialdehyde -> 
             4-hydroxybutyrate -> 4-hydroxybutyryl-CoA -> 
             4-hydroxybutyraldehyde -> 1,4-butanediol.
    """
    print("\nStarting Validation for Paper 1: 1,4-Butanediol...")
    
    sbml_files = ["models/e_coli_bdo_engineered.xml"] 
    if not all(os.path.exists(f) for f in sbml_files):
        print(f"  [SKIPPED] Model missing: {sbml_files[0]}")
        return

    # Seed set: Central carbon metabolites (assuming glucose feed)
    seeds = {"M_glc__D_c", "M_succ_c", "M_coa_c", "M_nadh_c", "M_atp_c", "M_h_c"} 
    target = "M_14bdo_c" 
    
    # Beta = 15. The engineered pathway from succinate is ~6 steps.
    results = run(sbml_paths=sbml_files, seed_set=seeds, targets=[target], beta=15, verbose=False)
    print_validation_summary("Paper 1: 1,4-BDO", target, results)


def validate_1_3_pdo():
    """
    Paper 2: 1,3-propanediol (1,3-PDO) in E. coli.
    Pathway: Glucose -> DHAP -> Glycerol-3-P -> Glycerol -> 
             3-hydroxypropionaldehyde (3-HPA) -> 1,3-propanediol.
    """
    print("\nStarting Validation for Paper 2: 1,3-Propanediol...")
    
    sbml_files = ["models/e_coli_pdo_engineered.xml"]
    if not all(os.path.exists(f) for f in sbml_files):
        print(f"  [SKIPPED] Model missing: {sbml_files[0]}")
        return

    # Seed set: Glucose + required cofactors 
    seeds = {"M_glc__D_c", "M_nadh_c", "M_h2o_c", "M_cobalt2_c", "M_atp_c"} 
    target = "M_13pdo_c" 
    
    # Beta = 10. The pathway from glucose via DHAP is relatively short.
    results = run(sbml_paths=sbml_files, seed_set=seeds, targets=[target], beta=10, verbose=False)
    print_validation_summary("Paper 2: 1,3-PDO", target, results)


def validate_artemisinin():
    """
    Paper 3: Artemisinic Acid in S. cerevisiae (Yeast).
    Pathway: Acetyl-CoA -> Mevalonate pathway -> FPP -> Amorphadiene ->
             Artemisinic alcohol -> Artemisinic aldehyde -> Artemisinic acid.
    """
    print("\nStarting Validation for Paper 3: Artemisinic Acid...")
    
    sbml_files = ["models/s_cerevisiae_artemisinic_engineered.xml"]
    if not all(os.path.exists(f) for f in sbml_files):
        print(f"  [SKIPPED] Model missing: {sbml_files[0]}")
        return

    # Seed set: Glucose, Acetyl-CoA, and essential cofactors for cytochromes
    seeds = {"M_glc__D_c", "M_accoa_c", "M_nadph_c", "M_o2_c", "M_nad_c", "M_atp_c"}
    target = "M_artemisinic_acid_c"
    
    # Beta = 25. High beta needed for long eukaryotic pathways and secondary metabolism.
    results = run(sbml_paths=sbml_files, seed_set=seeds, targets=[target], beta=25, verbose=False)
    print_validation_summary("Paper 3: Artemisinic Acid", target, results)


if __name__ == "__main__":
    print("=====================================================")
    print("   MetQuest Paper Validation Suite   ")
    print("=====================================================")
    
    if not os.path.exists("models"):
        print("\n[WARNING] 'models/' directory not found.")
        print("Please create a 'models' folder and add your engineered SBML files:\n"
              "  - models/e_coli_bdo_engineered.xml\n"
              "  - models/e_coli_pdo_engineered.xml\n"
              "  - models/s_cerevisiae_artemisinic_engineered.xml\n")
    else:
        validate_1_4_bdo()
        validate_1_3_pdo()
        validate_artemisinin()