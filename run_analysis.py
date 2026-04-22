import cobra
import os
from metquest_clone.run import run as mq_run

def get_serum_context(sbml_path):
    """
    Automates seed and target identification using cobrapy.
    """
    model = cobra.io.read_sbml_model(sbml_path)
    
    # 1. SEEDS: Human Serum Environment (Barra et al. 2020)
    # We exclude B1/B6 here to force the bacteria to use its own pathways
    serum_bases = ["glc__D", "ala__L", "arg__L", "asn__L", "asp__L", 
                   "fe2", "pi", "h2o", "co2", "o2", "nh4"]
    
    seeds = []
    for s in serum_bases:
        try:
            seeds.append(model.metabolites.get_by_id(f"{s}_e").id)
        except:
            pass # Skip if specific ion/nutrient isn't in this strain
            
    # 2. TARGETS: The vitamins we want to check reachability for
    # B1 (Thiamine Diphosphate) and B6 (Pyridoxal 5'-phosphate)
    potential_targets = ["thmpp_c", "pdx5p_c"] 
    targets = [m.id for m in model.metabolites if m.id in potential_targets]
    
    return set(seeds), targets

def main():
    sbml_file = "models/iYS854.xml" # S. aureus USA300
    
    if not os.path.exists(sbml_file):
        print(f"File {sbml_file} not found!")
        return

    # Phase 1: Automated ID Discovery
    print(f"[Automation] Mapping IDs for {sbml_file}...")
    seeds, targets = get_serum_context(sbml_file)
    
    print(f"[Automation] Seeds Found: {len(seeds)}")
    print(f"[Automation] Targets Identified: {targets}")

    # Phase 2: Run your MetQuest Clone
    # We pass the sbml_file as a list as per your mq_run signature
    results = mq_run(
        sbml_paths=[sbml_file], 
        seed_set=seeds, 
        targets=targets, 
        beta=15, 
        output_csv="eskape_results.csv"
    )

    # Phase 3: Keystone Evaluation
    for t in targets:
        if t not in results['scope']:
            print(f"!!! CRITICAL: Target {t} is UNREACHABLE in Serum context.")
        else:
            print(f"Success: Target {t} is reachable via {len(results['pathways'][t])} paths.")

if __name__ == "__main__":
    main()