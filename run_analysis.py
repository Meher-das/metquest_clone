import cobra
from .run import run as mq_run

def main():
    # 1. SETUP FILENAME AND MODEL ID
    # Your parser uses the filename or the internal SBML ID for prefixing
    sbml_path = "models/iYS854.xml"
    model_id = "iYS854" 

    # 2. SEED GENERATION (The "Shared Currency")
    # Your parser treats exchange reactions as unprefixed.
    # So we provide the IDs exactly as they appear in the <listOfReactions>
    # for the exchange reactions.
    seeds = {
        "glc__D_e", "ala__L_e", "arg__L_e", "asn__L_e", 
        "fe2_e", "pi_e", "h2o_e", "nh4_e"
    }

    # 3. TARGET GENERATION (The "Prefixed Internals")
    # Because these are cytoplasmic (_c), your parser will prefix them.
    # Logic: f"{model_id}::{internal_id}"
    targets = [
        f"{model_id}::thmpp_c", # Vitamin B1
        f"{model_id}::pdx5p_c"  # Vitamin B6
    ]

    # 4. EXECUTION
    # This calls your clone's high-level entry point
    results = mq_run(
        sbml_paths=[sbml_path], 
        seed_set=seeds, 
        targets=targets, 
        beta=15
    )

    # 5. KEYSTONE IDENTIFICATION (The "What if" analysis)
    # To identify a keystone, you would theoretically remove a reaction 
    # from the 'graph' returned in results and re-run guided_bfs.
    print(f"Reachable nodes in this context: {len(results['scope'])}")