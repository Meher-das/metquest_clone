import pandas as pd
import cobra
from run import run as mq_run

def validate_and_test_orthogonality(mq_results, host_model_path, literature_pathway=None):
    """
    mq_results: The dictionary returned by your run() function.
    host_model_path: Path to the SBML file of your host (e.g., iML1515.xml).
    literature_pathway: A list of reaction IDs from a paper to compare against.
    """
    pathways = mq_results["pathways"]
    all_metabolites = mq_results["graph"]
    host_model = cobra.io.read_sbml_model(host_model_path)
    
    report = []

    for target, pws in pathways.items():
        for i, pw in enumerate(pws):
            # 1. ORTHOGONALITY: Competition Index
            # We check how many 'seed' metabolites this path drains from essential biomass precursors
            precursors = [node for node in pw if node in mq_results["level_m"] and mq_results["level_m"][node] == 0]
            
            # Simple Orthogonality Check: Does it use high-demand precursors?
            # E.g., Acetyl-CoA (cpd00022), ATP (cpd00002), NADPH (cpd00005)
            high_demand = {'cpd00022', 'cpd00002', 'cpd00005', 'cpd00010'} 
            overlap = set(precursors).intersection(high_demand)
            ortho_score = "High" if not overlap else "Medium/Low (Competition Risk)"

            # 2. VALIDATION: Literature Matching
            match_score = 0
            if literature_pathway:
                # Calculate Jaccard Similarity between MetQuest path and Paper path
                set_mq = set(pw)
                set_lit = set(literature_pathway)
                intersection = set_mq.intersection(set_lit)
                match_score = (len(intersection) / len(set_lit)) * 100

            report.append({
                "Target": target,
                "Path_ID": i,
                "Length": len(pw),
                "Orthogonality": ortho_score,
                "Shared_Precursors": list(overlap),
                "Lit_Match_%": round(match_score, 2)
            })

    return pd.DataFrame(report)

# --- EXAMPLE USAGE ---
# 1. Run your existing MetQuest code
results = mq_run(["models/iYS854.xml"], {"cpd00031"}, ["cpd00041"])

# 2. Define a pathway from a paper (e.g., Yim et al. 1,4-BDO reactions)
paper_rxns = ["RXN_1", "RXN_2", "RXN_3", "RXN_4"]

# 3. Generate the Validation Report
df_report = validate_and_test_orthogonality(results, "iML1515.xml", literature_pathway=paper_rxns)
print(df_report)