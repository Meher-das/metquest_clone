import os
from xml.etree import ElementTree as ET
import time

def parse_sbml(files):

    # start_time = time.time()

    if isinstance(files, str):
        files = [files]

    graph, meta = {}, {}
    
    for path in files:
        # 1. Get the XML 'address' (namespace) and the Model ID
        root = ET.parse(path).getroot()
        ns = root.tag.split('}')[0] + '}' if '}' in root.tag else ''
        model_el = root.find(f"{ns}model")
        # Use the ID in the file, or the filename if ID is missing
        m_id = model_el.get('id') or os.path.basename(path).split('.')[0]

        # 2. Extract reactions
        rxn_list = model_el.find(f"{ns}listOfReactions")
        for rxn in rxn_list.findall(f"{ns}reaction"):
            r_raw = rxn.get('id')
            is_rev = rxn.get('reversible') == 'true'
            
            # Helper to pull metabolite IDs and apply prefix
            def get_m(tag):
                lst = rxn.find(f"{ns}{tag}")
                if lst is None: return []
                # Prefixing EVERY species since you said models are distinct
                return [f"{m_id}::{sr.get('species')}" for sr in lst.findall(f"{ns}speciesReference")]

            subs, pros = get_m('listOfReactants'), get_m('listOfProducts')
            
            # 3. Build the Bipartite Graph
            # Determine directions (always forward, plus reverse if applicable)
            directions = [("fwd", subs, pros)] + ([("rev", pros, subs)] if is_rev else [])
            
            for suffix, inputs, outputs in directions:
                node = f"{m_id}::{r_raw}_{suffix}"
                graph[node] = set(outputs) # Reaction -> Products
                for m in inputs:
                    graph.setdefault(m, set()).add(node) # Reactant -> Reaction
                
                meta[node] = {"model": m_id, "subs": inputs, "pros": outputs}
    # end_time = time.time()
    # print(f"Parsing time: {end_time - start_time:.2f} seconds")
    return graph, meta

# parse_sbml("/home/mandrin/Documents/MetQuest_Scratch/e_coli_core.xml")