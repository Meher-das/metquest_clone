"""
metquest.parser
Parse one or more SBML files into a directed bipartite graph.

Single model
------------
Each SBML reaction becomes one or two nodes in the graph:
  <rxn_id>_fwd   substrates → reaction → products
  <rxn_id>_rev   products   → reaction → substrates  (reversible only)

Community (multiple SBML files)
--------------------------------
Individual models are linked through a shared extracellular space.
Exchange reactions (id starts with "EX_", or one side is empty) expose
metabolites to that shared space. Two models that both exchange the same
species are connected via that species node — exactly as described in
Ravikrishnan et al. (2018).

Exchange reaction detection
----------------------------
A reaction is treated as an exchange if:
  • its id starts with "EX_"  (BiGG / COBRA convention), OR
  • one of its reactant / product lists is empty (boundary flux)

Compartment namespacing
-----------------------
When multiple models are loaded, each species id is prefixed with the
model id to avoid collisions (e.g. "iJO1366::glc__D_c").
Exchange species are NOT prefixed — they are the shared currency between
models and must map to the same node across files.
"""

import os
from xml.etree import ElementTree as ET


# ── XML helpers ───────────────────────────────────────────────────────────────

def _ns(root):
    """Extract XML namespace from root tag, e.g. '{http://...}'."""
    return root.tag.split('}')[0] + '}' if '}' in root.tag else ''


def _find_list(parent, ns, list_tag, item_tag):
    lst = parent.find(f"{ns}{list_tag}")
    return lst.findall(f"{ns}{item_tag}") if lst is not None else []


# ── exchange detection ────────────────────────────────────────────────────────

def _is_exchange(rxn_id, subs, pros):
    return rxn_id.startswith("EX_") or not subs or not pros


# ── single-file parser ────────────────────────────────────────────────────────

def _parse_sbml(path, graph, meta, model_prefix=None):
    """
    Parse one SBML file into graph / meta in-place.

    model_prefix : str | None
        When set, non-exchange species ids are prefixed as
        "<model_prefix>::<species_id>" to namespace multi-model graphs.
    """
    root = ET.parse(path).getroot()
    ns = _ns(root)
    model = root.find(f"{ns}model")
    if model is None:
        raise ValueError(f"No <model> element found in {path}")

    # Collect exchange species so we know which ids stay unprefixed
    exchange_species = set()
    for rxn in _find_list(model, ns, 'listOfReactions', 'reaction'):
        rid = rxn.get('id', '')
        subs = [sr.get('species','') for sr in _find_list(rxn, ns, 'listOfReactants', 'speciesReference')]
        pros = [sr.get('species','') for sr in _find_list(rxn, ns, 'listOfProducts',  'speciesReference')]
        if _is_exchange(rid, subs, pros):
            exchange_species.update(subs)
            exchange_species.update(pros)

    def _sid(species_id):
        """Apply model prefix to non-exchange species."""
        if model_prefix and species_id not in exchange_species:
            return f"{model_prefix}::{species_id}"
        return species_id

    # Build graph from reactions
    for rxn in _find_list(model, ns, 'listOfReactions', 'reaction'):
        rid_raw = rxn.get('id', '')
        rev     = rxn.get('reversible', 'false').lower() == 'true'
        subs    = [_sid(sr.get('species','')) for sr in _find_list(rxn, ns, 'listOfReactants', 'speciesReference') if sr.get('species')]
        pros    = [_sid(sr.get('species','')) for sr in _find_list(rxn, ns, 'listOfProducts',  'speciesReference') if sr.get('species')]

        # Skip reactions with both sides empty
        if not subs and not pros:
            continue

        rid = f"{model_prefix}::{rid_raw}" if model_prefix else rid_raw

        for direction, ins, outs in [("fwd", subs, pros)] + ([("rev", pros, subs)] if rev else []):
            node = f"{rid}_{direction}"
            graph.setdefault(node, set())
            for m in ins:
                graph.setdefault(m, set()).add(node)
            for m in outs:
                graph.setdefault(m, set())
                graph[node].add(m)
            meta[node] = {
                "substrates": ins,
                "products":   outs,
                "original_id": rid_raw,
                "direction":   direction,
                "model":       model_prefix or "",
                "exchange":    _is_exchange(rid_raw, subs, pros),
            }


# ── public API ────────────────────────────────────────────────────────────────

def parse_sbml(paths):
    """
    Parse one or more SBML files into a directed bipartite graph.

    Parameters
    ----------
    paths : str | list[str]
        Path(s) to SBML file(s) or a directory of .sbml / .xml files.
        A single file → no model prefixing.
        Multiple files → each model's non-exchange species are prefixed
        with the model's SBML id to prevent collisions.

    Returns
    -------
    graph : dict[str, set[str]]
    meta  : dict[str, dict]
        Keys per entry: substrates, products, original_id,
                        direction, model, exchange.
    """
    if isinstance(paths, str):
        paths = [paths]

    # Expand directories
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += [os.path.join(p, f) for f in sorted(os.listdir(p))
                      if f.endswith((".sbml", ".xml"))]
        else:
            files.append(p)

    if not files:
        raise ValueError(f"No SBML files found in: {paths}")

    graph, meta = {}, {}
    multi = len(files) > 1

    for fp in files:
        # Derive model prefix from SBML <model id="..."> when combining
        prefix = None
        if multi:
            root = ET.parse(fp).getroot()
            ns   = _ns(root)
            model_el = root.find(f"{ns}model")
            prefix = (model_el.get('id') or os.path.splitext(os.path.basename(fp))[0]) if model_el is not None \
                     else os.path.splitext(os.path.basename(fp))[0]

        _parse_sbml(fp, graph, meta, model_prefix=prefix)
    return graph, meta