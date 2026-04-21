"""
metquest.output
Three output formats: frozenset list, CSV, NetworkX DiGraph subgraphs.
"""
import csv, os

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


def _pathways(table, target, beta):
    if target not in table:
        return []
    return [pw for k in range(1, beta + 1)
            if k < len(table[target]) and table[target][k] is not None
            for pw in table[target][k]]


def _mets(pw, meta):
    mets = set()
    for rid in pw:
        mets.update(meta.get(rid, {}).get("substrates", []))
        mets.update(meta.get(rid, {}).get("products", []))
    return sorted(mets)


def get_pathway_list(table, targets, beta):
    if isinstance(targets, str):
        targets = [targets]
    return {t: _pathways(table, t, beta) for t in targets}


def get_nx_subgraphs(table, targets, beta, graph, meta):
    if not _HAS_NX:
        raise ImportError("networkx required: pip install networkx")
    if isinstance(targets, str):
        targets = [targets]
    result = {}
    for t in targets:
        graphs = []
        for pw in _pathways(table, t, beta):
            G = nx.DiGraph()
            for rid in pw:
                G.add_node(rid, node_type="reaction")
                for s in meta.get(rid, {}).get("substrates", []):
                    G.add_node(s, node_type="metabolite")
                    G.add_edge(s, rid)
                for p in meta.get(rid, {}).get("products", []):
                    G.add_node(p, node_type="metabolite")
                    G.add_edge(rid, p)
            graphs.append(G)
        result[t] = graphs
    return result


def write_csv(table, targets, beta, meta, sources, path):
    if isinstance(targets, str):
        targets = [targets]
    if isinstance(sources, str):
        sources = [sources]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["pathway_id", "source", "target", "size", "reaction_list", "metabolite_list"])
        src = ";".join(sorted(sources))
        pid = 1
        for t in targets:
            for pw in _pathways(table, t, beta):
                w.writerow([pid, src, t, len(pw), ";".join(sorted(pw)), ";".join(_mets(pw, meta))])
                pid += 1
    return path
