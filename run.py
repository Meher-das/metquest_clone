"""
metquest.run
High-level entry point.

    from metquest import run
    results = run(sbml_paths, seed_set, targets, beta=15, output_csv=None)

CLI:
    python -m metquest.run --sbml DIR --seeds cpd:C00031 --targets cpd:C00041 --beta 15
"""
import argparse, time
from parser              import parse_sbml
from guided_bfs          import guided_bfs
from dynamic_programming import enumerate_pathways
from output              import get_pathway_list, write_csv

try:
    from .output import get_nx_subgraphs
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


def run(sbml_paths, seed_set, targets, beta=15, output_csv=None, verbose=True):
    if isinstance(targets, str):
        targets = [targets]
    t0 = time.time()

    def log(msg):
        if verbose: print(msg)

    log("[MetQuest] Parsing SBML…")
    graph, meta = parse_sbml(sbml_paths)
    log(f"           {len(meta)} reaction nodes, {len(graph) - len(meta)} metabolite nodes")

    log("[MetQuest] Phase 1: Guided BFS…")
    scope, Rv, level_m, level_r = guided_bfs(graph, meta, seed_set)
    log(f"           Scope: {len(scope)}  |  Visited reactions: {len(Rv)}")

    missing = [t for t in targets if t not in scope]
    if missing:
        log(f"  [WARNING] Unreachable targets: {missing}")

    log(f"[MetQuest] Phase 2: DP enumeration (β={beta})…")
    table = enumerate_pathways(scope, Rv, seed_set, meta, beta)

    pathways = get_pathway_list(table, targets, beta)
    for t, pws in pathways.items():
        log(f"           {t}: {len(pws)} pathway(s)")

    subgraphs = {}
    if _HAS_NX:
        subgraphs = get_nx_subgraphs(table, targets, beta, graph, meta)

    if output_csv:
        write_csv(table, targets, beta, meta, list(seed_set), output_csv)
        log(f"[MetQuest] CSV → {output_csv}")

    log(f"[MetQuest] Done in {time.time()-t0:.2f}s")
    return {"graph": graph, "meta": meta, "scope": scope, "Rv": Rv,
            "level_m": level_m, "level_r": level_r,
            "table": table, "pathways": pathways, "subgraphs": subgraphs}


def _cli():
    p = argparse.ArgumentParser(prog="metquest")
    p.add_argument("--sbml",    required=True, nargs="+")
    p.add_argument("--seeds",   required=True, nargs="+")
    p.add_argument("--targets", required=True, nargs="+")
    p.add_argument("--beta",    type=int, default=15)
    p.add_argument("--csv",     default=None)
    args = p.parse_args()
    run(args.sbml, set(args.seeds), args.targets, args.beta, args.csv)


if __name__ == "__main__":
    _cli()
