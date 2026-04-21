"""
metquest.guided_bfs
Phase 1: Guided BFS from seed set S.

Fires reactions whose every substrate is in scope, expands scope with their
products, and re-triggers stuck reactions until no new metabolites appear.

Returns: scope, Rv, level_m, level_r
"""
from collections import deque


def guided_bfs(graph, reaction_meta, seed_set):
    scope = set(seed_set)
    Rv, stuck = set(), set()
    level_m = {m: 0 for m in seed_set}
    level_r = {}

    def try_fire(rid):
        subs = reaction_meta[rid]["substrates"]
        if all(s in scope for s in subs):
            if rid in Rv:
                return []
            Rv.add(rid)
            stuck.discard(rid)
            level_r[rid] = max((level_m.get(s, 0) for s in subs), default=0) + 1
            new = []
            for p in reaction_meta[rid]["products"]:
                if p not in scope:
                    scope.add(p)
                    level_m[p] = level_r[rid] + 1
                    new.append(p)
            return new
        stuck.add(rid)
        return []

    queue = deque()
    for rid in reaction_meta:
        for m in try_fire(rid):
            queue.append(m)

    while queue:
        met = queue.popleft()
        for rid in list(stuck):
            for m in try_fire(rid):
                queue.append(m)
        for nb in graph.get(met, set()):
            if nb in reaction_meta and nb not in Rv:
                for m in try_fire(nb):
                    queue.append(m)

    return scope, Rv, level_m, level_r
