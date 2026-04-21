"""
metquest.dynamic_programming
Phase 2: DP enumeration of all pathways (Algorithms 1 & 2, Ravikrishnan 2018).

Table[m][k] = set of frozensets, each a pathway of size k producing m from S.
None (⊥) means m cannot be produced in exactly k reactions.
"""
from itertools import product as iproduct


def _partitions(total, n, max_val):
    """All n-tuples of non-negative integers summing to total, each ≤ max_val."""
    if n == 1:
        if 0 <= total <= max_val:
            yield (total,)
        return
    for first in range(min(total, max_val) + 1):
        for rest in _partitions(total - first, n - 1, max_val):
            yield (first,) + rest


def _populate(rid, partition, ins, outs, table):
    if any(table[ins[i]][partition[i]] is None for i in range(len(ins))):
        return
    for combo in iproduct(*[table[ins[i]][partition[i]] for i in range(len(ins))]):
        pw = frozenset().union(*combo) | {rid}
        j = len(pw)
        for y in outs:
            if j < len(table[y]):
                if table[y][j] is None:
                    table[y][j] = set()
                table[y][j].add(pw)


def enumerate_pathways(scope, Rv, seed_set, reaction_meta, beta):
    table = {m: ([{frozenset()}] + [None] * beta if m in seed_set else [None] * (beta + 1))
             for m in scope}

    for rid in Rv:
        ins, outs = reaction_meta[rid]["substrates"], reaction_meta[rid]["products"]
        if all(s in seed_set for s in ins):
            pw = frozenset([rid])
            for y in outs:
                if y in table:
                    if table[y][1] is None:
                        table[y][1] = set()
                    table[y][1].add(pw)

    for k in range(2, beta + 1):
        for rid in Rv:
            ins, outs = reaction_meta[rid]["substrates"], reaction_meta[rid]["products"]
            if not any(y in table for y in outs):
                continue
            n = len(ins)
            if not n:
                continue
            for ell in range(k - 1, n * (k - 1) + 1):
                for p in _partitions(ell, n, k - 1):
                    _populate(rid, p, ins, outs, table)

    return table
