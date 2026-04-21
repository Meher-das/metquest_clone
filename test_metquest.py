"""
Tests using synthetic SBML written to temp files — no external data needed.

Single-model toy network
------------------------
  HEX1     (irrev): A + ATP -> B + ADP
  R2       (irrev): B -> C
  R3       (irrev): A + C -> D
  R4       (rev):   D <-> E
  R5       (irrev): B + C -> F
  EX_A     (exchange, rev): A <->      [boundary]

Seed = {A, ATP}. Expected scope includes B, C, D, E, F.
Pathway to D (size 3): {HEX1, R2, R3}
Pathway to F (size 3): {HEX1, R2, R5}

Community test
--------------
Two models sharing metabolite X via exchange reactions.
Model 1 produces X; model 2 consumes X to make Z.
"""
import os, sys, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from metquest.parser              import parse_sbml
from metquest.guided_bfs          import guided_bfs
from metquest.dynamic_programming import enumerate_pathways
from metquest.output              import write_csv


# ── SBML builders ─────────────────────────────────────────────────────────────

def _sbml(model_id, species, reactions):
    """
    Minimal SBML Level 2v4 string.

    species   : list of species ids
    reactions : list of (id, reversible, [substrates], [products])
    """
    def rxn_xml(rid, rev, subs, pros):
        def srefs(ids, tag):
            if not ids: return ""
            refs = "".join(f'<speciesReference species="{s}"/>' for s in ids)
            return f"<{tag}>{refs}</{tag}>"
        return (f'<reaction id="{rid}" reversible="{"true" if rev else "false"}">'
                f'{srefs(subs,"listOfReactants")}{srefs(pros,"listOfProducts")}'
                f'</reaction>')

    sp = "".join(f'<species id="{s}" compartment="c"/>' for s in species)
    rx = "".join(rxn_xml(*r) for r in reactions)
    return (f'<?xml version="1.0"?>'
            f'<sbml xmlns="http://www.sbml.org/sbml/level2/version4" level="2" version="4">'
            f'<model id="{model_id}">'
            f'<listOfCompartments><compartment id="c"/></listOfCompartments>'
            f'<listOfSpecies>{sp}</listOfSpecies>'
            f'<listOfReactions>{rx}</listOfReactions>'
            f'</model></sbml>')


TOY_SPECIES = ["A","B","C","D","E","F","ATP","ADP"]
TOY_REACTIONS = [
    ("HEX1", False, ["A","ATP"], ["B","ADP"]),
    ("R2",   False, ["B"],       ["C"]),
    ("R3",   False, ["A","C"],   ["D"]),
    ("R4",   True,  ["D"],       ["E"]),
    ("R5",   False, ["B","C"],   ["F"]),
    ("EX_A", True,  ["A"],       []),
]


def _write_sbml(content):
    """Write SBML string to a temp file, return path."""
    f = tempfile.NamedTemporaryFile(suffix=".sbml", mode='w', delete=False)
    f.write(content)
    f.close()
    return f.name


def _run(seed, beta=4, paths=None):
    if paths is None:
        p = _write_sbml(_sbml("toy", TOY_SPECIES, TOY_REACTIONS))
        paths = p
    graph, meta = parse_sbml(paths)
    scope, Rv, lm, lr = guided_bfs(graph, meta, seed)
    table = enumerate_pathways(scope, Rv, seed, meta, beta)
    return scope, Rv, table, meta


def _pws(table, target, beta=4):
    return [pw for k in range(1, beta+1)
            if table[target][k] is not None
            for pw in table[target][k]]


# ── parser tests ──────────────────────────────────────────────────────────────

def test_parse_single():
    p = _write_sbml(_sbml("toy", TOY_SPECIES, TOY_REACTIONS))
    graph, meta = parse_sbml(p)
    os.unlink(p)
    assert "HEX1_fwd" in meta
    assert "R4_fwd"   in meta
    assert "R4_rev"   in meta          # reversible → both directions
    assert "EX_A_fwd" in meta
    assert meta["EX_A_fwd"]["exchange"] is True
    assert meta["HEX1_fwd"]["exchange"] is False


def test_parse_reversible_swap():
    p = _write_sbml(_sbml("toy", TOY_SPECIES, TOY_REACTIONS))
    graph, meta = parse_sbml(p)
    os.unlink(p)
    assert meta["R4_fwd"]["substrates"] == ["D"]
    assert meta["R4_fwd"]["products"]   == ["E"]
    assert meta["R4_rev"]["substrates"] == ["E"]   # swapped
    assert meta["R4_rev"]["products"]   == ["D"]


def test_parse_exchange_no_products():
    """EX_A has empty products — detected as exchange regardless of prefix."""
    p = _write_sbml(_sbml("toy", TOY_SPECIES, TOY_REACTIONS))
    graph, meta = parse_sbml(p)
    os.unlink(p)
    assert meta["EX_A_fwd"]["exchange"] is True


# ── BFS tests ─────────────────────────────────────────────────────────────────

def test_scope():
    scope, *_ = _run({"A", "ATP"})
    assert {"B","C","D","E","F"} <= scope


def test_Rv():
    _, Rv, *_ = _run({"A", "ATP"})
    assert {"HEX1_fwd","R2_fwd","R3_fwd","R4_fwd","R5_fwd"} <= Rv


def test_unreachable_without_ATP():
    scope, *_ = _run({"A"})   # ATP missing → HEX1 stuck → nothing downstream
    assert "B" not in scope


def test_reversible_bfs():
    scope, Rv, *_ = _run({"E", "ATP"})
    assert "D" in scope and "R4_rev" in Rv


# ── DP pathway tests ──────────────────────────────────────────────────────────

def test_pathway_to_D():
    _, _, table, meta = _run({"A","ATP"})
    pws = _pws(table, "D")
    orig_ids = {frozenset(meta[r]["original_id"] for r in pw) for pw in pws}
    assert frozenset(["HEX1","R2","R3"]) in orig_ids


def test_pathway_to_F():
    _, _, table, meta = _run({"A","ATP"})
    pws = _pws(table, "F")
    orig_ids = {frozenset(meta[r]["original_id"] for r in pw) for pw in pws}
    assert frozenset(["HEX1","R2","R5"]) in orig_ids


# ── community test ────────────────────────────────────────────────────────────

def test_community_shared_exchange():
    """
    Model 1: P -> X  (internal), EX_X: X <-> (exchange)
    Model 2: EX_X: X <-> (exchange), X -> Z  (internal)
    Seed = {P}. Z should be reachable via the shared X exchange node.
    """
    m1 = _sbml("org1",
               ["P","X"],
               [("R_PtoX", False, ["P"], ["X"]),
                ("EX_X",   True,  ["X"], [])])
    m2 = _sbml("org2",
               ["X","Z"],
               [("EX_X",   True,  ["X"], []),
                ("R_XtoZ", False, ["X"], ["Z"])])

    p1 = _write_sbml(m1)
    p2 = _write_sbml(m2)
    try:
        # In multi-model graphs, non-exchange species are prefixed with model id.
        # Seed must use the prefixed id "org1::P", not bare "P".
        scope, Rv, table, meta = _run({"org1::P"}, beta=6, paths=[p1, p2])
        assert "org2::Z" in scope, "Z must be reachable via shared exchange metabolite X"
    finally:
        os.unlink(p1); os.unlink(p2)


def test_community_model_prefixing():
    """Non-exchange species must be namespaced; exchange species must not."""
    m1 = _sbml("org1", ["A","X"], [("R1", False, ["A"], ["X"]), ("EX_X", True, ["X"], [])])
    m2 = _sbml("org2", ["X","B"], [("EX_X", True, ["X"], []), ("R2", False, ["X"], ["B"])])
    p1, p2 = _write_sbml(m1), _write_sbml(m2)
    try:
        graph, meta = parse_sbml([p1, p2])
        assert "X" in graph                   # shared exchange species: no prefix
        assert "org1::A" in graph             # internal species: prefixed
        assert "org2::B" in graph
    finally:
        os.unlink(p1); os.unlink(p2)


# ── output test ───────────────────────────────────────────────────────────────

def test_csv_output():
    import csv as csvmod
    _, _, table, meta = _run({"A","ATP"})
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "out.csv")
        write_csv(table, ["D","F"], 4, meta, ["A","ATP"], path)
        rows = list(csvmod.DictReader(open(path)))
    assert len(rows) >= 2
    assert all(k in rows[0] for k in
               ["pathway_id","source","target","size","reaction_list","metabolite_list"])


if __name__ == "__main__":
    tests = [test_parse_single, test_parse_reversible_swap, test_parse_exchange_no_products,
             test_scope, test_Rv, test_unreachable_without_ATP, test_reversible_bfs,
             test_pathway_to_D, test_pathway_to_F,
             test_community_shared_exchange, test_community_model_prefixing,
             test_csv_output]
    for fn in tests:
        fn()
        print(f"PASS  {fn.__name__}")
    print(f"\nAll {len(tests)} tests passed.")
