"""Known-answer tests of the amendment-A4 split (recipe.nucfit_refold.split_group). Hand-made graphs; the expected
removals follow from the rule: remove from the first violating component, tiers Rfam seed + MDSA, then Rfam large,
then BRAliBASE; inside a tier most links inside the component, ties by case id."""
from recipe.nucfit_refold import components, split_group


def clique(names):
    return [(a, b, "exact") for i, a in enumerate(names) for b in names[i + 1:]]


A = [f"bA{i}" for i in range(4)]
B = [f"bB{i}" for i in range(4)]


def test_one_hub_between_two_families():
    # two BRAliBASE cliques of 4 joined only through the Rfam seed case h; limit 5 -> remove h only
    links = clique(A) + clique(B) + [("h", x, "identity") for x in A + B]
    bench = {x: "bralibase" for x in A + B} | {"h": "rfam"}
    removed, comps = split_group(A + B + ["h"], links, bench, 5, {"bralibase": 99, "rfam": 99})
    assert removed == ["h"]
    assert comps == [sorted(A), sorted(B)]


def test_tier_order_seed_before_large():
    # two bridges: Rfam large L (degree 8) and Rfam seed s (degree 2: one link into each family).
    # Tier 1 (seed) goes first although L has more links; s alone does not split (L still bridges), so L follows.
    links = clique(A) + clique(B) + [("L", x, "identity") for x in A + B] + [("s", "bA0", "exact"), ("s", "bB0", "exact")]
    bench = {x: "bralibase" for x in A + B} | {"L": "rfam_large", "s": "rfam"}
    removed, comps = split_group(A + B + ["L", "s"], links, bench, 5, {"bralibase": 99, "rfam": 99, "rfam_large": 99})
    assert removed == ["s", "L"]
    assert comps == [sorted(A), sorted(B)]


def test_bralibase_only_as_last_tier():
    # a BRAliBASE chain x0-x1-x2-x3-x4 (no other benchmark); limit 2: highest degree first (x1, x2, x3 have 2;
    # tie -> case id x1), then the rest {x0}, {x2, x3, x4} -> x3 (degree 2) -> pieces of <= 2
    xs = [f"x{i}" for i in range(5)]
    links = [(xs[i], xs[i + 1], "exact") for i in range(4)]
    bench = {x: "bralibase" for x in xs}
    removed, comps = split_group(xs, links, bench, 2, {"bralibase": 99})
    assert removed == ["x1", "x3"]
    assert comps == [["x0"], ["x2"], ["x4"]]


def test_benchmark_share_limit():
    # size limit not hit (9 <= 10) but one component holds 8 BRAliBASE cases > per-benchmark limit 5 -> remove hub
    links = clique(A) + clique(B) + [("h", x, "identity") for x in A + B]
    bench = {x: "bralibase" for x in A + B} | {"h": "mdsa"}
    removed, comps = split_group(A + B + ["h"], links, bench, 10, {"bralibase": 5, "mdsa": 99})
    assert removed == ["h"] and comps == [sorted(A), sorted(B)]


def test_nothing_to_do():
    links = clique(A)
    removed, comps = split_group(A, links, {x: "bralibase" for x in A}, 10, {"bralibase": 10})
    assert removed == [] and comps == [sorted(A)]


def test_components_deterministic_order():
    adj = {"a": {"b"}, "b": {"a"}, "c": set(), "d": {"e"}, "e": {"d"}}
    assert components({"a", "b", "c", "d", "e"}, adj) == [["a", "b"], ["d", "e"], ["c"]]
