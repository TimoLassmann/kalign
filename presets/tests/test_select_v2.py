"""Amendment v2 mode rule (recipe.select_v2)."""

import math

from recipe import select_v2 as S


def pt(key, cpu, sp, tc=None):
    return {"key": key, "cpu": cpu, "sp": sp, "tc": sp if tc is None else tc, "source": "stage1"}


def test_accurate_stops_where_doubling_gains_less_than_the_threshold():
    front = [pt("a", 1, 0.70), pt("b", 2, 0.71), pt("c", 4, 0.715), pt("d", 8, 0.716)]
    picks, _ = S.knee_modes(front, fastest_other_cpu=100)
    # c -> d gains 0.001 for one doubling (< 0.002): c is acceptable, b is not (b -> c gains 0.005)
    assert picks["accurate"]["key"] == "c"


def test_fast_respects_the_competitor_limit():
    front = [pt("a", 1, 0.70), pt("b", 2, 0.71), pt("c", 4, 0.715), pt("d", 8, 0.716)]
    picks, info = S.knee_modes(front, fastest_other_cpu=10)   # limit 2 s
    assert math.isclose(info["fast_limit_cpu"], 2.0)
    assert picks["fast"]["key"] == "b"


def test_fast_falls_back_to_cheapest_point_and_flags_it():
    front = [pt("a", 5, 0.70), pt("b", 6, 0.71)]
    picks, info = S.knee_modes(front, fastest_other_cpu=10)   # limit 2 s, nothing within
    assert picks["fast"]["key"] == "a" and info["flag"]


def test_quality_is_mean_of_sp_and_tc():
    assert S.q({"sp": 0.8, "tc": 0.6}) == 0.7
