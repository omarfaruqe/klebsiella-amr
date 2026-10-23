import pytest

from kamr.curation.breakpoints import Breakpoint, interpret_mic, normalize_unit

MERO = Breakpoint("meropenem", s_le=2, r_gt=8, applies=True, version="EUCAST test")
GENT = Breakpoint("gentamicin", s_le=2, r_gt=2, applies=True, version="EUCAST test")


@pytest.mark.parametrize(
    ("sign", "value", "expected"),
    [
        ("=", "2", "S"),
        ("", "0.5", "S"),
        ("=", "4", "I"),
        ("=", "8", "I"),
        ("=", "16", "R"),
        ("<=", "0.25", "S"),
        ("<=", "4", None),
        (">", "8", "R"),
        (">", "4", None),
        (">=", "16", "R"),
        (">=", "8", None),
        ("=", "not-a-number", None),
    ],
)
def test_interpret_meropenem(sign, value, expected):
    assert interpret_mic(sign, value, MERO) == expected


def test_no_intermediate_when_breakpoints_equal():
    assert interpret_mic("=", "2", GENT) == "S"
    assert interpret_mic("=", "4", GENT) == "R"


def test_normalize_unit():
    assert normalize_unit("mg/L") == "mg/l"
    assert normalize_unit("µg/mL") == "ug/ml"
    assert normalize_unit(None) == ""
