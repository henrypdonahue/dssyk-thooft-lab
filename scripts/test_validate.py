"""Tests for the snapshot runner's recipe parser and artifact comparison."""
import json
from pathlib import Path

import pytest

import validate

ROOT = Path(__file__).resolve().parent.parent


def _pair(tmp_path, old, new, suffix=".json"):
    a, b = tmp_path / f"old{suffix}", tmp_path / f"new{suffix}"
    if suffix == ".json":
        a.write_text(json.dumps(old))
        b.write_text(json.dumps(new))
    else:
        a.write_text(old)
        b.write_text(new)
    return a, b


def test_every_makefile_target_parses():
    for target in validate.TARGETS:
        recipes = validate.recipes(ROOT, target)
        assert recipes, target
        for module, argv in recipes:
            assert (ROOT / module).is_dir(), module
            assert argv


def test_identical_artifacts_compare_clean(tmp_path):
    data = {"a": [1.0, 2.0, 3.0], "b": {"c": 0.5, "d": "text"}}
    c = validate.compare_artifact(*_pair(tmp_path, data, data))
    assert c["max_rel"] == 0.0
    assert c["n_structure_changes"] == c["n_text_changes"] == 0


def test_runtime_fields_are_ignored(tmp_path):
    c = validate.compare_artifact(*_pair(tmp_path, {"x": 1.0, "runtime_s": 3.0},
                                         {"x": 1.0, "runtime_s": 9.0}))
    assert c["max_rel"] == 0.0


def test_numeric_text_and_structure_changes_are_reported(tmp_path):
    c = validate.compare_artifact(*_pair(
        tmp_path, {"x": [1.0, 2.0], "s": "p"}, {"x": [1.0, 2.1], "s": "q", "w": 1}))
    assert c["max_rel"] == pytest.approx(0.1 / 2.1)
    assert c["max_rel_at"] == "/x[1]"
    assert c["text_changes"] == ["/s"]
    assert c["structure_changes"] == ["/w"]


def test_roundoff_zeros_are_not_compared(tmp_path):
    """The dynamical histogram of an exactly conserved channel is roundoff
    (1e-22 against a total weight of 1); a change there is not a change."""
    old = {"total_weight": 1.0, "hist": [1e-22, 3e-23]}
    new = {"total_weight": 1.0, "hist": [4e-21, 1e-24]}
    assert validate.compare_artifact(*_pair(tmp_path, old, new))["max_rel"] == 0.0


def test_csv_compared_cellwise(tmp_path):
    c = validate.compare_artifact(*_pair(tmp_path, "n,v\n0,1.5\n1,2.5\n",
                                         "n,v\n0,1.5\n1,2.6\n", suffix=".csv"))
    assert c["kind"] == "csv"
    assert c["max_rel"] == pytest.approx(0.1 / 2.6)


def test_binary_files_are_flagged_not_compared(tmp_path):
    a, b = tmp_path / "a.png", tmp_path / "b.png"
    a.write_bytes(b"\x89PNG1")
    b.write_bytes(b"\x89PNG2")
    assert validate.compare_artifact(a, b)["kind"] == "binary"


def test_nan_and_inf_changes_are_reported(tmp_path):
    nan, inf = float("nan"), float("inf")
    for old, new in (({"x": 1.0}, {"x": nan}), ({"x": nan}, {"x": 1.0}),
                     ({"x": 1.0}, {"x": inf}), ({"a": [nan, 1.0, 2.0]}, {"a": [1.0, 1.0, 2.0]})):
        c = validate.compare_artifact(*_pair(tmp_path, old, new))
        assert c["n_text_changes"] == 1, (old, new)
    # an unchanged NaN is not a change, and does not hide a real one
    c = validate.compare_artifact(*_pair(tmp_path, {"a": [nan, 1.0, 2.0]}, {"a": [nan, 1.0, 5.0]}))
    assert c["n_text_changes"] == 0 and c["max_rel"] == pytest.approx(3.0 / 5.0)


def test_type_changes_are_structural(tmp_path):
    for old, new in (({"x": 1}, {"x": True}), ({"x": {"a": 1}}, {"x": [1]}),
                     ({"x": 1.0}, {"x": "1.0"})):
        c = validate.compare_artifact(*_pair(tmp_path, old, new))
        assert c["n_structure_changes"] == 1, (old, new)


def test_int_and_float_compare_as_numbers(tmp_path):
    c = validate.compare_artifact(*_pair(tmp_path, {"x": 2}, {"x": 2.0}))
    assert c["max_rel"] == 0.0 and c["n_structure_changes"] == 0


def test_runtime_key_patterns():
    for k in ("runtime_s", "gates_runtime_s", "runtime", "elapsed_s",
              "mean_seconds_per_instance", "seconds"):
        assert validate._is_runtime_key(k), k
    for k in ("mu0", "omega2", "lam"):
        assert not validate._is_runtime_key(k), k


def test_nested_records_and_length_changes(tmp_path):
    old = {"rows": [{"a": 1.0, "b": [1.0, 2.0]}, {"a": 2.0}]}
    new = {"rows": [{"a": 1.0, "b": [1.0, 2.2]}, {"a": 2.0}], "extra": 1}
    c = validate.compare_artifact(*_pair(tmp_path, old, new))
    assert c["max_rel_at"] == "/rows[0]/b[1]" and c["structure_changes"] == ["/extra"]
    c = validate.compare_artifact(*_pair(tmp_path, {"a": [1, 2]}, {"a": [1, 2, 3]}))
    assert c["n_structure_changes"] == 1


def test_unparseable_input(tmp_path):
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    a.write_text("{")
    b.write_text("{}")
    assert "unparseable" in validate.compare_artifact(a, b)["error"]


def test_unsupported_recipe_is_rejected(tmp_path):
    (tmp_path / "Makefile").write_text("test:\n\techo hi\n")
    with pytest.raises(ValueError):
        validate.recipes(tmp_path, "test")
