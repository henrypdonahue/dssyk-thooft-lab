"""The committed WALKTHROUGH.html is exactly what the generator produces
from the committed template and data (regenerate with
python scripts/walkthrough/make_walkthrough.py)."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def _generator():
    path = ROOT / "scripts" / "walkthrough" / "make_walkthrough.py"
    spec = importlib.util.spec_from_file_location("make_walkthrough", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_walkthrough_is_current():
    html = _generator().build_html()
    assert html == (ROOT / "WALKTHROUGH.html").read_text(), \
        "WALKTHROUGH.html is stale: run scripts/walkthrough/make_walkthrough.py"


def test_errata_count_matches_the_last_entry():
    gen = _generator()
    text = (ROOT / "ERRATA.md").read_text()
    assert f"**E{gen.errata_count()}." in text
    assert f"**E{gen.errata_count() + 1}." not in text
