"""Tests for the IDT 96-well plate order export."""

import io
import sys
import zipfile
from datetime import date
from pathlib import Path

import pytest

_repo_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_repo_root / "src"))

from probedesign.idt import (
    ProbeSet,
    build_plates,
    gene_label,
    parse_oligos_text,
    plate_filename,
    read_oligo_sets,
    write_plates_xlsx,
)
from probedesign.cli import main
from click.testing import CliRunner

openpyxl = pytest.importorskip("openpyxl")

REGRESSION = _repo_root / "tests" / "regression"


def make_set(name, n, seq="acgtacgtacgtacgtac"):
    return ProbeSet(name, [(f"{name}_{i + 1}", seq) for i in range(n)])


def read_xlsx(data):
    wb = openpyxl.load_workbook(io.BytesIO(data))
    return {ws.title: [list(r) for r in ws.iter_rows(values_only=True)] for ws in wb.worksheets}


# ── Gene labels ──────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name, gene", [
    ("Eaat1-cmm", "Eaat1"),
    ("Hsc70-4-cmm", "Hsc70-4"),
    ("TDP-43-CDS", "TDP-43"),
    ("TDP-43-CDS-DmelCO", "TDP-43"),
    ("Eaat1_FBgn0026439_constitutive", "Eaat1"),
    ("ACTB_ENSG00000075624_constitutive", "ACTB"),
    ("KRT19_exon", "KRT19"),
    ("futsch-RA", "futsch"),
    ("Hsc70-4", "Hsc70-4"),
    ("CDKN1A_32", "CDKN1A_32"),
    ("cmm", "cmm"),
])
def test_gene_label(name, gene):
    assert gene_label(name) == gene


# ── Parsing oligo files ──────────────────────────────────────────────────────

def test_parse_smfish_oligos():
    text = (REGRESSION / "CDKN1A_32" / "CDKN1A_32_oligos.txt").read_text()
    pset = parse_oligos_text(text)
    assert pset.name == "CDKN1A_32"
    assert pset.oligos[0] == ("CDKN1A_32_1", "catgggttctgacggacatc")
    assert len(pset.oligos) == len(text.strip().splitlines())


def test_parse_hcr_oligos_skips_header():
    text = (
        "pair\thalf\tstart\tGC%\tTm\tGibbs\tstrict\tfull_oligo\tname\n"
        "1\tP1\t80\t48.0\t70.1\t-31.0\tyes\tGAGGAGGGCAGCAAACGGaaacgtacgtac\tKRT19_HCRB1_01\n"
        "1\tP2\t53\t52.0\t71.0\t-30.5\tyes\tacgtacgtacgtaGAAGAGTCTTCCTTTACG\tKRT19_HCRB1_02\n"
    )
    pset = parse_oligos_text(text)
    assert pset.name == "KRT19"
    assert [n for n, _ in pset.oligos] == ["KRT19_HCRB1_01", "KRT19_HCRB1_02"]


def test_parse_rejects_bad_sequence():
    with pytest.raises(ValueError, match="not a DNA sequence"):
        parse_oligos_text("1\t5\t50\t60\t-23\tacgt\tX_1\n2\t9\t50\t60\t-23\tnot-dna\tX_2\n")


def test_read_oligo_sets_from_batch_zip():
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("A-cmm_oligos.txt", "1\t1\t50\t60\t-23\tacgt\tA-cmm_1\n")
        zf.writestr("A-cmm_seq.txt", "# header\n")
        zf.writestr("B_HCR_oligos.txt", "pair\thalf\tstart\tGC%\tTm\tGibbs\tstrict\tfull_oligo\tname\n"
                                        "1\tP1\t1\t50\t60\t-31\tyes\tacgt\tB_HCRB2_01\n")
        zf.writestr("batch_summary.tsv", "x\n")
    sets = read_oligo_sets("results.zip", buf.getvalue())
    assert [s.name for s in sets] == ["A-cmm", "B"]


# ── Plate layout and packing ─────────────────────────────────────────────────

def test_wells_fill_down_columns_and_sets_start_new_column():
    plates, warnings = build_plates([make_set("A-cmm", 45), make_set("B-cmm", 10)])
    assert len(plates) == 1
    wells = [w for w, _, _ in plates[0].wells]
    assert wells[:9] == ["A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2"]
    assert wells[44] == "E6"          # last of A (45 = 5 cols + 5)
    assert wells[45] == "A7"          # B starts on a fresh column
    assert warnings == []


def test_set_spilling_over_moves_to_new_plate():
    # 4 x 32 oligos (4 columns each): 3 fit on plate 1, the 4th moves whole to plate 2
    sets = [make_set(g, 32) for g in ("Eaat1-cmm", "Eaat2-cmm", "futsch-cmm", "Hsc70-4-cmm")]
    plates, warnings = build_plates(sets, min_oligos=24)
    assert [[s.name for s in p.sets] for p in plates] == [
        ["Eaat1-cmm", "Eaat2-cmm", "futsch-cmm"], ["Hsc70-4-cmm"]]
    assert [p.sheet_name for p in plates] == ["Eaat1 Eaat2 futsch", "Hsc70-4"]
    assert warnings == []


def test_min_oligos_rebalances_plates():
    # With a 48 minimum, 3 + 1 leaves 32 on plate 2, so it packs 2 + 2
    sets = [make_set(g, 32) for g in ("Eaat1-cmm", "Eaat2-cmm", "futsch-cmm", "Hsc70-4-cmm")]
    plates, warnings = build_plates(sets, min_oligos=48)
    assert [p.sheet_name for p in plates] == ["Eaat1 Eaat2", "futsch Hsc70-4"]
    assert [len(p.wells) for p in plates] == [64, 64]
    assert warnings == []


def test_under_minimum_warns_when_unavoidable():
    plates, warnings = build_plates([make_set("A", 16)], min_oligos=24)
    assert len(plates) == 1
    assert any("IDT needs at least 24" in w for w in warnings)


def test_set_larger_than_plate_is_split_over_whole_plates():
    plates, warnings = build_plates([make_set("Big-cds", 120), make_set("Small-cds", 16)])
    assert [len(p.wells) for p in plates] == [96, 40]
    assert plates[1].wells[24][0] == "A4"   # Small starts after the 24-oligo remainder
    assert [p.sheet_name for p in plates] == ["Big", "Big Small"]
    assert any("split over 2 plates" in w for w in warnings)


def test_duplicate_oligo_names_warn():
    _, warnings = build_plates([make_set("A", 24), make_set("A", 24)])
    assert any("Duplicate oligo names" in w for w in warnings)


def test_sheet_names_are_unique_and_limited():
    sets = [make_set(f"Gene{i:02d}-cmm", 8) for i in range(12)] + \
           [make_set(f"Gene{i:02d}-cmm", 8) for i in range(12)]
    plates, _ = build_plates(sets, min_oligos=0)
    names = [p.sheet_name for p in plates]
    assert all(len(n) <= 31 for n in names)
    assert len({n.lower() for n in names}) == len(names)
    assert names[0].endswith("+8")          # 12 genes don't fit in 31 chars


# ── File name and workbook ───────────────────────────────────────────────────

def test_plate_filename():
    plates, _ = build_plates([make_set("TDP-43-CDS", 48), make_set("TDP-43-CDS-DmelCO", 45),
                              make_set("Su(H) cmm", 8)], min_oligos=0)
    assert plate_filename(plates, on=date(2026, 11, 2)) == \
        "IDT_plate-ordering_20261102_TDP-43_Su-H.xlsx"


def test_write_plates_xlsx_round_trip():
    plates, _ = build_plates([make_set("Eaat1-cmm", 32), make_set("Eaat2-cmm", 3)])
    sheets = read_xlsx(write_plates_xlsx(plates))
    rows = sheets["Eaat1 Eaat2"]
    assert rows[0] == ["Well Position", "Name", "Sequence"]
    assert rows[1] == ["A1", "Eaat1-cmm_1", "acgtacgtacgtacgtac"]
    assert rows[33] == ["A5", "Eaat2-cmm_1", "acgtacgtacgtacgtac"]
    assert len(rows) == 1 + 35


# ── CLI ──────────────────────────────────────────────────────────────────────

def test_cli_idtplate(tmp_path):
    files = []
    for gene in ("Eaat1-cmm", "Eaat2-cmm"):
        f = tmp_path / f"{gene}_oligos.txt"
        f.write_text("".join(f"{i}\t1\t50\t60\t-23\tacgtacgtacgtacgtac\t{gene}_{i}\n"
                             for i in range(1, 33)))
        files.append(str(f))
    out_dir = tmp_path / "out"
    result = CliRunner().invoke(main, ["idtplate", *files, "-o", str(out_dir)])
    assert result.exit_code == 0, result.output
    written = list(out_dir.glob("IDT_plate-ordering_*_Eaat1_Eaat2.xlsx"))
    assert len(written) == 1
    assert list(read_xlsx(written[0].read_bytes())) == ["Eaat1 Eaat2"]


def test_cli_idtplate_bad_file(tmp_path):
    f = tmp_path / "x_oligos.txt"
    f.write_text("\n")
    result = CliRunner().invoke(main, ["idtplate", str(f), "-o", str(tmp_path)])
    assert result.exit_code != 0
    assert "no oligos found" in result.output
