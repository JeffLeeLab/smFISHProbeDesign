"""IDT 96-well plate order files for designed oligo sets.

Writes the Excel upload IDT accepts for plate orders: one worksheet per plate,
columns Well Position / Name / Sequence, wells filled down each column
(A1..H1, A2..). Each probe set starts on a fresh column and stays on one
plate. Plates are packed in input order and each holds at least IDT's
minimum number of oligos where possible.
"""

import io
import math
import re
import zipfile
from collections import Counter
from dataclasses import dataclass
from datetime import date
from pathlib import PurePosixPath
from typing import Iterable, List, Optional, Sequence, Tuple

ROWS = "ABCDEFGH"
N_COLUMNS = 12
PLATE_WELLS = len(ROWS) * N_COLUMNS
DEFAULT_MIN_OLIGOS = 24   # standard 96-well plates; IDT HOTPlates need 48
HEADER = ("Well Position", "Name", "Sequence")
SHEET_NAME_MAX = 31       # Excel limit
FILENAME_GENES_MAX = 100

# Name tokens that mark the targeted region rather than the gene,
# e.g. Eaat1-cmm, TDP-43-CDS, Eaat1_FBgn0026439_constitutive.
_REGION_TOKENS = {
    "cmm", "cds", "orf", "exon", "exons", "exonic", "intron", "introns",
    "intronic", "utr", "5utr", "3utr", "5'utr", "3'utr", "utr5", "utr3",
    "constitutive", "mrna", "premrna", "pre", "nascent", "transcript",
}
_ID_TOKEN = re.compile(
    r"(ENS[A-Z]*[GT]\d+(\.\d+)?|FBgn\d+|FBtr\d+|WBGene\d+|NM|NR|XM|XR|R[A-Z])"
)
_SEQ = re.compile(r"[ACGTUNacgtun]+")
_PROBE_NUMBER = re.compile(r"(_HCR[A-Za-z0-9]+)?_\d+$")
_SHEET_FORBIDDEN = re.compile(r"[\[\]:*?/\\]")


@dataclass
class ProbeSet:
    """One designed probe set: (oligo name, sequence) pairs in plate order."""

    name: str
    oligos: List[Tuple[str, str]]
    gene: str = ""

    def __post_init__(self):
        if not self.gene:
            self.gene = gene_label(self.name)

    @property
    def n_columns(self) -> int:
        return math.ceil(len(self.oligos) / len(ROWS))


@dataclass
class Plate:
    """One worksheet of the order: its probe sets and filled wells."""

    sheet_name: str
    sets: List[ProbeSet]
    wells: List[Tuple[str, str, str]]   # (well, oligo name, sequence)


def gene_label(set_name: str) -> str:
    """Gene part of a probe-set name: everything before the first region or ID token.

    'Eaat1-cmm' -> 'Eaat1', 'TDP-43-CDS-DmelCO' -> 'TDP-43',
    'Eaat1_FBgn0026439_constitutive' -> 'Eaat1'. A name with no such token is
    returned whole.
    """
    parts = re.split(r"([-_ ])", set_name)
    tokens = parts[0::2]
    for i in range(1, len(tokens)):
        tok = tokens[i]
        if tok.lower() in _REGION_TOKENS or _ID_TOKEN.fullmatch(tok):
            gene = "".join(parts[:2 * i - 1]).rstrip("-_ ")
            return gene or set_name
    return set_name


def set_name_from_oligo_names(names: Sequence[str]) -> str:
    """Probe-set name from its oligo names ('KRT19_12' / 'KRT19_HCRB1_03' -> 'KRT19')."""
    return _PROBE_NUMBER.sub("", names[0]) or names[0]


def parse_oligos_text(text: str, source: str = "oligos") -> ProbeSet:
    """Probe set from a ProbeDesign _oligos.txt or _HCR_oligos.txt.

    Sequence and name are the last two tab-separated columns in both formats;
    a header row (HCR) is skipped.
    """
    oligos = []
    for lineno, line in enumerate(text.splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise ValueError(f"{source}: line {lineno} is not tab-separated")
        seq, name = fields[-2].strip(), fields[-1].strip()
        if not _SEQ.fullmatch(seq):
            if not oligos:
                continue   # header row
            raise ValueError(f"{source}: line {lineno}: {seq!r} is not a DNA sequence")
        if not name:
            raise ValueError(f"{source}: line {lineno} has no oligo name")
        oligos.append((name, seq))
    if not oligos:
        raise ValueError(f"{source}: no oligos found")
    return ProbeSet(set_name_from_oligo_names([n for n, _ in oligos]), oligos)


def read_oligo_sets(filename: str, data: bytes) -> List[ProbeSet]:
    """Probe sets from one file: an oligos .txt, or a .zip of them (batch results)."""
    if filename.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as zf:
            members = [
                n for n in zf.namelist()
                if n.endswith("_oligos.txt") and not n.startswith("__MACOSX")
            ]
            if not members:
                raise ValueError(f"{filename}: no *_oligos.txt files inside")
            return [
                parse_oligos_text(zf.read(n).decode("utf-8-sig"), PurePosixPath(n).name)
                for n in members
            ]
    return [parse_oligos_text(data.decode("utf-8-sig"), filename)]


def _partition(pieces: List[ProbeSet], min_oligos: int) -> Optional[List[List[ProbeSet]]]:
    """Cut the ordered pieces into plates: fewest plates, then earlier plates fullest.

    Every plate fits in 12 columns and holds at least min_oligos. None if no
    cut does.
    """
    n = len(pieces)
    # best[i]: (plates, (-filled per plate...), plate end indices) for pieces[i:]
    best: List[Optional[tuple]] = [None] * n + [(0, (), ())]
    for i in range(n - 1, -1, -1):
        cols = filled = 0
        for j in range(i + 1, n + 1):
            cols += pieces[j - 1].n_columns
            filled += len(pieces[j - 1].oligos)
            if cols > N_COLUMNS:
                break
            if filled < min_oligos or best[j] is None:
                continue
            plates, key, ends = best[j]
            cand = (plates + 1, (-filled,) + key, (j,) + ends)
            if best[i] is None or cand[:2] < best[i][:2]:
                best[i] = cand
    if best[0] is None:
        return None
    groups, start = [], 0
    for end in best[0][2]:
        groups.append(pieces[start:end])
        start = end
    return groups


def _layout(group: List[ProbeSet]) -> List[Tuple[str, str, str]]:
    """Wells down each column, each set starting on a fresh column."""
    wells = []
    col = 0
    for pset in group:
        for k, (name, seq) in enumerate(pset.oligos):
            c, r = divmod(k, len(ROWS))
            wells.append((f"{ROWS[r]}{col + c + 1}", name, seq))
        col += pset.n_columns
    return wells


def _unique(items: Iterable[str]) -> List[str]:
    return list(dict.fromkeys(items))


def _join_limited(labels: List[str], limit: int, sep: str) -> str:
    """Join labels, dropping trailing ones as '+N' to stay within limit."""
    full = sep.join(labels)
    if len(full) <= limit:
        return full
    for k in range(len(labels) - 1, 0, -1):
        short = f"{sep.join(labels[:k])}{sep}+{len(labels) - k}"
        if len(short) <= limit:
            return short
    return labels[0][:limit]


def _sheet_name(group: List[ProbeSet], taken: set) -> str:
    """Worksheet name from the plate's genes; Excel-safe, unique, <= 31 chars."""
    genes = _unique(_SHEET_FORBIDDEN.sub("-", s.gene).strip("' ") for s in group)
    base = _join_limited([g for g in genes if g], SHEET_NAME_MAX, " ") or "Plate"
    name, n = base, 1
    while name.lower() in taken:
        n += 1
        suffix = f" {n}"
        name = base[:SHEET_NAME_MAX - len(suffix)] + suffix
    taken.add(name.lower())
    return name


def build_plates(
    sets: List[ProbeSet],
    min_oligos: int = DEFAULT_MIN_OLIGOS,
) -> Tuple[List[Plate], List[str]]:
    """Assign probe sets to plates in order, without splitting a set.

    Each set starts on a fresh column. Of the ways to cut the ordered sets into
    plates of at least min_oligos, picks the fewest plates, filling earlier
    plates first. A set larger than a plate is split over whole plates. If no
    cut meets the minimum, packs without it and warns about the short plates.

    Returns:
        (plates, warnings)
    """
    if not sets:
        raise ValueError("no probe sets to place")
    warnings = []
    pieces = []
    for pset in sets:
        if len(pset.oligos) > PLATE_WELLS:
            n_plates = math.ceil(len(pset.oligos) / PLATE_WELLS)
            warnings.append(
                f"{pset.name} has {len(pset.oligos)} oligos, more than one plate holds; "
                f"it is split over {n_plates} plates."
            )
            pieces.extend(
                ProbeSet(pset.name, pset.oligos[i:i + PLATE_WELLS], pset.gene)
                for i in range(0, len(pset.oligos), PLATE_WELLS)
            )
        else:
            pieces.append(pset)

    groups = _partition(pieces, min_oligos) or _partition(pieces, 0)
    taken: set = set()
    plates = [Plate(_sheet_name(g, taken), g, _layout(g)) for g in groups]

    for plate in plates:
        if len(plate.wells) < min_oligos:
            warnings.append(
                f"Plate '{plate.sheet_name}' has {len(plate.wells)} oligos; "
                f"IDT needs at least {min_oligos} per plate."
            )
    counts = Counter(name for plate in plates for _, name, _ in plate.wells)
    dupes = [n for n, c in counts.items() if c > 1]
    if dupes:
        warnings.append(f"Duplicate oligo names: {', '.join(dupes[:5])}"
                        + (" ..." if len(dupes) > 5 else ""))
    return plates, warnings


def plate_filename(plates: List[Plate], on: Optional[date] = None) -> str:
    """IDT_plate-ordering_<YYYYMMDD>_<gene>_<gene>....xlsx"""
    genes = _unique(
        re.sub(r"[^A-Za-z0-9.+-]+", "-", s.gene).strip("-.")
        for plate in plates for s in plate.sets
    )
    part = _join_limited([g for g in genes if g], FILENAME_GENES_MAX, "_") or "probes"
    return f"IDT_plate-ordering_{(on or date.today()):%Y%m%d}_{part}.xlsx"


def write_plates_xlsx(plates: List[Plate]) -> bytes:
    """The order workbook as .xlsx bytes, one worksheet per plate."""
    try:
        from openpyxl import Workbook
    except ImportError as e:
        raise ImportError(
            "IDT plate export needs openpyxl: pip install openpyxl"
        ) from e

    wb = Workbook()
    wb.remove(wb.active)
    for plate in plates:
        ws = wb.create_sheet(plate.sheet_name)
        ws.append(HEADER)
        for row in plate.wells:
            ws.append(row)
        longest = max((len(n) for _, n, _ in plate.wells), default=10)
        seq_len = max((len(s) for _, _, s in plate.wells), default=20)
        ws.column_dimensions["A"].width = 14
        ws.column_dimensions["B"].width = max(14, longest + 2)
        ws.column_dimensions["C"].width = max(14, seq_len + 2)
    buf = io.BytesIO()
    wb.save(buf)
    return buf.getvalue()
