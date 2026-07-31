"""
Pinned ribosomal-protein panel, shared by every script that assigns gene class.

WHY THIS EXISTS
---------------
The panel was previously a regular expression. Regexes on gene symbols do not
track curation: `^(RPL|RPS|MRPL|MRPS)\\d` silently admitted ten RPS6K* kinases
(ribosomal protein S6 kinases, which are signalling enzymes) and RPS19BP1,
carried the obsolete symbol MRPS36, and omitted seven genuine members whose
approved symbols do not start with a ribosomal prefix at all: AURKAIP1, CHCHD1,
DAP3, FAU, GADD45GIP1, PTCD3 and UBA52. Patching the exceptions one at a time
was replacing a curation problem with a maintenance problem.

The panel is now a static committed CSV built from three HGNC curated gene
groups, with retrieval date, source URLs and a SHA-256 recorded alongside. The
SHA-256 is written into every matched-control cache so a null can be tied to
the exact panel it was drawn under.

Membership is resolved by **Ensembl gene ID**, not symbol. MRPS36 is the reason:
a stale symbol matches nothing and disappears silently, whereas a stale
identifier can be detected. Symbols are retained in the CSV for readability and
are cross-checked on load.

Usage:
    from _ribosomal_panel import ribosomal_symbols, panel_provenance
    RIBO = ribosomal_symbols(table_s1)      # set of upper-case symbols
    if sym.upper() in RIBO: ...
"""

import hashlib
import json
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
if _HERE.parent.name == "analysis":
    _DATA = _HERE.parent.parent / "data"
else:
    _DATA = _HERE.parent.parent.parent / "notebooks" / "data"

PANEL_PATH = _DATA / "ribosomal_panel.csv"
PROVENANCE_PATH = _DATA / "ribosomal_panel_provenance.json"

if not PANEL_PATH.exists():
    raise FileNotFoundError(
        f"Ribosomal panel not found at {PANEL_PATH}. It is a committed data "
        f"file, not a generated one; restore it from the repository.")

PANEL = pd.read_csv(PANEL_PATH)
PANEL_SHA256 = hashlib.sha256(PANEL_PATH.read_bytes()).hexdigest()

_expected = json.loads(PROVENANCE_PATH.read_text()) if PROVENANCE_PATH.exists() else {}
if _expected.get("panel_sha256") not in (None, PANEL_SHA256):
    raise RuntimeError(
        f"{PANEL_PATH.name} hashes to {PANEL_SHA256[:16]} but its provenance "
        f"file records {_expected['panel_sha256'][:16]}. The panel has been "
        f"edited without updating its provenance; refusing to proceed.")


def _strip_version(s):
    return s.astype(str).str.split(".").str[0]


def ribosomal_symbols(table_s1: pd.DataFrame) -> set:
    """Upper-case symbols in the annotation table that are panel members.

    Resolution is by Ensembl gene ID. Any panel member that cannot be resolved
    that way falls back to its approved symbol, and the fallback is reported
    rather than applied silently.
    """
    ens = set(PANEL["ensembl_gene_id"].astype(str))
    tab_ens = _strip_version(table_s1["ensembl_id"])
    hit = table_s1.loc[tab_ens.isin(ens), "gene_symbol"].astype(str).str.upper()
    resolved = set(hit)

    unresolved = ens - set(tab_ens)
    if unresolved:
        miss = PANEL.loc[PANEL["ensembl_gene_id"].isin(unresolved),
                         "approved_symbol"].astype(str).str.upper()
        by_symbol = set(miss) & set(
            table_s1["gene_symbol"].astype(str).str.upper())
        if by_symbol:
            print(f"    [ribosomal panel] {len(by_symbol)} member(s) not "
                  f"resolvable by Ensembl ID, matched by symbol instead: "
                  f"{sorted(by_symbol)}")
            resolved |= by_symbol
        still = set(miss) - by_symbol
        if still:
            print(f"    [ribosomal panel] {len(still)} member(s) absent from "
                  f"the annotation table: {sorted(still)}")
    return resolved


def panel_provenance() -> dict:
    """Recorded in every output so a result ties to the panel that produced it."""
    return {"ribosomal_panel": PANEL_PATH.name,
            "ribosomal_panel_sha256": PANEL_SHA256,
            "ribosomal_panel_n": int(len(PANEL)),
            "ribosomal_panel_source": "HGNC gene groups 728, 729, 646",
            "ribosomal_panel_retrieved": _expected.get("retrieved_utc_date")}
