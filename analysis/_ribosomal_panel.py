"""Shared loader for the pinned ribosomal gene panel.

Resolves genes by Ensembl ID with a reported-symbol fallback, and refuses to
load if the panel CSV's hash disagrees with its provenance file, so an analysis
cannot silently run against an edited panel.

Inputs:   data/ribosomal_panel.csv, data/ribosomal_panel_provenance.json
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
