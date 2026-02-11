#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find interface residues ("hotspots") on a target protein chain(s) given a complex structure
(PDB or mmCIF). Designed for RFdiffusion PPI binder design.

Key features:
- Robust parsing with gemmi (PDB + mmCIF) with insertion codes (icode) support
- Heavy-atom interface scoring with residue-type weights (biological hotspot prior)
- Hotspot selection: spread or patch
- Outputs:
  - ranked_interface_residues.tsv
  - selected_hotspots.txt (RFdiffusion-ready string)
  - hotspots.pml (PyMOL script)

Usage example (kept compatible with the original script):
python find_interface_hotspots.py --pdb complex.cif --target_chains A --partner_chains B --select 6 --mode spread --outdir output/hotspot_out
"""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import gemmi
import numpy as np


# -----------------------------
# Types
# -----------------------------
ResKey = Tuple[str, int, str]  # (chain_id, resseq, icode)
ChainRanges = Dict[str, Tuple[int, int]]


class AtomRec(dict):
    """A light-weight record for an atom.

    Keys:
      chain (str), resseq (int), icode (str), resname (str), atomname (str),
      element (str), xyz (np.ndarray shape=(3,)), reskey (ResKey)
    """
    pass


# -----------------------------
# Optional KDTree (kept)
# -----------------------------
def try_import_ckdtree():
    """Try importing scipy's cKDTree for fast neighbor search.

    Kept for backwards compatibility: if scipy isn't available, we fall back
    to a pure-numpy grid bucketing search.
    """
    try:
        from scipy.spatial import cKDTree  # type: ignore
        return cKDTree
    except Exception:
        return None


def build_kdtree(coords: np.ndarray):
    cKDTree = try_import_ckdtree()
    if cKDTree is None:
        return None, None
    return cKDTree(coords), "scipy"


def grid_neighbor_search(query_coords: np.ndarray, ref_coords: np.ndarray, cutoff: float) -> List[List[int]]:
    """Pure numpy grid bucketing neighbor search (fallback if scipy not installed).

    Returns a list of neighbor indices for each query point.
    """
    cell = float(cutoff)
    grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)

    ref_cell = np.floor(ref_coords / cell).astype(int)
    for idx, c in enumerate(ref_cell):
        grid[(int(c[0]), int(c[1]), int(c[2]))].append(idx)

    q_cell = np.floor(query_coords / cell).astype(int)
    neighbors: List[List[int]] = []
    cutoff2 = cutoff * cutoff

    for qc, qxyz in zip(q_cell, query_coords):
        cand: List[int] = []
        qcx, qcy, qcz = int(qc[0]), int(qc[1]), int(qc[2])

        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    key = (qcx + dx, qcy + dy, qcz + dz)
                    if key in grid:
                        cand.extend(grid[key])

        if not cand:
            neighbors.append([])
            continue

        cand_idx = np.asarray(cand, dtype=int)
        cand_coords = ref_coords[cand_idx]
        d2 = np.sum((cand_coords - qxyz) ** 2, axis=1)
        mask = d2 <= cutoff2
        neighbors.append(cand_idx[mask].tolist())

    return neighbors


# -----------------------------
# Robust parsing (gemmi-only)
# -----------------------------
def _normalize_cli_path(p: str) -> str:
    """Normalize a CLI path string.

    This makes the script tolerant to users pasting a Python raw-string literal
    *directly* into the shell, e.g.:
      --pdb r"C:\\path with spaces\\file.cif"

    In real shells, you normally just quote it:
      --pdb "C:\\path with spaces\\file.cif"
    """
    s = p.strip()

    # If user pasted Python raw literal: r"..." or r'...'
    # Note: shell removes quotes, so argv could look like: rC:\path\file.cif
    if len(s) >= 3 and s[0] in ("r", "R") and s[1] == '"' and s[-1] == '"':
        s = s[2:-1]
    elif len(s) >= 3 and s[0] in ("r", "R") and s[1] == "'" and s[-1] == "'":
        s = s[2:-1]
    elif len(s) >= 3 and s[0] in ("r", "R") and re.match(r"^[rR][A-Za-z]:[\\/]", s):
        # common Windows argv: rC:\...
        s = s[1:]

    # Strip outer quotes if any remain
    if (len(s) >= 2) and ((s[0] == s[-1] == '"') or (s[0] == s[-1] == "'")):
        s = s[1:-1]

    return s


def _icode_to_str(icode: str) -> str:
    """Convert gemmi icode char to a safe string ('' if empty)."""
    if not icode:
        return ""
    if icode == "\0" or icode == " ":
        return ""
    return str(icode)


def parse_structure_atoms(path: str, keep_hetatm: bool = False) -> Tuple[List[AtomRec], Dict[ResKey, np.ndarray], Dict[ResKey, str], ChainRanges]:
    """Parse atoms from a structure file (PDB or mmCIF) using gemmi only.

    Returns:
      atoms: List[AtomRec] (heavy atoms only)
      residue_ca: ResKey -> CA xyz
      residue_name: ResKey -> 3-letter residue name (as in file)
      chain_ranges: chain_id -> (min_resseq, max_resseq) ignoring insertion codes
    """
    atoms: List[AtomRec] = []
    residue_ca: Dict[ResKey, np.ndarray] = {}
    residue_name: Dict[ResKey, str] = {}
    chain_resseqs: Dict[str, List[int]] = defaultdict(list)

    st = gemmi.read_structure(path)
    # Keep only the first model (typical for AF/ESMFold; still safe for multi-model)
    model = st[0]

    # Ensure altlocs are handled consistently: keep a single conformer
    # (similar to "altLoc in (' ', 'A')" logic)
    try:
        st.remove_alternative_conformations()
        model = st[0]
    except Exception:
        # If gemmi fails for some edge case, proceed without removal;
        # we will still filter altloc per-atom below.
        pass

    for chain in model:
        ch = chain.name if chain.name else "_"
        for res in chain:
            if res is None or not res.name:
                continue
            if (not keep_hetatm) and str(res.het_flag) == 'H':
                # Exclude HETATM residues (ligands, waters, etc.) by default
                continue

            resseq = int(res.seqid.num)
            icode = _icode_to_str(res.seqid.icode)
            rk: ResKey = (ch, resseq, icode)

            residue_name[rk] = res.name
            chain_resseqs[ch].append(resseq)

            for atom in res:
                # altloc safety: keep only '' or 'A'
                alt = getattr(atom, "altloc", "")
                if alt not in ("", " ", "A", "\0"):
                    continue

                # filter hydrogens
                el = atom.element.name.upper() if getattr(atom, "element", None) else ""
                if not el:
                    # fallback: infer from atom name
                    el = atom.name.strip()[:1].upper() if atom.name else ""
                if el.startswith("H") or atom.name.strip().upper().startswith("H"):
                    continue

                pos = atom.pos
                xyz = np.array([pos.x, pos.y, pos.z], dtype=np.float32)

                atoms.append(AtomRec(
                    chain=ch,
                    resseq=resseq,
                    icode=icode,
                    resname=res.name,
                    atomname=atom.name.strip(),
                    element=el,
                    xyz=xyz,
                    reskey=rk
                ))

                if atom.name.strip() == "CA":
                    residue_ca[rk] = xyz

    chain_ranges: ChainRanges = {}
    for ch, resseqs in chain_resseqs.items():
        if resseqs:
            chain_ranges[ch] = (int(min(resseqs)), int(max(resseqs)))

    return atoms, residue_ca, residue_name, chain_ranges


# -----------------------------
# Residue mapping + weights
# -----------------------------
_resname_to_aa1_cache: Dict[str, str] = {}


def resname_to_aa1(resname3: str) -> str:
    """Map a residue 3-letter (or component) name to a 1-letter amino-acid code.

    Uses gemmi's internal chemical component table to handle non-standard residues.
    Returns 'X' if unknown or not an amino-acid.
    """
    key = (resname3 or "").strip().upper()
    if key in _resname_to_aa1_cache:
        return _resname_to_aa1_cache[key]

    aa1 = "X"
    try:
        info = gemmi.find_tabulated_residue(key)
        if info and info.found() and info.is_amino_acid():
            # one_letter_code may be lower-case for modified residues; normalize to upper.
            olc = (info.one_letter_code or "X").upper()
            # gemmi uses 'X' when unknown; keep it.
            aa1 = olc if len(olc) == 1 else "X"
    except Exception:
        aa1 = "X"

    _resname_to_aa1_cache[key] = aa1
    return aa1


def residue_weight(aa1: str) -> float:
    """Residue-type weights for hotspot prior."""
    a = (aa1 or "X").upper()
    if a in {"W", "Y", "F"}:
        return 1.2
    if a in {"L", "I", "M"}:
        return 1.1
    if a in {"G", "P", "A"}:
        return 0.9
    return 1.0


# -----------------------------
# ID safety (RFdiffusion + PyMOL)
# -----------------------------
def _validate_chain_for_rfdiffusion(chain: str) -> None:
    """RFdiffusion hotspot tokens typically assume 1-char chain IDs (PDB-style).

    We enforce this to avoid ambiguous tokens such as 'AA100' or 'A1100'.
    """
    if len(chain) != 1:
        raise ValueError(
            f"RFdiffusion hotspot tokens assume 1-character chain IDs, but got chain '{chain}'. "
            f"Please rename chains to single characters (e.g., in PyMOL: alter all, chain='A'; sort) "
            f"or export as a PDB with 1-char chain IDs."
        )
    if chain in {" ", "\t", "\n"}:
        raise ValueError("Invalid chain ID for RFdiffusion.")


def residue_id_for_outputs(chain: str, resseq: int, icode: str) -> str:
    """Build a RFdiffusion-style residue token, e.g. 'A100' or 'A100A'."""
    _validate_chain_for_rfdiffusion(chain)
    if icode:
        if len(icode) != 1:
            raise ValueError(f"Invalid insertion code '{icode}' (expected 1 char).")
        return f"{chain}{resseq}{icode}"
    return f"{chain}{resseq}"


def pymol_resi_token(resseq: int, icode: str) -> str:
    """PyMOL resi token uses insertion code like '95A' if present."""
    return f"{resseq}{icode}" if icode else str(resseq)


def pymol_chain_token(chain: str) -> str:
    """Quote chain ID for PyMOL selection safely."""
    # Escape double quotes for safety
    ch = chain.replace('"', r'\"')
    return f'"{ch}"'


# -----------------------------
# Exclude residues parsing
# -----------------------------
_EXCL_RE_SIMPLE = re.compile(r"^([A-Za-z]+)(-?\d+)([A-Za-z]?)$")


def parse_exclude_residues(spec: Optional[str]) -> Set[ResKey]:
    """Parse --exclude_residues into a set of ResKey.

    Accepts:
      - PDB-like compact tokens: A10, A10A, B50
      - Safer tokens for multi-char chain IDs: AA:10, AA:10A

    Notes:
      - insertion code (icode) is optional and assumed to be at most 1 character.
      - If chain IDs contain digits or special chars, please use the ':' form.
    """
    out: Set[ResKey] = set()
    if not spec:
        return out

    for raw in spec.split(","):
        tok = raw.strip()
        if not tok:
            continue

        chain = ""
        resseq_s = ""
        icode = ""

        if ":" in tok:
            # AA:10A or AA:10
            chain, rest = tok.split(":", 1)
            chain = chain.strip()
            rest = rest.strip()
            m = re.match(r"^(-?\d+)([A-Za-z]?)$", rest)
            if not m:
                raise ValueError(f"Invalid exclude token '{tok}'. Expected like 'AA:10' or 'AA:10A'.")
            resseq_s = m.group(1)
            icode = m.group(2) or ""
        else:
            m = _EXCL_RE_SIMPLE.match(tok)
            if not m:
                raise ValueError(f"Invalid exclude token '{tok}'. Expected like 'A10' or 'A10A' (or 'AA:10').")
            chain = m.group(1)
            resseq_s = m.group(2)
            icode = m.group(3) or ""

        resseq = int(resseq_s)
        out.add((chain, resseq, icode))

    return out


# -----------------------------
# Scoring + selection
# -----------------------------
def score_interface(
    atoms: Sequence[AtomRec],
    residue_ca: Dict[ResKey, np.ndarray],
    residue_name: Dict[ResKey, str],
    target_chains: Sequence[str],
    partner_chains: Sequence[str],
    cutoff: float = 4.5,
    sigma: float = 2.0,
    exclude: Optional[Set[ResKey]] = None,
):
    """Compute per-target-residue interface scores based on heavy-atom contacts.

    Geometry score per residue:
      sum(exp(-(d/sigma)^2)) over atom-atom contacts within cutoff.

    Weighted score per residue:
      weighted_score = geometry_score * residue_weight(aa1)

    Also tracks:
      n_contacts (atom-atom pairs)
      n_partner_res (unique partner residues contacted)
      min_dist
    """
    exclude_set: Set[ResKey] = exclude or set()

    target_atoms = [a for a in atoms if a["chain"] in target_chains]
    partner_atoms = [a for a in atoms if a["chain"] in partner_chains]

    if len(target_atoms) == 0:
        raise ValueError(f"No atoms found for target chains: {target_chains}")
    if len(partner_atoms) == 0:
        raise ValueError(f"No atoms found for partner chains: {partner_chains}")

    t_coords = np.vstack([a["xyz"] for a in target_atoms]).astype(np.float32, copy=False)
    p_coords = np.vstack([a["xyz"] for a in partner_atoms]).astype(np.float32, copy=False)

    # neighbor search
    tree, mode = build_kdtree(p_coords)
    if tree is not None:
        neigh = tree.query_ball_point(t_coords, r=float(cutoff))
        neighbor_mode = f"KDTree({mode})"
    else:
        neigh = grid_neighbor_search(t_coords, p_coords, cutoff=float(cutoff))
        neighbor_mode = "grid-fallback"

    score_geom: Dict[ResKey, float] = defaultdict(float)
    n_contacts: Dict[ResKey, int] = defaultdict(int)
    min_dist: Dict[ResKey, float] = defaultdict(lambda: float("inf"))
    partner_res_sets: Dict[ResKey, Set[ResKey]] = defaultdict(set)

    p_reskeys = [a["reskey"] for a in partner_atoms]

    for i, t_atom in enumerate(target_atoms):
        rk = t_atom["reskey"]
        if rk in exclude_set:
            continue

        idxs = neigh[i]
        if not idxs:
            continue

        idx_arr = np.asarray(idxs, dtype=int)
        diffs = p_coords[idx_arr] - t_coords[i]
        dists = np.sqrt(np.sum(diffs * diffs, axis=1))

        n_contacts[rk] += int(len(idxs))
        md = float(np.min(dists))
        if md < min_dist[rk]:
            min_dist[rk] = md

        score_geom[rk] += float(np.sum(np.exp(- (dists / float(sigma)) ** 2)))

        for j in idxs:
            partner_res_sets[rk].add(p_reskeys[j])

    rows: List[dict] = []
    for rk, geom_sc in score_geom.items():
        ch, rs, ic = rk
        res3 = residue_name.get(rk, "UNK")
        aa1 = resname_to_aa1(res3)
        w = residue_weight(aa1)
        sc = float(geom_sc * w)

        rows.append({
            "chain": ch,
            "resseq": int(rs),
            "icode": ic,
            "resname": res3,
            "score": sc,  # weighted final score
            "n_contacts": int(n_contacts.get(rk, 0)),
            "n_partner_res": int(len(partner_res_sets.get(rk, set()))),
            "min_dist": float(min_dist.get(rk, float("inf"))),
            "has_ca": rk in residue_ca,
        })

    rows.sort(key=lambda r: (r["score"], r["n_contacts"]), reverse=True)
    return rows, neighbor_mode


def select_hotspots(
    rows: Sequence[dict],
    residue_ca: Dict[ResKey, np.ndarray],
    k: int = 6,
    mode: str = "spread",
    min_sep: float = 7.0,
    patch_radius: float = 12.0,
    exclude: Optional[Set[ResKey]] = None,
) -> List[dict]:
    """Pick K hotspot residues from ranked rows.

    mode:
      - spread: greedy pick ranked residues but enforce CA-CA >= min_sep between selected
      - patch: pick top residue then pick next best within patch_radius of that residue
    """
    exclude_set = exclude or set()

    candidates = []
    for r in rows:
        rk = (r["chain"], r["resseq"], r["icode"])
        if rk in exclude_set:
            continue
        if not r.get("has_ca", False):
            continue
        candidates.append(r)

    if not candidates:
        return []

    if mode == "patch":
        seed = candidates[0]
        seed_key = (seed["chain"], seed["resseq"], seed["icode"])
        seed_ca = residue_ca[seed_key]
        near: List[dict] = []
        for r in candidates:
            rk = (r["chain"], r["resseq"], r["icode"])
            ca = residue_ca[rk]
            d = float(np.linalg.norm(ca - seed_ca))
            if d <= float(patch_radius):
                near.append(r)
        near.sort(key=lambda r: (r["score"], r["n_contacts"]), reverse=True)
        return near[: int(k)]

    # spread mode
    selected: List[dict] = []
    selected_ca: List[np.ndarray] = []

    for r in candidates:
        rk = (r["chain"], r["resseq"], r["icode"])
        ca = residue_ca[rk]
        ok = True
        for sca in selected_ca:
            if float(np.linalg.norm(ca - sca)) < float(min_sep):
                ok = False
                break
        if ok:
            selected.append(r)
            selected_ca.append(ca)
            if len(selected) >= int(k):
                break

    # if not enough due to separation, fill without separation
    if len(selected) < int(k):
        for r in candidates:
            if r in selected:
                continue
            selected.append(r)
            if len(selected) >= int(k):
                break

    return selected


def make_rf_hotspot_string(selected: Sequence[dict]) -> str:
    toks: List[str] = []
    for r in selected:
        toks.append(residue_id_for_outputs(r["chain"], int(r["resseq"]), str(r["icode"])))
    return "ppi.hotspot_res=[" + ",".join(toks) + "]"


def make_contig_suggestion(chain_ranges: ChainRanges, target_chains: Sequence[str], binder_len: str = "80-120") -> Optional[str]:
    segs: List[str] = []
    for ch in target_chains:
        if ch not in chain_ranges:
            continue
        mn, mx = chain_ranges[ch]
        segs.append(f"{ch}{mn}-{mx}")
    if not segs:
        return None
    return "contigmap.contigs=[" + " ".join(segs) + f"/0 {binder_len}" + "]"


def write_pymol_script(out_pml: str, pdb_basename: str, selected: Sequence[dict], target_chains: Sequence[str]) -> None:
    """Write a small PyMOL script to visualize selected hotspots."""
    chain_to_resi: Dict[str, List[str]] = defaultdict(list)
    for r in selected:
        chain_to_resi[r["chain"]].append(pymol_resi_token(int(r["resseq"]), str(r["icode"])))

    lines: List[str] = []
    # Quote filename to handle spaces
    lines.append(f'load "{pdb_basename}"')
    lines.append("hide everything")
    lines.append("show cartoon")
    lines.append("color gray70, all")

    # highlight target chains
    for ch in target_chains:
        lines.append(f"color slate, chain {pymol_chain_token(ch)}")

    # hotspots selection
    sels: List[str] = []
    for ch, resis in chain_to_resi.items():
        resi_expr = "+".join(resis)
        sels.append(f"(chain {pymol_chain_token(ch)} and resi {resi_expr})")

    if sels:
        lines.append("select hot, " + " or ".join(sels))
        lines.append("show sticks, hot")
        lines.append("color yellow, hot")
        lines.append("set stick_radius, 0.25, hot")
        lines.append("zoom hot, 12")

    lines.append("set cartoon_transparency, 0.15")
    lines.append("bg_color white")

    with open(out_pml, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


# -----------------------------
# CLI
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Complex structure file (.pdb or .cif/.mmcif).")
    ap.add_argument("--target_chains", required=True, help="Target chain(s) to bind, e.g. 'A' or 'A,C'")
    ap.add_argument("--partner_chains", required=True, help="Partner chain(s) currently binding target, e.g. 'B' or 'B,D'")
    ap.add_argument("--cutoff", type=float, default=4.5, help="Heavy-atom contact cutoff in Angstrom. Default=4.5")
    ap.add_argument("--sigma", type=float, default=2.0, help="Score smoothing sigma in Angstrom. Default=2.0")
    ap.add_argument("--top", type=int, default=60, help="How many top interface residues to report. Default=60")
    ap.add_argument("--select", type=int, default=6, help="How many hotspots to select for RFdiffusion. Default=6")
    ap.add_argument("--mode", choices=["spread", "patch"], default="spread",
                    help="Hotspot selection mode. spread=spatially separated; patch=cluster around top residue.")
    ap.add_argument("--min_sep", type=float, default=7.0, help="Min CA-CA separation for spread mode. Default=7.0")
    ap.add_argument("--patch_radius", type=float, default=12.0, help="Patch radius for patch mode. Default=12.0")
    ap.add_argument("--binder_len", default="80-120", help="Binder length range for contig suggestion. Default=80-120")
    ap.add_argument("--outdir", default="hotspot_out", help="Output directory.")
    ap.add_argument(
        "--exclude_residues",
        default="",
        help="Comma-separated residues to exclude (set score=0 / do not select). Example: A10,A25,B50 or AA:10,AA:10A",
    )

    args = ap.parse_args()

    pdb_path = _normalize_cli_path(args.pdb)

    target_chains = [c.strip() for c in args.target_chains.split(",") if c.strip()]
    partner_chains = [c.strip() for c in args.partner_chains.split(",") if c.strip()]

    os.makedirs(args.outdir, exist_ok=True)

    exclude_set = parse_exclude_residues(args.exclude_residues) if args.exclude_residues else set()

    atoms, residue_ca, residue_name, chain_ranges = parse_structure_atoms(pdb_path)

    rows, neighbor_mode = score_interface(
        atoms, residue_ca, residue_name,
        target_chains=target_chains,
        partner_chains=partner_chains,
        cutoff=float(args.cutoff),
        sigma=float(args.sigma),
        exclude=exclude_set,
    )

    # write ranked table (kept the original column order for compatibility)
    tsv_path = os.path.join(args.outdir, "ranked_interface_residues.tsv")
    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("rank\tchain\tresseq\ticode\tresname\tscore\tn_contacts\tn_partner_res\tmin_dist\n")
        for i, r in enumerate(rows[: int(args.top)], start=1):
            f.write(
                f"{i}\t{r['chain']}\t{r['resseq']}\t{r['icode']}\t{r['resname']}\t"
                f"{r['score']:.3f}\t{r['n_contacts']}\t{r['n_partner_res']}\t{r['min_dist']:.3f}\n"
            )

    selected = select_hotspots(
        rows, residue_ca, k=int(args.select), mode=str(args.mode),
        min_sep=float(args.min_sep), patch_radius=float(args.patch_radius),
        exclude=exclude_set,
    )

    hot_str = make_rf_hotspot_string(selected)
    contig_sug = make_contig_suggestion(chain_ranges, target_chains, binder_len=str(args.binder_len))

    # write selected
    sel_path = os.path.join(args.outdir, "selected_hotspots.txt")
    with open(sel_path, "w", encoding="utf-8") as f:
        f.write("# RFdiffusion-ready\n")
        f.write(hot_str + "\n")
        if contig_sug:
            f.write(contig_sug + "\n")

    # write PyMOL script
    pml_path = os.path.join(args.outdir, "hotspots.pml")
    pdb_basename = os.path.basename(pdb_path)
    write_pymol_script(pml_path, pdb_basename, selected, target_chains)

    # print summary
    print(f"[OK] Parsed structure: {pdb_path}")
    print(f"[OK] Neighbor search: {neighbor_mode} | cutoff={args.cutoff}Å")
    print(f"[OK] Output dir: {args.outdir}")
    if exclude_set:
        # Show in a stable order
        excl_sorted = sorted(list(exclude_set), key=lambda x: (x[0], x[1], x[2]))
        excl_str = ", ".join([f"{c}{r}{i}" for (c, r, i) in excl_sorted])
        print(f"[OK] Excluded residues: {excl_str}")

    print("\nChain ranges detected:")
    for ch, (mn, mx) in sorted(chain_ranges.items()):
        print(f"  Chain {ch}: {mn}-{mx}")

    print("\nTop selected hotspots:")
    for r in selected:
        rid = residue_id_for_outputs(r["chain"], int(r["resseq"]), str(r["icode"]))
        print(f"  {rid} {r['resname']}  score={r['score']:.2f}  contacts={r['n_contacts']}  min_dist={r['min_dist']:.2f}")

    print("\nPaste into RFdiffusion:")
    print(f"  '{hot_str}'")
    if contig_sug:
        print(f"  '{contig_sug}'")

    print("\nFiles written:")
    print(f"  - {tsv_path}")
    print(f"  - {sel_path}")
    print(f"  - {pml_path}")


if __name__ == "__main__":
    main()
