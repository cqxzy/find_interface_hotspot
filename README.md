# find_interface_hotspot
interface-hotspot finder for RFdiffusion PPI binder design using gemmi + distance-based contact scoring with residue-type weighting and optional exclusions.


# RFdiffusion Interface Hotspot Finder (PPI)

A lightweight, single-file Python script to identify and rank interface “hotspot” residues from a target–partner complex for **RFdiffusion PPI binder design**.

It uses **gemmi** for robust structure parsing (**PDB / mmCIF**, including insertion codes), computes distance-based contact scores, applies simple **residue-type weighting** (aromatic/hydrophobic bonus; small/structural penalty), and optionally **excludes user-specified residues** (e.g., active sites, glycosylation sites).

---

## Installation

This script only depends on **gemmi** and **numpy**.

```bash
pip install gemmi numpy
````

> Optional acceleration: if `scipy` is available, the script will use a KD-tree implementation automatically; otherwise it falls back to a pure numpy/grid search (no extra install required).

---

## Usage

### Typical (same CLI style as the original script)

```bash
python find_interface_hotspots.py \
  --pdb input/fold_2026_02_05_19_50dele1hri_2_model_1.cif \
  --target_chains A \
  --partner_chains B \
  --select 6 \
  --mode spread \
  --outdir output/hotspot_out
```

### Exclude specific residues

Exclude residues by chain + residue number (optionally with insertion code):

```bash
python find_interface_hotspots.py \
  --pdb input/complex.cif \
  --target_chains A \
  --partner_chains B \
  --select 6 \
  --mode spread \
  --outdir output/hotspot_out \
  --exclude_residues A10,A25,B50
```

You may also use a more explicit syntax (useful when chain IDs are multi-character in mmCIF):

* `CHAIN:RESSEQ` or `CHAIN:RESSEQICODE`
  Example: `AA:10`, `AA:10A`

```bash
python find_interface_hotspots.py \
  --pdb input/complex.cif \
  --target_chains AA \
  --partner_chains BB \
  --select 6 \
  --mode spread \
  --outdir output/hotspot_out \
  --exclude_residues AA:10,AA:25,BB:50A
```

### Windows path with spaces

Use quotes:

```bash
python find_interface_hotspots.py --pdb "C:\Users\...\file.cif" --target_chains A --partner_chains B --select 6 --mode spread --outdir output\hotspot_out
```

> Note: `r"..."` is a Python string literal, not a standard command-line syntax. Quoted paths are recommended.

---

## What the script does

1. Loads **PDB/mmCIF** with `gemmi` (handles insertion codes).
2. Builds interface contacts between `--target_chains` and `--partner_chains`.
3. Scores each residue by a geometric contact function and multiplies by a residue-type weight.
4. Optionally removes residues listed in `--exclude_residues`.
5. Selects top hotspots using:

   * `--mode top`: pick the highest scoring residues
   * `--mode spread`: pick high scorers while spatially spreading them across the interface (to avoid clustering)
6. Writes outputs to `--outdir`.

---

## Scoring

### Geometric contact score

Contacts are scored based on inter-atomic distances using a smooth decay (Gaussian-like):

* Higher score for closer contacts
* Diminishing contributions with distance

### Residue-type weights (biological heuristic)

The final residue score is:

**final_score = geometric_score × residue_weight**

Weights:

* **Trp (W), Tyr (Y), Phe (F)**: `1.2` (aromatic hotspots)
* **Leu (L), Ile (I), Met (M)**: `1.1` (hydrophobic core)
* **Gly (G), Pro (P), Ala (A)**: `0.9` (often structural / lower binding contribution)
* Others: `1.0`

---

## Outputs

The script writes a small set of files to `--outdir` (names may vary slightly by your version, but typically include):

### 1) Hotspot list (RFdiffusion-friendly tokens)

A plain text line or file containing residues formatted like:

* `A100` (chain + residue number)
* If insertion code is present, it will be preserved in a safe form where applicable.

These tokens are intended for RFdiffusion hotspot inputs (e.g., `ppi.hotspot_res=[A54,A97,...]`).

### 2) Ranked table (debug/inspection)

A TSV/CSV-style table with per-residue information such as:

* chain
* residue number (resseq)
* insertion code (icode)
* residue name / one-letter code
* geometric score
* residue weight
* final score
* selection flag (whether it was chosen as a hotspot)

### 3) PyMOL selection snippet

A convenience selector expression you can paste into PyMOL to highlight hotspots on the structure.

---

## Notes / Caveats

### 1) RFdiffusion chain ID assumptions

RFdiffusion hotspot tokens typically assume **single-character chain IDs** (e.g., `A100`).

* Internally, the script safely supports multi-character chain IDs (common in mmCIF).
* If your structure uses multi-character chain IDs, the script can still score and select hotspots, but it may refuse or warn when exporting RFdiffusion-style compact tokens due to ambiguity.

**Tip:** If needed, rename chains to single letters in a preprocessing step before running RFdiffusion.

### 2) Insertion codes (icode)

Insertion codes are handled explicitly as part of the residue identity.

* Example: residue `100A` is distinct from `100`.

### 3) Non-standard residues

Non-standard amino acids are mapped to 1-letter codes using gemmi’s residue tables when possible.
Unknown residues are labeled as `X` and use default weight `1.0`.

### 4) Performance

* With `scipy` installed, KD-tree acceleration may be used automatically.
* Without it, a grid/numpy fallback is used (slower but dependency-free).

### 5) Excluding residues

`--exclude_residues` is enforced before hotspot selection.
Use it to remove:

* catalytic residues
* glycosylation sites
* known off-limits regions
* residues involved in other essential interactions

This work is done at [Yang Lab](https://jieyang-lab.com/) at [UVA](https://www.virginia.edu/) school of medicine, [Department of Biochemistry-Molecular Genetics](https://med.virginia.edu/bmg/), under the supervision of [Prof.Jie Yang](https://med.virginia.edu/faculty/faculty-listing/wfw7nc/)
