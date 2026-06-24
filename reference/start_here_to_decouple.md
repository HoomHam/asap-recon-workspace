# Decoupling Plan: ASAP → ASAP + XeCS

**Created:** 2026-06-24  
**Status:** ✅ EXECUTED 2026-06-24 (see "Execution Result" below)  
**Architecture plan (full detail):** `workspace/Repo_Architecture.md`  
**Read handoffs too:** `handoff-report.md` (root), `workspace/handoff-report.md`,
`2026_XeCS_Recon/workspace/handoffs/XeCS_Handoff_2026-06-24.md`

---

## ✅ Execution Result (2026-06-24)

Done under Hooman's directive **"nothing disappears; label truly-extra `_delete`."**

- **`2026_XeCS_Recon/` created** (git init, commit `c20b0d4`, 45 tracked files):
  `recon/` operators, `pipeline/` 4D, `workspace/compare/`, `workspace/lustig/`,
  `workspace/reference/` (8 CS docs), `workspace/data/config.json`, CLAUDE.md×2 (local,
  per global ignore), AGENTS.md×2, fresh `.venv` (py3.11.15, same pins).
- **Self-tests PASS:** `recon/selftest.py`, `recon/selftest_4d.py` (run from /tmp).
- **Wiring = `.pth`, not sys.path edits** (zero file edits): `xecs_recon.pth` in both
  venvs points at `2026_XeCS_Recon/recon`. ASAP stayers resolve operators from XeCS.
- **No deletions.** Redundant ASAP originals moved to `workspace/helpers/_delete/`
  (Phases 11) — verified byte-identical to XeCS first. Remove later when XeCS is daily.
- **ASAP internal moves** (Phase 9): `codes/2023_Faraz_Recon_HH/`→`faraz/`;
  `codes/2025_*` + `2025-09-24_ACR`→`archive/codes/`. `kasap.c` KEPT in `codes/`.
- **Deviations:** (1) `cg_tune.py` → `recon/` not `compare/` (it's in the `cs_recon`
  import chain; compare/ would break `selftest_4d`). (2) Phase 11 `git rm` → move to
  `_delete/`. (3) CS reference docs kept in BOTH repos (cross-linked from ASAP docs).
- **Commits:** ASAP workspace `8a07d2d`; XeCS `c20b0d4`. Root repo untouched (no commit).
- **Still TODO (Phase 14):** instruction docs `Auto_CS_Static/4D/BART/Lustig_Recon.md`
  (XeCS) + `Auto_Faraz_Recon.md` (ASAP). Lustig `.venv_lustig` recreate in XeCS.

*(Original plan kept verbatim below for reference.)*

---

## What This Does

Splits the current monolithic workspace into two independent repos matching
Hooman's two-repo system (main repo + workspace, per project):

```
BEFORE:
  2026_ASAP_Recon/          ← Steve/Kento repo (read-only)
    workspace/              ← everything: Steve pipeline + CS code + comparison tools

AFTER:
  2026_ASAP_Recon/          ← Steve/Kento repo (read-only) — UNCHANGED
    workspace/              ← ASAP-specific only: Steve pipeline, Faraz, ASAP analysis

  2026_XeCS_Recon/          ← NEW: general Xe-129 CS operator library (git init)
    workspace/              ← CS-specific: BART, Lustig, comparison, 4D pipeline
```

Cross-repo calls via `sys.path`. Shared data via `config.json` (no copies).

---

## HARD RULES — READ BEFORE TOUCHING ANYTHING

1. **NEVER** `git commit` or `git push` in `2026_ASAP_Recon/` root directory
2. `2026_ASAP_Recon/workspace/` git is separate — commits there are fine
3. `2026_XeCS_Recon/` gets a **fresh `git init`** — not a fork, not a submodule
4. Move files with **`cp`**, not `git mv` — they cross repo boundaries
5. Delete from ASAP workspace ONLY after verifying XeCS copies are correct
6. Never touch `main.py`, `tyger_recon.py`, `convert_siemens_to_mrd.py` in root

---

## Phase 1 — Create XeCS Repo Skeleton

```bash
cd /Users/hoomham/Hooman/Work/Codes
mkdir -p 2026_XeCS_Recon/{recon,pipeline}
mkdir -p 2026_XeCS_Recon/workspace/{compare,lustig,outputs,data,instructions,reference,handoffs,archive}
cd 2026_XeCS_Recon
git init
```

Write `2026_XeCS_Recon/CLAUDE.md` with:
- Model: Fable 5
- Purpose: general Xe-129 CS recon library (asap_recon.py + cs_recon*.py are the main code)
- Two-repo system: this is the CS main repo; ASAP workspace calls into it via sys.path
- Data paths: never hardcode — always read from `workspace/data/config.json`
- Git: all commits allowed here

Write `2026_XeCS_Recon/workspace/CLAUDE.md` with:
- Git: commit/push freely (this .git is 2026_XeCS_Recon's own)
- CS comparison scripts live in `compare/`; they import from `../recon/` and `../../recon/`
- ASAP data accessed via `data/config.json` — do not copy data here
- Run with `workspace/.venv/bin/python` (created in Phase 2)

---

## Phase 2 — Create XeCS Python Environment

The existing arm64 venv (Python 3.11.15, finufft 2.5.1, sigpy 0.1.27) cannot be
moved — absolute paths are baked in. Create a fresh identical one in XeCS:

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon

# Use the same Python binary as the existing venv
/usr/bin/python3.11 -m venv workspace/.venv 2>/dev/null || \
  /opt/homebrew/Caskroom/miniforge/base/bin/python3.11 -m venv workspace/.venv 2>/dev/null || \
  python3.11 -m venv workspace/.venv

# Install exact same packages as ASAP workspace/helpers/.venv
workspace/.venv/bin/pip install \
  finufft==2.5.1 \
  sigpy==0.1.27 \
  numpy==2.4.6 \
  scipy==1.17.1 \
  matplotlib==3.10.9 \
  PyWavelets==1.9.0 \
  nibabel==5.4.2 \
  h5py==3.16.0 \
  numba==0.65.1 \
  tqdm==4.68.2 \
  pyMapVBVD==0.6.1 \
  pillow==12.2.0
```

Verify:
```bash
workspace/.venv/bin/python -c "import finufft, sigpy, numpy; print('venv OK')"
```

---

## Phase 3 — Copy CS Operator Library → XeCS/recon/

Source: `/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/helpers/recon/`

```bash
SRC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/helpers/recon
DST=/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/recon

cp $SRC/asap_recon.py      $DST/
cp $SRC/cs_recon.py        $DST/
cp $SRC/cs_recon_4d.py     $DST/
cp $SRC/binning.py         $DST/
cp $SRC/surrogates.py      $DST/
cp $SRC/selftest.py        $DST/
cp $SRC/selftest_4d.py     $DST/
```

---

## Phase 4 — Copy 4D CS Pipeline → XeCS/pipeline/

```bash
SRC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/helpers/recon
DST=/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/pipeline

cp $SRC/cine_4d.py          $DST/
cp $SRC/export_4d.py        $DST/
cp $SRC/slice_video.py      $DST/
cp $SRC/nav_movie.py        $DST/
cp $SRC/surrogate_compare.py $DST/
cp $SRC/diaphragm_check.py  $DST/
cp $SRC/kernel_check.py     $DST/
```

---

## Phase 5 — Copy CS Comparison Tools → XeCS/workspace/compare/

```bash
SRC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/helpers/recon
DST=/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/workspace/compare

cp $SRC/bart_compare.py         $DST/
cp $SRC/bartio.py               $DST/
cp $SRC/lustig_compare.py       $DST/
cp $SRC/wavelet_twoway.py       $DST/
cp $SRC/tv_threeway.py          $DST/
cp $SRC/slice_matched_compare.py $DST/
cp $SRC/z_register_compare.py   $DST/
cp $SRC/zreg_sixway_montage.py  $DST/
cp $SRC/resolution_sweep.py     $DST/
cp $SRC/cg_tune.py              $DST/
cp $SRC/metrics_v2.py           $DST/
```

---

## Phase 6 — Copy Lustig One-Shot → XeCS/workspace/lustig/

```bash
SRC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/helpers/lustig_oneshot
DST=/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/workspace/lustig

cp -r $SRC/. $DST/
# Note: .venv_lustig (torch+torchkbnufft) is inside lustig_oneshot/ — copy it too
# It uses absolute paths so must be recreated if it breaks:
#   DST/.venv_lustig/bin/pip install torch torchkbnufft
```

---

## Phase 7 — Copy CS Reference Docs → XeCS/workspace/reference/

```bash
SRC=/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/reference
DST=/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/workspace/reference

cp $SRC/CS_Implementation.md                      $DST/
cp $SRC/Dynamic_4D_CS_Implementation.md           $DST/
cp $SRC/4D_CS_Theory_Limitations.md               $DST/
cp $SRC/Physics_Notes.md                          $DST/
cp $SRC/Compressed_Sensing_Dynamic_Imaging.md     $DST/
cp $SRC/Compressed_Sensing_Dynamic_Imaging_Literature.md $DST/
cp $SRC/MR_Registration_in_k-Space.md             $DST/
cp $SRC/MR_Registration_in_k-Space_and_Motion.md  $DST/
```

These stay in ASAP workspace/reference/ — do NOT move:
- `Final_Report_CS_Comparison.md`
- `BART_Comparison.md`
- `Lustig_CS_Baseline.md`, `Lustig_CS_Tuning.md`
- `Recon_Overview_Steve.md`, `Recon_Overview_Faraz.md`
- `Recon_Comparison_StaticGas.md`
- `Tyger_Setup.md`

---

## Phase 8 — Wire the Cross-Repo Interlink

### 8a — Data config (XeCS reads ASAP data)

Write `/Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon/workspace/data/config.json`:

```json
{
  "asap_workspace": "/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace",
  "recon_io_static": "/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/v3_fov250/recon_io",
  "recon_io_dyn_025JC": "/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/v3_dyn_025JC/recon_io_dyn",
  "recon_io_dyn_027JC": "/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/v3_dyn/recon_io_dyn"
}
```

### 8b — XeCS compare scripts: read config + import operators

Add this header block to every script in `XeCS/workspace/compare/` that needs
data paths or CS operators. Insert after existing imports, before any `import
asap_recon` or hardcoded paths:

```python
import sys, json, pathlib

# XeCS repo root (2 levels up from workspace/compare/)
_XECS = pathlib.Path(__file__).parents[2]
# CS operators
sys.path.insert(0, str(_XECS / 'recon'))
# ASAP data paths
_cfg = json.loads((_XECS / 'workspace' / 'data' / 'config.json').read_text())
RECON_IO      = pathlib.Path(_cfg['recon_io_static'])
RECON_IO_DYN  = pathlib.Path(_cfg['recon_io_dyn_025JC'])
```

Do the same for scripts in `XeCS/pipeline/` that need data:

```python
import sys, json, pathlib

_XECS = pathlib.Path(__file__).parents[1]          # pipeline/ is 1 level from XeCS root
sys.path.insert(0, str(_XECS / 'recon'))
_cfg  = json.loads((_XECS / 'workspace' / 'data' / 'config.json').read_text())
```

### 8c — ASAP workspace scripts: import XeCS operators

For `workspace/helpers/recon/compare_baseline.py` and `workspace/helpers/recon/cs_montage.py`
(the two ASAP-side scripts that import from asap_recon/cs_recon), add:

```python
import sys, pathlib

# 4 parents up from workspace/helpers/recon/ reaches Codes/
_CODES = pathlib.Path(__file__).parents[4]
sys.path.insert(0, str(_CODES / '2026_XeCS_Recon' / 'recon'))
```

---

## Phase 9 — Internal ASAP Workspace Moves (git mv)

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace

# Faraz MATLAB fork: codes/ → faraz/
git mv codes/2023_Faraz_Recon_HH faraz

# Archive old code folders
mkdir -p archive/codes
git mv codes/2025_CS          archive/codes/2025_CS
git mv codes/2025_Xe129_CS    archive/codes/2025_Xe129_CS
git mv codes/2025-09-24_ACR   archive/codes/2025-09-24_ACR
git mv codes/AGENTS.md        archive/codes/AGENTS.md
rmdir codes   # should be empty now

git add -A
git commit -m "restructure: archive old codes/, move Faraz fork to faraz/"
```

---

## Phase 10 — Verify XeCS Copies Are Correct

Before deleting anything from ASAP workspace:

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon

# Self-tests must pass
workspace/.venv/bin/python recon/selftest.py
workspace/.venv/bin/python recon/selftest_4d.py

# Spot-check an operator import
workspace/.venv/bin/python -c "from asap_recon import recon; print('operator OK')"

# Spot-check a compare script can find config + operators
workspace/.venv/bin/python -c "
import sys, json, pathlib
_XECS = pathlib.Path('.')
sys.path.insert(0, str(_XECS / 'recon'))
cfg = json.loads((_XECS / 'workspace/data/config.json').read_text())
import asap_recon
print('compare imports OK')
print('recon_io:', cfg['recon_io_static'])
"
```

---

## Phase 11 — Delete Moved Files From ASAP Workspace

Only run this after Phase 10 passes.

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace

RECON=helpers/recon

# Operators → moved to XeCS/recon/
git rm $RECON/asap_recon.py $RECON/cs_recon.py $RECON/cs_recon_4d.py
git rm $RECON/binning.py $RECON/surrogates.py
git rm $RECON/selftest.py $RECON/selftest_4d.py

# 4D pipeline → moved to XeCS/pipeline/
git rm $RECON/cine_4d.py $RECON/export_4d.py $RECON/slice_video.py
git rm $RECON/nav_movie.py $RECON/surrogate_compare.py
git rm $RECON/diaphragm_check.py $RECON/kernel_check.py

# Comparison tools → moved to XeCS/workspace/compare/
git rm $RECON/bart_compare.py $RECON/bartio.py $RECON/lustig_compare.py
git rm $RECON/wavelet_twoway.py $RECON/tv_threeway.py
git rm $RECON/slice_matched_compare.py $RECON/z_register_compare.py
git rm $RECON/zreg_sixway_montage.py $RECON/resolution_sweep.py
git rm $RECON/cg_tune.py $RECON/metrics_v2.py

# Lustig one-shot → moved to XeCS/workspace/lustig/
git rm -r helpers/lustig_oneshot/

git commit -m "restructure: remove CS code moved to 2026_XeCS_Recon"
```

---

## Phase 12 — Update AGENTS.md + CLAUDE.md Files

### Files to update in ASAP workspace:

**`workspace/helpers/recon/AGENTS.md`** — remove all entries for moved files;
keep only: dump_inputs*.py, compare_baseline.py, steve_kernel_numpy.py,
faraz_montage.py, faraz_zoom_check.py, convert_calib.py, cs_montage.py.
Update pitfalls section: note XeCS repo for CS operators.

**`workspace/CLAUDE.md`** — update helpers/ folder description; remove
lustig_oneshot/ reference; add `faraz/` folder row; update active reference
docs table (remove CS docs that moved to XeCS).

**`workspace/helpers/recon/README.md`** — update file list.

### Files to create in XeCS:

**`2026_XeCS_Recon/recon/AGENTS.md`** — entry points (selftest, asap_recon,
cs_recon, cs_recon_4d), contracts (venv path, operator adjoint test),
pitfalls (DCF required for PDHG/FISTA; CG converges without it).

**`2026_XeCS_Recon/workspace/compare/AGENTS.md`** — entry points for each
comparison script, data path convention (always via config.json), z-affine
registration required before cross-pipeline comparison.

---

## Phase 13 — First XeCS Git Commit

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_XeCS_Recon

git add recon/ pipeline/
git add workspace/compare/ workspace/lustig/ workspace/reference/
git add workspace/data/config.json
git add CLAUDE.md workspace/CLAUDE.md
git add recon/AGENTS.md workspace/compare/AGENTS.md
git commit -m "initial: Xe-129 CS recon library + workspace extracted from ASAP workspace"
```

---

## Phase 14 — Write Instruction Docs (Style: Auto_Steve_Recon.md)

Reference style at: `workspace/instructions/Auto_Steve_Recon.md`
(build map → call graph → data provenance → bash run guide → comparison playbook)

### In ASAP workspace `instructions/` (write these):

| File | Entry point |
|------|-------------|
| `Auto_Faraz_Recon.md` | `workspace/faraz/spiral_human_20240227.m` (MATLAB) |

### In XeCS workspace `instructions/` (write these):

| File | Entry point |
|------|-------------|
| `Auto_CS_Static_Recon.md` | `2026_XeCS_Recon/recon/cs_recon.py` |
| `Auto_CS_4D_Recon.md` | `2026_XeCS_Recon/pipeline/cine_4d.py` |
| `Auto_BART_Recon.md` | `2026_XeCS_Recon/workspace/compare/bart_compare.py` + `~/bin/bart-src/bart` |
| `Auto_Lustig_Recon.md` | `2026_XeCS_Recon/workspace/lustig/run_lustig.py` |

---

## What STAYS in ASAP workspace/helpers/recon/ — Do NOT Move

| File | Why it stays |
|------|-------------|
| `dump_inputs.py` | Imports Steve's loaders from root repo |
| `dump_inputs_dyn.py` | Same |
| `compare_baseline.py` | Compares FINUFFT vs Steve — ASAP-anchored |
| `steve_kernel_numpy.py` | Validates Steve's GPU kernel |
| `faraz_montage.py` | Faraz figure tools |
| `faraz_zoom_check.py` | Faraz zoom bug documentation |
| `convert_calib.py` | Calibration format converter |
| `cs_montage.py` | Our CS vs Steve vs Faraz (cross-repo, ASAP side) |

---

## Cross-Repo Diagram (Final State)

```
2026_ASAP_Recon/                          2026_XeCS_Recon/
  workspace/helpers/recon/                  recon/
    compare_baseline.py ──sys.path──────>   asap_recon.py
    cs_montage.py       ──sys.path──────>   cs_recon.py
    dump_inputs.py                          cs_recon_4d.py
    steve_kernel_numpy.py                   binning.py
    faraz_*.py                              surrogates.py
                                          pipeline/
                                            cine_4d.py ──sys.path──> ../recon/
                                          workspace/compare/
                                            bart_compare.py ──sys.path──> ../../recon/
                                            wavelet_twoway.py             + config.json → ASAP data
                                            zreg_sixway_montage.py
                                          workspace/lustig/
                                            run_lustig.py
                                          workspace/data/
                                            config.json ──paths──> ASAP workspace/data/
```

---

## Verification Checklist

- [ ] `2026_XeCS_Recon/workspace/.venv/bin/python recon/selftest.py` → PASS
- [ ] `2026_XeCS_Recon/workspace/.venv/bin/python recon/selftest_4d.py` → PASS
- [ ] `compare_baseline.py` can import `asap_recon` from XeCS via sys.path
- [ ] `bart_compare.py` can read `config.json` and find recon_io/ data
- [ ] ASAP workspace `git status` clean after Phase 11 commit
- [ ] XeCS `git log` shows initial commit with all expected files
- [ ] `workspace/faraz/` exists with Faraz's MATLAB code
- [ ] `workspace/archive/codes/` exists with old library code
- [ ] `workspace/helpers/lustig_oneshot/` is gone from ASAP workspace
- [ ] `workspace/codes/` is gone from ASAP workspace
