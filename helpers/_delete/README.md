# _delete/ — Redundant Duplicates (safe to remove)

**Created 2026-06-24 by the ASAP→XeCS decouple.**

These files were moved OUT of the ASAP workspace because their canonical home is now
the new `2026_XeCS_Recon` repo. They are kept here (not deleted) per the "nothing
disappears" rule — every byte was verified identical to the XeCS copy before the move.

## What's here

| Path | Canonical home now |
|------|--------------------|
| `recon/` (25 files: operators, 4D pipeline, CS comparison) | `2026_XeCS_Recon/recon/`, `2026_XeCS_Recon/pipeline/`, `2026_XeCS_Recon/workspace/compare/` |
| `lustig_oneshot/` | `2026_XeCS_Recon/workspace/lustig/` |

## Why it's safe to delete this folder

- The ASAP-side scripts that still need these operators
  (`helpers/recon/compare_baseline.py`, `cs_montage.py`, `faraz_montage.py`,
  `faraz_zoom_check.py`) now import them from XeCS via an `xecs_recon.pth` in
  `helpers/.venv/lib/python3.11/site-packages/` (points at `2026_XeCS_Recon/recon`).
  Verified: imports resolve to the XeCS copies, all stayers compile.
- The XeCS self-tests pass (`recon/selftest.py`, `recon/selftest_4d.py`).

## Note on `lustig_oneshot/.venv_lustig`
That venv has absolute paths baked in; moving it broke its shebangs. The working
Lustig flow is in `2026_XeCS_Recon/workspace/lustig/` (recreate its `.venv_lustig`
with `torch torchkbnufft` if needed). This copy is preserved for reference only.

**To finish the cleanup later:** confirm XeCS is in daily use, then `git rm -r` this
folder from the workspace.
