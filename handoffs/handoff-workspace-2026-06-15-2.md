# ASAP Recon — Session Handoff

**Date:** 2026-06-12 (updated in place; previous phase spec archived as `workspace/handoffs/handoff-2026-06-11.md`)
**Phase closed:** Steve-vs-Faraz comparison — FINALIZED (do not reopen) — **and CS build steps 1–2 done** (static single-bin CS working, λ sweep produced, λ verdict pending).
**Next phase:** metric fixes + λ freeze → **CS-vs-CS comparison (MATLAB Lustig + BART)** → temporal 4D.

---

## Project nature (read before anything)

Benign engineering: comparison/reimplementation of MRI image-reconstruction
code on **physical phantom data** (ACR phantom). No human data, no PHI.
Hooman's model preference for this project: **Fable 5**. See root `CLAUDE.md`
Session Notes.

## Where everything stands

| Artifact | Path | Contains |
|----------|------|----------|
| CS layer (working) | `workspace/helpers/recon/cs_recon.py` | L1-wavelet (FISTA) + TV (PDHG) on our finufft operator, DCF-preconditioned; λ-as-threshold parameterization |
| CS montage tool | `workspace/helpers/recon/cs_montage.py` | CS vs Faraz-corrected vs Steve-equiv, shared [0,1] color axis; `--t-rel` selects λ |
| CS theory + lessons | `workspace/reference/CS_Implementation.md` | Deep doc: objective, both λ-parameterization failures, conditioning trap, priors, metric blind spots |
| Node for helpers/recon | `workspace/helpers/recon/AGENTS.md` | Contracts + pitfalls — **read before touching recon code** |
| Sweep numbers | `workspace/data/v3_fov250/recon_io/cs_sweep_metrics.json` | maxeig, t_ref, all per-(reg, t) metrics |
| Sweep/verdict figures | same folder: `cs_sweep_sheet.png`, `montage_cs_t0.003.png`, `montage_cs_t0.01.png` | For Hooman's λ verdict |
| Comparison truth (closed phase) | `workspace/reference/Recon_Comparison_StaticGas.md` | Steve-vs-Faraz, estimator-objectives capstone |
| Physics primer | `workspace/reference/Physics_Notes.md` | §12 = CS motivation |
| Prior handoffs | `workspace/handoffs/handoff-2026-06-{09,10,11}.md` | History; 06-11 has the original CS spec |

**Headline results this session:**
1. CS works: wavelet t = 0.003–0.01 beats **both** bars simultaneously — SNR 32.8–47.2 (> Steve 28.7) and lowfreq-CV 0.084–0.085 (< Faraz 0.093). The handoff thesis confirmed on the fully-sampled case. TV never matched it (best lfCV 0.124).
2. Two instructive failures, both measured, both in `CS_Implementation.md` §3–4: λ in objective units = silent no-op; unweighted AᴴA stalls gradient solvers ~1000× (density spread 1.8e6) while CG cuts through. Both produced plausible-looking images — the eye-over-metrics rule held twice.
3. DCF returned as a *preconditioner* (w = 1/|AAᴴ1| in the data term) — exactly the D3 per-sample-weight slot; bin weights will multiply into the same W.
4. Implementation is **ours assembled from primitives**: FINUFFT operators (our validated wrappers), sigpy only for wavelet/FD linops + L1 prox + FISTA/PDHG solvers. NOT sigpy's packaged MRI apps, NOT BART. Matters for the comparison goals below.

## Decisions this session

- **D4 (plan change, Hooman):** standalone undersampling experiment SHELVED. It returns naturally with diaphragm binning in the 4D phase. Do not run it as a separate test.
- **D5:** λ is parameterized as soft-threshold in coefficient units relative to p99 of |W·x_cg|; never sweep λ in objective units (contract in `AGENTS.md`).
- **D6:** gradient solvers always get the DCF-weighted operator; recompute density weights per sampling mask when undersampling eventually returns.
- λ choice between t0.003 and t0.01: **Hooman's visual verdict, still pending** — montages rendered for exactly this.

## The roadmap (Hooman's ordering, 2026-06-12)

### Step 1 — Fix metrics, freeze λ
Known-broken metrics (`cg_tune.metrics`): SNR inflated when priors zero the
background (σ_bg collapses — t0.1 wavelet scored SNR 68 on a visibly worse
image); extents threshold-fragile (CG reads 250 mm because noise clears the
mask threshold). Fix candidates: noise σ from a fixed corner ROI instead of
adaptive background mask; extents from a fixed absolute threshold calibrated
on the ACR truth; add an edge-sharpness metric (resolution-insert profile)
so SNR can't be won by blurring. Then Hooman picks λ from the montages and
it gets frozen in `AGENTS.md`.

### Step 2 — CS-vs-CS comparison (the new core goal)
Compare our CS against two external implementations, same raw data:

1. **MATLAB Lustig (SparseMRI toolbox, nonlinear-CG)** — deliberately
   included *even though partially circular*: this is the CS Hooman used
   before this project, so it benchmarks **where he was standing before vs
   what was built here**. Personal-baseline comparison, not an independence
   test. (Expect adaptation work: SparseMRI is 2D-Cartesian-oriented;
   3D spiral needs its NUFFT hooks or slice-wise treatment — scope decision
   for that session.)
2. **BART `pics` (C-based)** — the independent benchmark from the
   Uecker/Lustig ecosystem. Handles 3D non-Cartesian natively (`-t`
   trajectory input, scaled to its [-N/2, N/2] convention — ours is
   grid-index, conversion needed), ℓ1-wavelet (`-l1` + `-R W:…`) and TV.
   Also the on-ramp for 4D: `pics` has first-class temporal dims and
   temporal-TV/low-rank regularizers. Install via homebrew.

Caveat recorded: sigpy itself is Lustig-lab (Frank Ong), so a sigpy-app
comparison would share lineage with our build — that's why BART is the
independence test, and the MATLAB comparison is framed as personal-baseline.

Comparison protocol: same acq/trajectory inputs (`recon_io/`), align with
`orient_to_match`, fixed metrics from step 1, montage + Hooman's eye.

### Step 3 — Temporal 4D
Diaphragm binning → bins as W_b weight vectors (architecture unchanged,
decision D3) + temporal regularizer (temporal TV and/or low-rank across
bins). **This is where undersampled recon actually gets dealt with.** Needs
dynamic data; b-map output stage (port Steve's `calcb`) becomes relevant
before multi-coil.

## Hard constraints carried forward
- Run everything with `workspace/helpers/.venv/bin/python` (finufft 2.5.1, sigpy 0.1.27, arm64).
- NEVER git commit/push in repo root; all git in `workspace/` (`HoomHam/asap-recon-workspace`) — everything through commit `9d00c94` is pushed.
- Never modify Faraz's code folders; scanner data stays in gitignored `workspace/data/`.
- `selftest.py` must pass after any operator change (operator untouched this session).
- Read `workspace/helpers/recon/AGENTS.md` pitfalls before touching recon code — λ no-op trap, conditioning trap, `cs_montage.py`↔`cs_sweep_metrics.json` coupling.

## Open loops (not blocking, don't lose)
1. **λ verdict** — Hooman, from `montage_cs_t0.003.png` vs `montage_cs_t0.01.png`.
2. **Tell Faraz** about the ×1.205 `resizing` zoom bug (numbers in `zoom_check_metrics.json`; his script `~/Hooman/Work/Analysis/2025-09-24_ACR/spiral_gpt_ACR_20250915.m:445-452`).
3. **Steve's rebuttal** to `workspace/for_steve.md` — pending; answer from comparison doc + benchmarks.
4. **b-map stage** (Steve's `calcb` port) — needed before multi-coil / 4D.
5. Colab certification of `steve_kernel_numpy` vs GPU output — optional.
6. Parked solver alternatives: ADMM-with-inner-CG (statistically cleaner than DCF-weighted objective); undecimated wavelets if block texture survives the chosen λ.
7. Untracked stray in workspace: `archive/fable_5_asap_06-12-2026.txt` — not from this session's work; Hooman to keep or delete.

## Suggested skills
- `/handoff update` as comparison milestones land
- `/code-review` on `cs_recon.py` before the BART comparison session (fresh eyes on the normalization math)

## Quick start for the cold session
```bash
cd workspace/helpers/recon
../.venv/bin/python selftest.py                                        # sanity (~1 min)
# read AGENTS.md (this folder) + reference/CS_Implementation.md §3-4 & §8
# step 1: fix cg_tune.metrics (SNR corner-ROI, fixed-threshold extents), rerun
#         cs_recon.py, get Hooman's lambda verdict, freeze it
# step 2: brew install bart; wire trajectory conversion (grid-index -> BART units)
```
