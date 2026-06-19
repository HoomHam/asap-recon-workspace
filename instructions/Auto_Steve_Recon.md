# Auto-Steve-Recon — Build Map, Run Guide, and CS-Comparison Playbook

What this is: a **one-command pipeline** that runs Steve Kadlecek's ASAP spiral
recon on a dataset using the **Tyger cloud GPU** (the M4 Mac has no NVIDIA GPU),
then writes images + diagnostic figures. Built on top of Steve's code (read-only in
the root repo) and Kento's headless Tyger entrypoint, with a DIAPHRAGM navigator
ported in and a custom post-processor.

Two ways to invoke:
- **Agent:** say *"analyze 25JC with Steve's implementation"* → `asap-recon` agent
  (`workspace/.claude/agents/asap-recon.md`) runs the script and reports.
- **Bash yourself:** see [Run it in bash](#run-it-in-bash).

---

## 1. Component map — what is where, what calls what

```
                ┌─────────────────────────── LOCAL (M4 Mac, miniforge py3.13) ──────────────────────────┐
 you / agent ──▶│ workspace/pipeline/asap_run.py   (ORCHESTRATOR, 6 stages)                              │
                │   1 resolve_dataset()  ── globs workspace/data/xe/<datatype>/*/*<token>*/               │
                │   2 get_params()       ── workspace/pipeline/param_gui.py  (tkinter)  → params.json     │
                │   3 convert()          ── root: convert_siemens_to_mrd.convert_siemens_to_mrd()         │
                │                              .dat + traj .npy + pneumotach  →  runs/<id>/input.mrd       │
                │   4 submit_to_tyger()  ── ~/bin/tyger run exec -f recon_codespec.yml --logs ───────────┐ │
                │   5 publish()          ── tyger_deploy/plot_recon.py montages → workspace/outputs/...  │ │
                │   6 post_process       ── workspace/pipeline/post_process.py  (mat + figs)             │ │
                └────────────────────────────────────────────────────────────────────────────────────┘ │ │
                                                                                                        │ │
                ┌─────────────────────────── CLOUD (Tyger GPU node, container) ──────────────────────────┘ │
                │ ghcr.io/hoomham/xe-tyger-recon:<sha>   (built by GitHub Actions on the fork)              │
                │   ENTRYPOINT: tyger_recon.py  reconstruct_from_mrd(input.mrd → output.mrd)                │
                │     reads MRD header params + NdArrays                                                    │
                │     results.calcb()              ── b-matrix (B0/off-res + coil phase)                    │
                │     _diaphragm_navigator()       ── low-res nav loop (only if binning==DIAPHRAGM)         │
                │     results.dyn_recon()           ── binned gridding recon (GPDYN/DPDYN)                  │
                │     _write_results_to_mrd()       ── images + nav arrays → output.mrd ────────────────────┘
                └──────────────────────────────────────────────────────────────────────────────────────────┘
```

### Files

| File | Where | Role |
|------|-------|------|
| `workspace/pipeline/asap_run.py` | workspace | Orchestrator. 6 stages: resolve → params → convert → submit → publish → post_process. |
| `workspace/pipeline/param_gui.py` | workspace | tkinter param picker (MS, IS, nbins, griddx, bindt, gplb, dplb, freqfilter, **binning**). |
| `workspace/pipeline/post_process.py` | workspace | Reads `output.mrd` → `recon.mat` + `fig/` (slice videos, diaphragm.gif, navigator.gif, resp_traces.png). |
| `workspace/pipeline/recon_codespec.yml` | workspace | Tyger job spec. **Pins the container image `:<sha>`** (immutable; Tyger caches mutable tags). |
| `convert_siemens_to_mrd.py` | root (read-only) | `.dat` (mapvbvd) + traj `.npy` + pneumotach → `input.mrd`, bakes params into the MRD header. |
| `tyger_recon.py` | **fork branch** `diaphragm-recon` | Cloud entrypoint. Reads `input.mrd`, runs recon, writes `output.mrd`. Has the ported DIAPHRAGM navigator. |
| `results.py` | root (read-only) | Steve's recon: `calcb` (b-matrix), `dyn_usimg_recon` (per-interleave low-res), `dyn_recon` (binned). |
| `raw.py`, `gtypes.py`, `recon.py` | root (read-only) | Raw load/binning, types/param defaults, CUDA gridding kernel. |
| `tyger_deploy/plot_recon.py` | root (read-only) | Static gas/dissolved montage PNGs (stage 5). |
| `Dockerfile`, `requirements.txt` | root (on fork branch) | Container build. `mrd-fork` pinned to `b6b6d18` (schema must match the local writer). |
| `.github/workflows/build-image.yml` | root (on fork branch) | CI: builds amd64 image, pushes `ghcr.io/hoomham/xe-tyger-recon:diaphragm-recon` + `:<sha>`. |

### How the container gets built (so you know the "build")
1. Edit `tyger_recon.py` / `requirements.txt` / `Dockerfile` on fork branch
   `diaphragm-recon` (NOT on MEDCAP origin — root repo is push-forbidden).
2. `git push hooman diaphragm-recon` → GitHub Actions builds the amd64 image and
   pushes it to `ghcr.io/hoomham/xe-tyger-recon` with a mutable `diaphragm-recon`
   tag **and** an immutable `:<commit-sha>` tag.
3. Put the new `:<sha>` in `workspace/pipeline/recon_codespec.yml` (use the SHA tag —
   Tyger GPU nodes cache the mutable tag and won't re-pull it).
4. The ghcr package is **public** so Tyger can pull it anonymously.

---

## 2. Data provenance — what lives in the MRD files

This is the fast path for "where is X". The MRD files are the only channel between
local and cloud; everything is an NdArray discriminated by a `meta` key, plus header
`user_parameters`.

### `runs/<id>/input.mrd` (written by `convert_siemens_to_mrd`)
| meta key | array | meaning |
|----------|-------|---------|
| `gas_phase_trajectory` | (3, …) | k-space traj for gas phase (from `<seqname>_gp.npy`) |
| `dissolved_phase_trajectory` | (3, …) | dissolved-phase traj (`_dp.npy`; may be absent → gas-only) |
| `dynamic_acquisition` | (chan, samp, lines) | raw dynamic spiral FIDs |
| `reference_acquisition` | (chan, samp, lines) | breath-hold reference (absent for single-`.dat`) |
| `pneumotach` | (2, N) | row0 = time, row1 = pressure |

Header `user_parameters`: `MS, IS, nbins, gplb, dplb, freqfilter, griddx, bindt,
killpts, TR, TE, DPoff, dtdyn, dtspec, numspec, nusimg, binning`.

### `output.mrd` (written by `tyger_recon._write_results_to_mrd`)
| meta key | array | meaning |
|----------|-------|---------|
| `gas_phase_image` | (nbins, Z, Y, X) | binned gas recon (GPDYN) |
| `dissolved_phase_image` | (nbins, Z, Y, X) complex | dissolved recon (DPDYN) |
| `nav_coronal` | (nframes, z, x) | low-res navigator coronal projection (apex-up, diaphragm-down) |
| `nav_diaphragm_z` | (nframes,) | tracked diaphragm z per navigator frame (high-z edge) |
| `nav_time` | (nframes,) | mean interleave time per navigator frame |
| `nav_volume` | (ntotalilvs,) | interpolated navigator respiratory waveform (binning surrogate) |
| `nav_ilvtime` | (ntotalilvs,) | time of each interleave |
| `nav_volmeastime` | (nmeas,) | times of the accepted diaphragm measurements |

`recon.mat` mirrors all of the above (gas_phase, dissolved_phase_*, diaphragm_pos
[image-derived], nav_*, fid_signal, pneumo_*).

---

## 3. Run it in bash

```bash
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3.13   # arm64; has mrd + mapvbvd

# 0. (once) make sure you're logged in to Tyger
~/bin/tyger login /Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml

# 1. run — GUI param picker pops up (pick binning: SIGNAL / PNEUMOTACH / DIAPHRAGM)
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon
$PY workspace/pipeline/asap_run.py 25JC

# variants
$PY workspace/pipeline/asap_run.py 25JC --no-gui          # gvar defaults, no GUI
$PY workspace/pipeline/asap_run.py 25JC --params p.yml     # params from yaml
$PY workspace/pipeline/asap_run.py 25JC --slice 50         # montage slice index
```

Result → `workspace/outputs/25JC/<timestamp>/` : `output.mrd`, `recon.mat`,
`fig/{axial,coronal,sagittal,diaphragm,navigator}.gif`, `resp_traces.png`,
`tyger.log`, `manifest.json`.

### Re-process figures only (no cloud run; output.mrd already exists)
```bash
$PY workspace/pipeline/post_process.py workspace/outputs/25JC/<timestamp>
```

### Recover a finished-but-not-downloaded job (buffer TTL ~1 hr)
```bash
# outputBufferId is in runs/<id>/tyger.log
~/bin/tyger buffer read <outputBufferId> -o workspace/outputs/.../output.mrd
```

### Verify DIAPHRAGM actually ran (not silent SIGNAL fallback)
```bash
grep "binned dynamic recon" workspace/outputs/25JC/<timestamp>/tyger.log
# want: [main] binned dynamic recon with DIAPHRAGM binning...
```

---

## 4. CS-comparison playbook (for next-session "why is Steve's X vs our CS")

When comparing our CS recon (`workspace/helpers/recon/`) to Steve's, start here.
Map of "if the question is about… look at…":

| Question touches | Look at |
|------------------|---------|
| Undersampling / interleaves / nusimg | `input.mrd` header `nusimg`; `raw.py` binning; `nav_volume` (which interleave→bin) |
| Respiratory binning / why bins differ | `output.mrd` `nav_volume` + `nav_diaphragm_z`; `tyger_recon._diaphragm_navigator` |
| B0 / off-resonance / coil phase | `results.calcb` (b-matrix `self.b`); **navigator uses LOW-res b (MS=IS+4)** then full-res — confirm CS uses correct-res B0 (see handoff B0 note) |
| Gridding kernel / DCF | `recon.py` `cudarecon` (Gaussian active, KB commented); CS DCF in `helpers/recon/` |
| Trajectory / k-space | `input.mrd` `gas_phase_trajectory`; `raw.py` `traj.load` |
| Image orientation / z-sign | navigator GPDYN is `np.flip(axis=0)` vs final recon (`results.py` dyn_usimg_recon vs dyn_recon) — known gotcha |
| Steve's numbers vs ours | `recon.mat` (`gas_phase`, `nav_*`) is the comparison artifact; load and diff |

**Provenance rule for an agent:** Steve's diagnostic quantities (diaphragm track,
binning waveform) are computed **in `tyger_recon.py` on the cloud** (ported from
`main.py:calcLVcb`) and **exported into `output.mrd` / `recon.mat`**. `post_process.py`
only *displays* them — it does not recompute. So to explain a Steve-side number,
read the `nav_*` arrays first, then the algorithm in `tyger_recon._diaphragm_navigator`
/ `results.py`.

---

## 5. Should this be a graphify graph?

Short answer: **this curated doc + the AGENTS.md intent layer is the better primary
tool**; graphify is a useful *secondary* for open-ended exploration.

- For targeted "where is X / what calls Y / where does this number come from"
  questions (your navigator example), a precise hand-written map (sections 1–4
  above) + `workspace/pipeline/AGENTS.md` answers in one read, cheaply and
  authoritatively. Graphify would re-derive it fuzzily and cost more tokens.
- Graphify pays off when the question is broad/structural across many files you
  *haven't* curated ("trace every consumer of the b-matrix across both repos"), or
  when you want community/cluster views. There is already a `graphify-out/` at the
  repo root.
- Recommendation: keep **this doc as the source of truth**, point the CS agent here
  first (and to `AGENTS.md`), and only reach for graphify for exploratory
  cross-repo sweeps. If a graphify query reveals something this doc lacks, add it
  here so the next agent gets it for free.
```
