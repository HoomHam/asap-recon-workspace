---
status: DRAFT — read-only survey, nothing run or changed (2026-09-14, Opus 5 session, requested by System session system-f3 for Hooman)
purpose: document our Siemens .dat → MRD → recon path before the 2026-09-15 Oxford gas-imaging meeting (raw-data exchange, their modified ISMRMRD)
confidence tags: READ = seen in code/doc · MEASURED = checked on a real file this session · INFERRED = reasoning, not verified
sharing: internal. Contains Hooman-only paths. Strip §1 paths + §5.4 before sending outside Penn.
---

# Siemens raw .dat → MRD → Tyger recon: how our pipeline actually works

## TL;DR

1. **No "MRD1".** Nothing in the repo, canon, handoffs or vault uses that term. We write **MRD v2**,
   the yardl binary stream, through the **MEDCAP fork** `mrd-python` pinned at `b6b6d18`. We do **not**
   write ISMRMRD 1.x HDF5. If Oxford's "modified ISMRMRD" is v1/HDF5, the two formats cannot read each other at all.
2. We **don't use** `siemens_to_ismrmrd`, `ismrmrd-python`, or any XSL parameter map. The converter is one
   Python script: `pymapvbvd` reads the TWIX file, and the result is written through the `mrd` fork.
3. The file is **not a spec-style raw stream**. It has no `Acquisition` records and no `encoding`. Each
   acquisition goes in whole, as one `NdArray`, and a **private `meta` key** says what it is. The recon
   parameters travel as header `userParameters`. Stock MRD v2 has no `meta` on arrays, so its readers reject the schema.
4. A lot of the physics lives **outside the file**: which interleaves are gas vs dissolved, the Xe
   frequency, bandwidth, FOV, and trajectory units. Anyone else reading our MRD needs `raw.py` too.
5. **These are human-subject data.** The `.dat` headers and the raw folders carry PHI. Our MRD leaves
   PHI out by construction. The IRB attestation covers sending data to Tyger, **not** to Oxford.

---

## 1. Getting the raw data

| Item | What we have | Tag / source |
|---|---|---|
| Scanner | Siemens **MAGNETOM Avanto**, 1.494 T; software **syngo MR D13** (a VD-line version, not VB17) | MEASURED: ASCII header of one 2024-01 dynamic `.dat` (`SoftwareVersions`, `ManufacturersModelName`, `flMagneticFieldStrength`) |
| File container | **VD/VE multi-RAID** layout. The file starts `uint32 0, uint32 nMeas=1`, then a MeasID/FileID entry (0x39=57, 0x1A8D=6797) that matches the name `meas_MID00057_FID06797_…` | MEASURED (`xxd` first 16 bytes). Only one file checked, so VD13 across 2022–2024 is **not verified** |
| Multi-RAID handling | `mapVBVD` returns a list when there are several measurements. We take the **last** one (the image scan; earlier entries are adjustment scans) | READ `convert_siemens_to_mrd.py:45-47` |
| Nucleus / frequency | Xe-129 at 17.666 MHz, **hard-coded**, not read from the header | READ `raw.py:29` |
| How files leave the scanner | **Not documented anywhere** in the repo, canon, handoffs or vault (TWIX export? who? which media?) | OPEN — ask Hooman |
| File naming | Standard TWIX: `meas_MID<5>_FID<5>_<protocol/seqname>.dat`. Sequences are Faraz's `fa_spiral_dyn_fancy_v2_20230131`, `…_v3_20230821`, `…_v3_20240130` | MEASURED (folder listing); READ `workspace/data/xe/human/traj/seqnames.txt` |
| Session folder contents | `*.dat` (dynamic, breath-hold reference, sometimes other nuclei, e.g. a C-13 calibration), `DICOMs/`, a per-session `.txt`, a `pneumotach_<n>_<nn>_<date>_<time>` binary, `Spectra.png` | MEASURED (one 2024-01 folder) |
| Master raw store | `/Volumes/HoomHamExt/_5t_images_roundtrip/Images/<YYYY-MM-DD>/<ID>/` (124 date folders; audit: 136 session folders → 104 unique sessions) | MEASURED `ls`; READ `dyn_recon.py:66`, ledger 2026-09-12/13 |
| Laptop copy (small subset) | `workspace/data/xe/human/<date>/<ID>/` = `gtypes.gvar.basefolder` + `/human` | READ `gtypes.py:27` |
| Trajectories | `workspace/data/xe/human/traj/<seqname>_{gp,dp}.npy`, shape `(nilv·npts_full, 3)` float64: v2 = (327680, 3), v3 = (425984, 3) = 832 × 512. They come from Faraz's MATLAB `KSpaceCoor` `.mat` via a 2-line `traj/convert.py`. `seqnames.txt` rows are `<vendor> <seqname>` | MEASURED shapes; READ `traj/convert.py:6-7` |
| Trajectory units | Inferred **cycles/mm**: `raw.py` multiplies by FOV = 350 mm to get Δk units. Max \|k\| ≈ 0.108 /mm for v3 | INFERRED `raw.py:60`, FOV hard-coded `raw.py:27` |
| Missing trajectories | 20 sessions on 2022 sequences have **no trajectory file** anywhere, so they can't be reconstructed | facts **F51** |
| Backup | Ext is **not backed up**. Is `_5t_images_roundtrip` the only copy of the raw `.dat`? | OPEN — flag to Hooman |

Selecting the files (the dynamic scan and its reference):
- `convert_siemens_to_mrd._classify_dat_files`: after the optional `seqname` substring filter, **largest `.dat` = dynamic, second-largest = breath-hold reference** (`convert_siemens_to_mrd.py:11-35`).
- GUI (`main.py:390-445`): matches `seqnames.txt` entries against file names, builds the trajectory paths `traj/<seqname>_gp.npy` / `_dp.npy` (`:436-445`), and finds the pneumotach by the substring `pneumotach` (`:425`). A `rawdata.job0` match switches to a Bruker path (`:400-401`) that the MRD converter does not handle.
- Batch (`workspace/pipeline/dyn_recon.py`): skips `._*` shadow files (`:96-99`), keeps `.dat` names that match a known seqname and takes the largest (`:122-132`), and **passes `--seqname`** to convert (`:195`). It also stages a symlink directory without the AppleDouble `._*.dat` files (`:170-188`, fact F16).

## 2. The conversion

### Tools / libraries
| Library | Version | Where pinned | Role |
|---|---|---|---|
| `pyMapVBVD` | 0.6.1 locally; unpinned in container | `requirements.txt:7` | TWIX reader (`twix.image.unsorted()`, `twix.hdr.MeasYaps`) |
| `mrd-python` (**MEDCAP/mrd-fork**, `python/` subdir) | `2026.6.12` @ commit `b6b6d184f59f…` (installed from `@dev`, recorded in `direct_url.json`) | `requirements.txt:9` | MRD v2 yardl binary writer/reader |
| numpy / scipy / numba / llvmlite | 2.4.6 / 1.17.1 / 0.65.1 / 0.47.0 | `requirements.txt:1-6` | recon (container). The numba pin is why the pipeline works again (F42) |
| `siemens_to_ismrmrd`, `ismrmrd` (v1), XSL/XML parameter maps | **not used** | — | — |

### Entry points and exact commands
```bash
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python3.13   # arm64 miniforge base: the only local env with mrd + mapvbvd (F15)

# (a) standalone converter. -i is a DIRECTORY, not a .dat (convert_siemens_to_mrd.py:229-230)
$PY convert_siemens_to_mrd.py -i <session_dir> -o input.mrd \
    --seqname fa_spiral_dyn_fancy_v3_20240130 \
    --gp-traj workspace/data/xe/human/traj/fa_spiral_dyn_fancy_v3_20240130_gp.npy \
    --dp-traj workspace/data/xe/human/traj/fa_spiral_dyn_fancy_v3_20240130_dp.npy \
    --pneumotach <session_dir>/pneumotach_… \
    --binning DIAPHRAGM            # SIGNAL | PNEUMOTACH | DIAPHRAGM; also --ms --is --nbins --griddx --bindt --gplb --dplb --freqfilter

# (b) the production path: resolve → convert → Tyger → plot → post-process, one session per call
$PY workspace/pipeline/dyn_recon.py --date <YYYY-MM-DD> --id <ID> \
    --methods s,p,d --codespec workspace/pipeline/recon_codespec_40d23a4.yml

# (c) GUI: main.py "Convert + Recon" → xe_dyn_raw.mrd, then tyger_recon.py LOCALLY (needs CUDA; fails on the Mac)
```
- CLI definition: `convert_siemens_to_mrd.py:226-258`. Function: `:115-224`. `killpts` is fixed at 2 from the CLI (`:257`).
- GUI call: `main.py:246-275`. It passes the GUI-filtered `dat_files` and chooses PNEUMOTACH only when exactly one pneumotach file exists (`:253-254`).
- `dyn_recon` convert call: `dyn_recon.py:191-207`. Tyger submit (never `--ttl`, F44): `buffer create` → `buffer write -i` → `run create -f <codespec> -b …` → poll → `buffer read -o -p 4` (`dyn_recon.py:244-306`, F43).
- **Pitfall:** `dyn_recon`'s default `--codespec` is `pipeline/recon_codespec.yml` (`:48`), which still points at the broken image `d136eb1` (card C13). Always pass `recon_codespec_40d23a4.yml`.
- **Stale doc:** `reference/Tyger_Setup.md:83` shows `--input /path/to/data.dat` (the code wants a directory) and gives the basefolder as `data/tygerinputs` (it is `data/xe`, `gtypes.py:27`).

### What `_read_twix` pulls from the header (`convert_siemens_to_mrd.py:37-88`)
| MRD user param | TWIX key (MeasYaps) | Unit written | Line |
|---|---|---|---|
| `TR` | `alTR[0]` | s (µs × 1e-6) | :57, :173 |
| `TE` | `alTE[0]` | s | :61, :174 |
| `DPoff` | `sWipMemBlock.adFree[2]` | ppm (dissolved-phase offset; `raw.py:110`) | :65, :175 |
| `dtdyn` | `sRXSPEC.alDwellTime[1]`, falling back to `[0]` | s (ns × 1e-9) | :69-74, :176 |
| `dtspec` | `sWipMemBlock.alFree[12]` | s (µs × 1e-6) | :84, :177 |
| `numspec` | `alFree[10] × alFree[11]` (raw product; the ÷20 happens recon-side) | count | :77-80, :182 |

The `sWipMemBlock` indices are **specific to Faraz's `fa_spiral_dyn` WIP sequence**. They don't carry over to another group's sequence. Every read is wrapped in `try/except: pass`, so a missing key **silently becomes 0** (`:56-86`).
The raw data comes from `twix.image.unsorted()` with `flagRemoveOS=False`, as complex64 `(samples, channels, lines)`. A 2-D result gets a singleton channel axis (`:48-52`). The mapVBVD end-of-file UserWarning is harmless (F33).

## 3. What our MRD file contains (the "MRD1" question)

### Format identity
- The protocol is `Mrd`: `header: Header?`, then `data: stream<StreamItem>`. It's written with `mrd.BinaryMrdWriter` (`convert_siemens_to_mrd.py:222-224`).
  The binary reader checks magic bytes, then the format version, then compares the **embedded schema JSON** string against its own. Any difference raises `RuntimeError: Invalid schema` (fork `mrd/_binary.py:78-89`; READ).
  So even two commits of the **same fork** can't read each other's files (F3 is why `b6b6d18` is pinned).
- **Fork vs stock MRD v2** (READ: fork schema from the installed `mrd/protocols.py`; stock from `ismrmrd/mrd` main `model/mrd_protocol.yml` and `mrd_intermediate.yml`):

| | MEDCAP fork @ b6b6d18 (ours) | stock ismrmrd/mrd main |
|---|---|---|
| Generic arrays | `ndArrayUint16…ndArrayComplexDouble` → `NdArray<T>{head: NdArrayHeader, data, meta: map<string, vector<ArrayMetaValue{string\|int64\|float64}>>}` | `arrayComplexFloat: Array<complexfloat>`, no header and **no meta** |
| Other StreamItem cases | acquisition, waveformUint32, image*, acquisitionBucket, reconData, imageArray, Pulseq (`pulseqDefinitions, blocks, rfEvent, arbitraryGradient, trapezoidalGradient, adcEvent, shape`) | adds `acquisitionPrototype`; Pulseq cases named `pulseq*` |
| Header.encoding | `vector<EncodingType>` (may be empty) | same concept |

  So a stock MRD v2 tool can't open our file, and an ISMRMRD v1 (HDF5 + XML header) tool certainly can't.

### Stream items we write (`convert_siemens_to_mrd.py:198-220`), in order
| meta key (value `'1'`) | StreamItem | dtype | shape | Source |
|---|---|---|---|---|
| `gas_phase_trajectory` | NdArrayDouble | float64 | (nilv·npts_full, 3) | `<seq>_gp.npy` (:201-204) |
| `dissolved_phase_trajectory` | NdArrayDouble | float64 | same | `<seq>_dp.npy`, optional (:205-208) |
| `pneumotach` | NdArrayDouble | float64 | (2, N): row 0 = t [s] from the first packet, row 1 = pressure | `_parse_pneumotach` (:91-112): 54-byte packets with magic `0xA6 0x20`, ms timestamp at bytes 2-5, P = −20 + 90·u16/65535 at bytes 32-33 |
| `reference_acquisition` | NdArrayComplexFloat | complex64 | **(channels, samples, lines)** | second-largest `.dat` (:214-217) |
| `dynamic_acquisition` | NdArrayComplexFloat | complex64 | (channels, samples, lines) | largest `.dat` (:218-220) |

- `NdArray.head` is left at defaults: no FOV, position or orientation (only `data=` and `meta=` are passed).
- The dynamic array matches Steve's direct mapVBVD read bit-for-bit (F49, n=1 session).

### Header fields we fill in (`convert_siemens_to_mrd.py:155-195`)
| Field | Value |
|---|---|
| `sequence_parameters.t_r` / `t_e` | ms (for information only; the recon ignores them) (:158-161) |
| `acquisition_system_information.receiver_channels`, `system_vendor='Siemens'` | :163-166 |
| `user_parameters.user_parameter_double` | `TR, TE, DPoff, dtdyn, dtspec` (from TWIX) + `griddx, bindt` (recon) (:172-180) |
| `user_parameters.user_parameter_long` | `numspec` + `nusimg=32` (hard-coded) + `killpts` + `MS, IS, nbins, gplb, dplb, freqfilter` (:181-191) |
| `user_parameters.user_parameter_string` | `binning` (:192-194) |
| **Not set** | `version`, `encoding` (empty: no encodedSpace/reconSpace/trajectory type/FOV/matrix), `experimental_conditions` (no H1 or Xe resonance frequency), `measurement_information`, and `subject_information` / `study_information` (on purpose: PHI) |

The userParameters are the **only way** recon settings reach the headless cloud job (`pipeline/AGENTS.md:27-28`).

### Xenon-specific content and what is *not* in the file
| Physics | Where it actually lives |
|---|---|
| Gas vs dissolved interleaves | **Not in the MRD.** `raw.py:263-287` infers the pattern (gas-only / gas-dissolved / gas-dissolved-dissolved) from peaks in the FFT of the FID-intensity sequence, then slices every 2nd or 3rd interleave |
| Spectroscopy prefix | `numspec = numspec_raw / 20 · nsmpperusimg / npts` ("TEMP KLUGE … 20 MAGIC COOKIE"), `raw.py:255-257` |
| TR rounding fix | 22.3 ms → 22.26 ms, `raw.py:247-248` |
| Per-channel noise normalization | Gaussian fit to the real-part histogram, `raw.py:249-253` |
| Samples per interleave | inferred from the \|k\|² periodicity of the trajectory, not the header (`raw.py:63-69`, F18); `killpts` trimmed on the sample axis (`tyger_recon.py:206-208`, F37) |
| Spiral BW, FOV, Xe frequency, spectral BW | hard-coded `raw.py:26-29` (1/10 µs, 350 mm, 17.666 MHz, 1/60 µs) |
| RBC/TP split | computed recon-side (`results.py` 298-327); unreliable in 60/86 split sessions (F47) |
| xyz trajectory calibration | not used in this path. Trajectories are Faraz's precomputed `KSpaceCoor`, not measured xyz calibrations |

## 4. Downstream

| Stage | Code | Reads | Writes / lands |
|---|---|---|---|
| Cloud recon | `tyger_recon.py` (container entrypoint, `Dockerfile:36`) | `BinaryMrdReader`: tells items apart by type + meta key (`:157-172`); userParams → `gvar` + meta dict (`:175-191`); gp trajectory required (`:195-196`), dp optional since `472fbc9` (`:197-198`) | `output.mrd` with the **input header passed through** (`:139-141`) |
| Raw → recon objects | `raw.load_from_arr(..., 'mrd_siemens', meta)` `raw.py:201-246` | transposes to internal `(samples, channels, lines)` (`:207-216`) | — |
| Binning | `tyger_recon.py:217-234`; DIAPHRAGM navigator `:23-115` (z-flip = a real change from Steve, F49) | — | — |
| GPU gridding | `recon.py` `cudarecon` (Gaussian kernel; KB commented out, F17) via `results.py` | — | — |
| Output items | `_write_results_to_mrd` `tyger_recon.py:117-141` | — | `gas_phase_image` float32 (nbins, Z, Y, X); `dissolved_phase_image` complex64, `rbc_tp_separated` `'1'` = aRBC + i·aTP, `'0'` = unsplit (`:125-131`); DIAPHRAGM only: `nav_coronal, nav_diaphragm_z, nav_time, nav_volume, nav_ilvtime, nav_volmeastime` (`:108-115`, `:132-138`) |
| Montages | `tyger_deploy/plot_recon.py:41-48` | output.mrd | `output_gp.png` / `_dp.png` |
| Post-process | `workspace/pipeline/post_process.py:55-69` (images), `:127-134` (re-reads `dynamic_acquisition` + `pneumotach` from input.mrd) | both MRDs | `recon.mat` (`:451-452`), `signal_pneumo.npz`, `fig/*.gif`, `resp_traces.png` |
| Output location (batch) | `dyn_recon.py:67` | — | `/Volumes/HoomHamExt/AIkill_Dynamic/<date>_<ID>/{s,p,d}/` (output.mrd, input.mrd, tyger.log with run/buffer ids, codespec.yml, pngs, recon.mat, fig/). Scratch: `/Volumes/HoomHamExt/Work/Codes/2026_ASAP_Recon/pipeline_runs/` (`:52`) |
| Output location (single / GUI) | `asap_run.py` → `workspace/outputs/<dataset>/<timestamp>/`; `main.py:251-252` → `<session_dir>/xe_dyn_{raw,recon}.mrd` | — | — |
| Paths that skip MRD | Faraz MATLAB (`mapVBVD.m`); `helpers/recon/dump_inputs_dyn.py` (copies `_read_twix`/`_parse_pneumotach` to avoid the `mrd` dependency); XeCS CS pipeline starts from those dumps | `.dat` directly | — |

**Container / environment:** `nvidia/cuda:12.6.2-devel-ubuntu24.04` (devel image needed for libnvvm; Ubuntu 24.04 because the fork needs Python 3.12), venv `/opt/venv` (`Dockerfile:8-32`). Only `gtypes/raw/recon/results/tyger_recon.py` are copied in (`:34`). The fork's GitHub Actions build a linux/amd64 image on each push to `diaphragm-recon`, tagged `ghcr.io/hoomham/xe-tyger-recon:{diaphragm-recon,<sha>}` (`.github/workflows/build-image.yml`, fork branch). Codespecs must pin `:<sha>` (F2). Current: `recon_codespec_40d23a4.yml`, cluster `tep-centralus-1`, nodePool `gpunp`, timeout 43200 s. `tyger_deploy/recon_codespec.yml` still points at `ghcr.io/medcap/xe-tyger-recon:latest` (Kento's upstream).

## 5. Pitfalls, retractions, open questions, and the Oxford exchange

### 5.1 Canon rows that bear on this path
| Fact | Status | Point |
|---|---|---|
| F3 | VALID | mrd fork must be pinned `b6b6d18`, otherwise `RuntimeError: Invalid schema` |
| F15 | VALID | only miniforge base has `mrd`; the arm64 `helpers/.venv` does not |
| F16 | VALID | `._*.dat` AppleDouble shadows on exfat get picked as the reference; filter `._*` |
| F18 / F37 | VALID | npts from \|k\|² periodicity; nusimg = 32 hard-coded; killpts trims samples only |
| F33 | VALID | mapVBVD EOF warning is harmless |
| F46 | VALID | 2024-01 8-ch sessions are gas-only **acquisitions** (the ledger's "RBC/TP params absent" was wrong, RETRACTED 2026-09-13) |
| F47 / F48 | VALID | stored aRBC/aTP not trustworthy in split sessions; use the whitened total-DP magnitude |
| F49 | VALID | MRD raw == mapVBVD raw bit-for-bit (n=1); the DIAPHRAGM z-flip is our change, not an MRD artefact |
| F51 | VALID | 20 sessions on 2022 sequences have no trajectory |
| F13 | RETRACTED (2026-07-10) | "missing `_dp.npy` → gas-only, not fatal" was false then. It **has been true since `472fbc9`** (`raw.py:62`, `tyger_recon.py:197`). `pipeline/AGENTS.md:36-37` states it without that history |
| F1 / F41 → F42, F40 → F43 | RETRACTED | numba pin fixed the pipeline; download speed varies |

### 5.2 Code-level traps
- Every TWIX header read swallows errors, so a different sequence gives zeros silently (`convert_siemens_to_mrd.py:56-86`). Guard before trusting TR/TE/DPoff on non-Faraz data.
- "Largest `.dat` = dynamic" fails if the folder holds bigger unrelated scans. Always pass `--seqname`.
- `seq.t_r` is in ms while userParam `TR` is in s. Only the userParam is used.
- The meta-key convention and array axis order `(channels, samples, lines)` are ours alone. Nothing in the file documents them.

### 5.3 Open questions (for Hooman)
1. Who exports `.dat` off the Avanto, and how (TWIX via the scanner console? USB/network?). Undocumented.
2. Is `_5t_images_roundtrip` on Ext the **only** copy of the raw `.dat`? (Ext is not backed up.)
3. Pneumotach binary: whose device/firmware is it (Spinhance?), and is there a written packet spec? The parser is reverse-engineered offsets.
4. Software version across the cohort: only one file checked (syngo MR D13).
5. Trajectory units and origin: confirm cycles/mm and which Faraz script generated `KSpaceCoor`.
6. What did "MRD1" mean in the request? v1 ISMRMRD, or "our MRD"?

### 5.4 Swapping data with a group on a different MRD variant
**Pin down first (Tuesday):** their exact format (ISMRMRD v1 HDF5 or MRD v2; library + commit; schema/XSD), how they encode gas vs dissolved (separate datasets? `idx.contrast`? flags?), how they carry trajectories (per-readout `traj` in `Acquisition` or separate), their scanner software (VB/VE/XA), and whether they want raw k-space or reconstructed images.

**What breaks if we send our `input.mrd` as-is:**
- A v1 reader can't parse it (different container). A stock v2 reader rejects the schema. A reader on the same fork at a different commit also rejects it (F3).
- Even when parsed, it lacks per-readout headers, `encoding`, and resonance frequency. Gas/dissolved assignment and scaling hacks live in `raw.py`.

**Realistic options** (a decision for Hooman, not done):
- **(a)** Export spec-style data: one `Acquisition` per interleave, with `trajectory`, `idx.contrast` (0 = gas, 1 = dissolved), `sample_time_us`, `center_frequency`, plus `encoding.trajectory=spiral`. This means a new writer, and the interleave pattern must be applied at export time (logic from `raw.py:263-287`).
- **(b)** Keep our format and write a small translator into their variant, run in its **own env or container**. Don't install their `ismrmrd`/`mrd` into miniforge base, or the `b6b6d18` pin used by our converter will break (F3/F15).
- **(c)** Neutral bundle: `.npy`/HDF5 arrays + a JSON sidecar holding §3's tables. Lowest effort, and easy to review.
- In every case send a **README** covering axis order, units, meta keys, and the hard-coded constants in `raw.py:26-29`. Do a round trip on **one phantom dataset** before any human data.

**Human-subject data (read before sharing anything):**
- `.dat` headers carry PHI: name, patient ID, date of birth (`regulatory/IRB_Cloud_Recon_Letter_2026-06-26.md:87-91`).
- Session folders also hold `DICOMs/`, per-session text files, and time-stamped pneumotach file names.
- Our MRD omits all of that by construction (letter `:93-108`). But that attestation covers **Tyger cloud recon**. Sending data to an outside group (Oxford) needs its own IRB/DUA clearance.
- **Never send raw `.dat` or session folders.**
- **Canon conflict (A5):** root `CLAUDE.md` and `workspace/CLAUDE.md` both say "phantom data … No PHI", but this cohort pipeline runs on human sessions. Hooman to decide which statement gets corrected.
