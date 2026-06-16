# Handoff Report — Tyger GPU Recon Setup
**Date:** 2026-06-12  
**Branch:** `dev` (MEDCAP/asap_recon)  
**Session focus:** Setting up Steve/Kento's ASAP recon to run on Tyger cloud GPU with Hooman's own data.

> **Note:** This handoff covers only the Tyger/main-repo work. Hooman's personal CS recon work lives separately in `workspace/` and has its own context — do NOT mix them. If you need workspace context, ask Hooman to direct you there.

---

## What Was Done This Session

### 1. Switched to `dev` branch
```bash
git fetch origin && git checkout dev && git pull
```
The `dev` branch is Kento's Tyger port of Steve's ASAP recon. New files vs `main`:
- `Dockerfile` — container definition for cloud GPU
- `tyger_recon.py` — Docker entrypoint (runs on Tyger, requires NVIDIA GPU)
- `convert_siemens_to_mrd.py` — local tool: Siemens .dat → MRD format
- `tyger_deploy/recon_codespec.yml` — Tyger job spec (cluster, image, buffers)
- `tyger_deploy/plot_recon.py` — (added late session) plots gas/dissolved phase from output MRD

### 2. Installed tyger CLI
Binary at `~/bin/tyger`. Source: `~/Downloads/tyger_darwin_arm64.tar.gz`.  
Required `xattr -d com.apple.quarantine` to bypass macOS Gatekeeper.  
Full install steps: `workspace/reference/Tyger_Setup.md` §2.

### 3. Logged into Tyger
```bash
tyger login -f /Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml
tyger login status  # → spinhance.tyger.cloud, role: owner
```
Cluster: `tep-centralus-1`. NodePool: `gpunp`. Image: `ghcr.io/medcap/xe-tyger-recon:latest`.

### 4. Installed Python dependencies (miniforge base, arm64)
```bash
/opt/homebrew/Caskroom/miniforge/base/bin/pip install Pillow pymapvbvd numba
/opt/homebrew/Caskroom/miniforge/base/bin/pip install "mrd-python @ git+https://github.com/MEDCAP/mrd-fork.git@dev#subdirectory=python"
```
Use `python` from miniforge (`/opt/homebrew/Caskroom/miniforge/base/bin/python3`, Python 3.13, arm64).

### 5. Ran first Tyger recon job (run 254 — succeeded)
Input: `workspace/data/tygerinputs/xe_dyn_raw.mrd` (pre-existing MRD)  
Output: `workspace/data/tygerinputs/myfirstrecon.mrd`  
The download timed out mid-transfer; recon itself succeeded. Re-downloaded via:
```bash
tyger buffer read <outputBufferId> > output.mrd  # buffer TTL = 1 hr from job end
```

### 6. Fixed `main.py` for local GUI use
The GUI was broken for Hooman's machine. Two bugs fixed (NOT committed — forbidden in this repo):

**Bug 1** — missing `/` in `ID_callback` (lines 396/398):
```python
# Before (broken):
g_dir.trajdirname = g.basefolder + datatype + '/traj/'
# After (fixed):
g_dir.trajdirname = g.basefolder + '/' + datatype + '/traj/'
```
Same fix on line 398 for `datadirname`.

**Bug 2** — `listdir` not filtering hidden files or non-directories. All menu `listdir` calls now use:
```python
[f for f in listdir(path) if not f.startswith('.')]
```

### 7. Reorganized data folder
GUI expects `basefolder/datatype/date/subjectID/`. Created `human/` level:
```
workspace/data/tygerinputs/
  human/
    2024-12-09/
      027JC/         ← patient data
    traj/            ← trajectory files (seqnames.txt, *.npy)
  xe_dyn_raw.mrd     ← input MRD (at root, not inside human/)
  myfirstrecon.mrd   ← output from run 254
```
`basefolder` is set in `gtypes.py` line 27:
```python
basefolder = '/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/tygerinputs'
```

---

## How to Run a New Tyger Recon Job

### Step 0 — Check Tyger login (token may expire)
```bash
tyger login status
# If expired:
tyger login -f /Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml
```

### Step 1 — Convert new .dat file to MRD (if starting from raw Siemens data)
```bash
python convert_siemens_to_mrd.py --input /path/to/scan.dat --output workspace/data/tygerinputs/xe_dyn_raw.mrd
```

### Step 2 — Submit Tyger job
```bash
cat workspace/data/tygerinputs/xe_dyn_raw.mrd | tyger run exec -f tyger_deploy/recon_codespec.yml --logs > workspace/data/tygerinputs/output.mrd
```

### Step 3 — If download times out
```bash
tyger run list          # find runId
tyger run show <runId>  # confirm status=Succeeded, get outputBufferId
tyger buffer read <outputBufferId> > workspace/data/tygerinputs/output.mrd
# WARNING: buffer expires 1 hour after job finishes — download promptly
```

### Step 4 — Plot results
```bash
python tyger_deploy/plot_recon.py workspace/data/tygerinputs/output.mrd
# Saves output_gp.png (gas phase) and output_dp.png (dissolved phase) alongside the MRD
```

---

## Next Session Tasks

### Task 1 — Logical input/output folder structure
Current: input and output MRDs sit at `tygerinputs/` root alongside subject data. Messy.  
Proposed (discuss with Steve/Kento first):
```
workspace/data/
  input/
    human/
      2024-12-09/027JC/
    traj/
  output/
    <runId or date>/
      output.mrd, output_gp.png, output_dp.png
```

### Task 2 — Make main.py fixes survive git pulls
Current bug fixes to `main.py` are NOT committed (forbidden — main repo is upstream MEDCAP mirror).  
Every `git pull` that touches `main.py` will wipe them.

Options (discuss with Kento/Steve):
- **Option A:** Submit the fixes as a PR to MEDCAP/asap_recon `dev` branch (clean, permanent)
- **Option B:** Keep a patch file in `workspace/` and re-apply after each pull
- **Option C:** Move `basefolder` and directory filtering into `gtypes.py` config so it's isolated

For Option C: the `basefolder` default in `gtypes.py` line 27 is already Hooman-specific. The `.DS_Store` filter and directory-only filter should be proposed upstream.

---

## Critical Context for a Blind Agent

**Repo rules (read CLAUDE.md):**
- NEVER `git commit` or `git push` in `/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/`
- `git pull` allowed
- All Hooman's personal work is in `workspace/` (separate git repo)

**This session's work is in the main repo `dev` branch — not workspace.**

**Full setup reference:** `workspace/reference/Tyger_Setup.md`  
**Memory (persists across sessions):** `~/.claude/projects/-Users-hoomham-Hooman-Work-Codes-2026-ASAP-Recon/memory/tyger_setup.md`

**main.py is currently fixed on disk but the fixes are not committed.** A `git pull` that modifies `main.py` will break the GUI again. The two fixes needed are documented in §6 above and in `workspace/reference/Tyger_Setup.md` §7.

---

## Suggested Skills

- `/handoff update <note>` — amend this doc mid-session
- Standard `git pull` to sync `dev` before starting work
- Check `tyger login status` before any job submission
