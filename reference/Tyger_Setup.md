# Tyger Recon Setup — Hooman's Machine

Setup done 2026-06-12. Covers: pulling dev branch, installing tyger CLI + Python deps, fixing main.py for local data, running a Tyger job.

---

## 1. Pull the dev branch

```bash
cd /Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon
git fetch origin
git checkout dev
git pull
```

---

## 2. Install tyger CLI

Binary lives at: `/Users/hoomham/Downloads/tyger_darwin_arm64.tar.gz`

```bash
tar -xzf ~/Downloads/tyger_darwin_arm64.tar.gz -C /tmp/
xattr -d com.apple.quarantine /tmp/tyger   # remove Gatekeeper block
cp /tmp/tyger ~/bin/tyger
chmod +x ~/bin/tyger
tyger --version
```

> If `~/bin` is not in PATH, add `export PATH="$HOME/bin:$PATH"` to `~/.zshrc`.

---

## 3. Install Python dependencies

Into miniforge base (arm64, Python 3.13):

```bash
/opt/homebrew/Caskroom/miniforge/base/bin/pip install Pillow pymapvbvd
/opt/homebrew/Caskroom/miniforge/base/bin/pip install "mrd-python @ git+https://github.com/MEDCAP/mrd-fork.git@dev#subdirectory=python"
```

> `numba` also needed for `main.py` (local GUI). Install same way if missing:
> `pip install numba`

---

## 4. Log in to Tyger

```bash
tyger login -f /Users/hoomham/Hooman/Work/Spinhance/Tyger/LOGIN_FILE.yml
tyger login status   # verify: should show spinhance.tyger.cloud, role owner
```

---

## 5. Data folder structure

GUI (`main.py`) expects:
```
basefolder/
  <datatype>/        ← e.g. "human"
    <date>/          ← e.g. "2024-12-09"
      <subjectID>/
    traj/
      seqnames.txt
      *.npy
```

Current basefolder (set in `gtypes.py`):
```
/Users/hoomham/Hooman/Work/Codes/2026_ASAP_Recon/workspace/data/tygerinputs
```

Add new subjects by creating `human/<date>/<ID>/` inside `tygerinputs/`.

---

## 6. Run a Tyger recon job

**Step 1 — Convert .dat → .mrd** (if starting from raw Siemens data):
```bash
python convert_siemens_to_mrd.py --input /path/to/data.dat --output xe_dyn_raw.mrd
```

**Step 2 — Submit to Tyger and get results**:
```bash
cat xe_dyn_raw.mrd | tyger run exec -f tyger_deploy/recon_codespec.yml --logs > output.mrd
```

> The codespec uses cluster `tep-centralus-1`, nodePool `gpunp`, image `ghcr.io/medcap/xe-tyger-recon:latest`.
> If download times out (large output), check run status and re-download:
> ```bash
> tyger run show <runId>
> tyger buffer read <outputBufferId> > output.mrd
> ```
> Buffer TTL is 1 hour from job completion — download promptly.

---

## 7. Local GUI (main.py)

Useful for browsing data and triggering conversions. Run:
```bash
python main.py
```

### Bugs fixed in dev branch (2026-06-12)
- Missing `/` between `basefolder` and `datatype` in `ID_callback` (lines 396/398)
- `listdir` calls not filtering `.DS_Store` and other hidden/non-directory entries

These fixes are in `main.py`. If the repo is reset from upstream, re-apply them.
