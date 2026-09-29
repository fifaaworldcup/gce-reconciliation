# Getting the Fermi data: what you run, what you send me

The split: **Fermitools only has to touch the data five times.** Those five
commands turn 5 GB of raw photons into two small FITS files. Send me those two
and I can do the template fitting, the likelihood, TS, ΔlnL and residuals in
pure Python here — no Fermitools needed on my side.

You are the hands. Paste me whatever breaks.

---

## Step 1 — conda + Fermitools (your Mac, your Terminal, ~20 min)

```bash
# miniforge, the conda that works cleanly on Apple Silicon
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"

# do NOT pin the python version. current fermitools needs 3.11 or 3.12 and
# pinning 3.9 makes the solve fail outright:
#   fermitools 2.4.0/2.5.1 -> python_abi 3.11
#   fermitools 2.5.3       -> python_abi 3.12
#   python 3.9             -> conflicts with both
conda create -y -n fermi -c conda-forge -c fermi fermitools
conda activate fermi

# the four tools we actually need
which gtselect gtmktime gtbin gtltcube gtexpcube2
python -c "import gt_apps; print('fermitools OK')"
```

**We do not need Fermipy for this.** The pure-Python fitter only wants a counts
cube and an exposure cube, and those come from the `gt*` command-line tools.
Skipping Fermipy also sidesteps a version fight — Fermipy 1.2 predates these
Python versions. Fermipy is only needed later, for the Level-2 iterative
pipeline.

Note this installs **Fermitools 2.5.3**, not the 2.2 in Project Conventions §7.
2.2 is not what conda serves any more. Flag the version drift to Afeefa and
record whatever version you actually get, since Conventions §12 wants the
software version logged with each run.

If the solve still fails on Apple Silicon, force x86 under Rosetta:

```bash
CONDA_SUBDIR=osx-64 conda create -y -n fermi -c conda-forge -c fermi fermitools
conda activate fermi && conda config --env --set subdir osx-64
```

## Step 2 — download the photons (FSSC web form, ~30 min of waiting)

https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi

| Field | Value |
|---|---|
| Object/coordinates | `0, 0` |
| Coordinate system | **Galactic** |
| Search radius | `30` degrees |
| Observation dates | `239557417` to `760320000` |
| Time system | **MET**, not Gregorian |
| Energy range | `500, 50000` MeV |
| **Spacecraft data** | **tick this — nothing works without it** |

**Use MET, not Gregorian.** Fermi science data begins at MET 239557417, which
is 2008-08-04 *15:43:37*, not midnight. Entering the date as `2008-08-04
00:00:00` is 15.7 hours too early and the form rejects the whole query with
"Start time occurs before data start MET". MET sidesteps the whole issue.

The two numbers above are our Conventions window: MET 239557417 is the first
moment of science data (MJD 54682.65), and MET 760320000 is MJD 60710
(2025-02-04). Note that is **16.5 years, not the 16.7 the Conventions claim** -
worth raising with Afeefa, though it changes nothing about the selection.

It emails you a list of `..._PH##.fits` files plus one `..._SC00.fits`. Download
them all into one folder. Several GB.

```bash
ls *_PH*.fits > events.txt      # gtselect wants a list file
```

## Step 3 — the five Fermitools commands

Run these from inside the `fermi data` folder.

```bash
conda activate fermi

# file lists. use *PH*, not *_PH* - the underscore may not be there, and a
# glob that matches nothing makes zsh abort the whole command.
ls *PH*.fits > events.txt
wc -l events.txt          # expect 7
SC=$(ls *SC00.fits); echo "spacecraft file: $SC"

# 1. event selection (Conventions section 1)
#    keep the FULL 0.5-50 GeV range here so we don't have to redownload
#    when we widen the analysis later.
gtselect infile=@events.txt outfile=gap1_sel.fits \
  ra=266.404988 dec=-28.936178 rad=30 \
  tmin=239557417 tmax=760320000 emin=500 emax=50000 \
  zmax=90 evclass=128 evtype=3

# 2. good time intervals
gtmktime scfile="$SC" \
  filter="DATA_QUAL>0 && LAT_CONFIG==1 && ABS(ROCK_ANGLE)<52" \
  roicut=no evfile=gap1_sel.fits outfile=gap1_gti.fits

# 3. counts cube, 1-10 GeV  <-- send me this
#    binned 1000-10000 in 8 bins, NOT 500-50000. binning the full range in 8
#    puts edges at 889 / 1581 / ... so 1 GeV and 10 GeV fall mid-bin and the
#    band Afeefa asked for can't be cut out cleanly. 8 bins across one decade
#    also matches Di Mauro's binning.
gtbin evfile=gap1_gti.fits scfile=NONE outfile=gap1_ccube.fits \
  algorithm=CCUBE nxpix=400 nypix=400 binsz=0.1 \
  coordsys=GAL xref=0 yref=0 axisrot=0 proj=CAR \
  ebinalg=LOG emin=1000 emax=10000 enumbins=8

# 4. livetime cube (slow - 1 to 3 hours. start it and walk away)
#    zmax MUST be repeated here. gtltcube does not inherit it from gtselect,
#    and a mismatch silently biases the exposure.
gtltcube evfile=gap1_gti.fits scfile="$SC" \
  outfile=gap1_ltcube.fits dcostheta=0.025 binsz=1 zmax=90

# 5. exposure cube  <-- send me this too
#    cmap=none, NOT cmap=gap1_ccube.fits. passing a cmap makes gtexpcube2
#    inherit that file's geometry and ignore the nxpix/nypix below, so we'd
#    get a 400x400 exposure map when we asked for 440x440.
gtexpcube2 infile=gap1_ltcube.fits cmap=none outfile=gap1_expcube.fits \
  irfs=P8R3_SOURCE_V3 evtype=3 \
  nxpix=440 nypix=440 binsz=0.1 coordsys=GAL \
  xref=0 yref=0 axisrot=0 proj=CAR \
  emin=1000 emax=10000 enumbins=8
```

`ra=266.404988 dec=-28.936178` is galactic (0,0) in J2000 — gtselect wants
celestial coordinates even though our ROI is galactic.

The exposure map is deliberately wider than the counts cube (440 vs 400, i.e.
44° vs 40°) so the PSF can pull flux in from just outside the ROI.

### Sanity checks before you send anything

```bash
python -c "
from astropy.io import fits
for f in ('gap1_ccube.fits','gap1_expcube.fits'):
    h=fits.open(f); print(f); h.info()
    d=h[0].data; print('  shape',d.shape,'sum %.4g'%d.sum(),'\n')
"
```

Expect the counts cube at `(8, 400, 400)` with a few hundred thousand to a few
million counts, and the exposure cube at `(8, 440, 440)` or `(9, 440, 440)`
with values around 10^10–10^11 cm² s. If the counts cube is empty or the
exposure is all zeros, something upstream failed — send me the terminal output
rather than the files.

## Step 4 — send me two files

```
gap1_ccube.fits     ~5 MB    the data
gap1_expcube.fits   ~50 MB   how long we looked, per direction and energy
```

Drop them in the chat. That's all I need — **not** the raw photons, not the
livetime cube.

## Also grab, whenever convenient

```bash
wget https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits
wget https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/iso_P8R3_SOURCE_V3_v1.txt
```

Plus the FL16Y catalog FITS from wherever the collaboration is hosting it —
that one I don't have a URL for, worth asking Afeefa.

---

## What I do with those two files

A binned Poisson likelihood in pure Python:

```
predicted counts  =  exposure x PSF-convolved( sum of N_k x template_k )
lnL               =  sum over bins of  n*ln(mu) - mu
```

Maximise over the normalisations, once with HESTIA as the GCE and once with
Coleman, from the identical background fit. That gives TS, best-fit norms,
ΔlnL, ΔAIC and residual maps — Afeefa's step 5, in full.

**What this is not:** it is not Fermipy, and it is not the Di Mauro Level-2
iterative source-finding. It is a clean Level-1 template comparison. That is
exactly what "prove the pipeline works before we scale up" asks for, and it
has a real advantage — when the Fermipy version does run later, we have an
independent implementation to check it against. Two codebases agreeing is
worth a lot more than one codebase running.

Level 2 still needs Fermitools, and still has to run on your Mac.
