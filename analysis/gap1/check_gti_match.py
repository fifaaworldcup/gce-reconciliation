"""
Decide whether a livetime cube built from one event file can be reused with
another.

gtltcube's output depends on the spacecraft file and on the good-time
intervals in the event file, and on nothing else - not the energies, not the
event count. So a cube built from a 5-50 GeV download is valid for a
500 MeV-500 GeV download of the same time range, and rebuilding it would cost
40 minutes for a byte-identical result.

"Should be identical" is not the same as "is identical" though, and a silently
mismatched livetime cube biases every exposure and therefore every
normalisation, without failing. So compare the GTIs properly and refuse on any
difference rather than eyeballing the start and stop times.

    python check_gti_match.py old_gti.fits new_gti.fits
"""

import sys

import numpy as np
from astropy.io import fits


def gti_of(path):
    with fits.open(path) as hd:
        g = None
        for e in hd[1:]:
            if e.name.upper() == "GTI":
                g = e.data
                break
        if g is None:
            raise SystemExit(f"{path} has no GTI extension")
        start = np.asarray(g["START"], float)
        stop = np.asarray(g["STOP"], float)
        h0 = hd[0].header
        h1 = hd[1].header
        meta = {}
        for k in ("TSTART", "TSTOP", "ONTIME", "LIVETIME"):
            meta[k] = h1.get(k, h0.get(k))
    return start, stop, meta


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    a, b = sys.argv[1], sys.argv[2]
    sa, ea, ma = gti_of(a)
    sb, eb, mb = gti_of(b)

    print(f"  {a}: {len(sa):,} intervals, "
          f"{sa.min():.1f} - {ea.max():.1f} MET, "
          f"live {np.sum(ea - sa):,.1f} s")
    print(f"  {b}: {len(sb):,} intervals, "
          f"{sb.min():.1f} - {eb.max():.1f} MET, "
          f"live {np.sum(eb - sb):,.1f} s")
    print()

    ok = True
    if len(sa) != len(sb):
        print(f"  DIFFER: {len(sa):,} vs {len(sb):,} intervals")
        ok = False
    else:
        # exact equality is the right test: both come from the same gtmktime
        # filter over the same spacecraft file, so any difference at all means
        # something about the selection changed
        ds = np.abs(sa - sb).max()
        de = np.abs(ea - eb).max()
        if ds > 0 or de > 0:
            print(f"  DIFFER: interval edges move by up to "
                  f"{max(ds, de):.6g} s")
            ok = False
        else:
            print("  interval tables are identical")

    for k in ("TSTART", "TSTOP", "ONTIME", "LIVETIME"):
        va, vb = ma.get(k), mb.get(k)
        if va is None or vb is None:
            continue
        same = abs(float(va) - float(vb)) < 1e-6
        print(f"  {k:9} {float(va):>18,.3f}  vs {float(vb):>18,.3f}  "
              f"{'ok' if same else 'DIFFER'}")
        if not same:
            ok = False

    print()
    if ok:
        print("REUSABLE - the livetime cube built from the first file is valid "
              "for the second")
    else:
        print("NOT REUSABLE - rebuild the livetime cube with gtltcube")
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
