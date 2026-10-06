"""Which halo is in which Arrow file.

The file names lie: within each pair the two are swapped, so the file named
G11 holds the halo at G1.2's AHF centre and G31/G32 likewise, 850+ kpc apart.
Prints the dark-matter centroid of each file against every AHF centre in the
data readme so the assignment can be checked rather than trusted.

  python3 whichcentre.py *.arrow
"""
import sys, os, numpy as np, pyarrow as pa, pyarrow.ipc as ipc
class JS(pa.ExtensionType):
    def __init__(self): super().__init__(pa.string(), "JuliaLang.Symbol")
    def __arrow_ext_serialize__(self): return b""
    @classmethod
    def __arrow_ext_deserialize__(cls, st, s): return cls()
try: pa.register_extension_type(JS())
except pa.lib.ArrowKeyError: pass
CEN = {"G1.1":(47.3092422,48.8026602,50.0027344),
       "G2.1":(48.8151992,46.7045156,53.6050664),
       "G3.1":(46.3911250,50.7461367,47.9353516),
       "G1.2":(46.7936523,49.0552930,49.8811797),
       "G2.2":(48.7157578,47.0623281,53.3371875),
       "G3.2":(46.7533242,50.3214297,47.7923984)}
for p in sys.argv[1:]:
    with pa.memory_map(p,"r") as s:
        t = ipc.open_file(s).read_all()
    pt = t["ptype"].combine_chunks().storage.to_numpy(zero_copy_only=False)
    dm = pt == "dm"
    xyz = np.stack([t[f"Coordinates{i}"].to_numpy(zero_copy_only=False)
                    for i in (1,2,3)],axis=1)[dm]
    med = np.median(xyz,axis=0)
    best = sorted(CEN, key=lambda k: np.linalg.norm(med-np.array(CEN[k])))
    d0 = np.linalg.norm(med-np.array(CEN[best[0]]))*1000/0.677
    d1 = np.linalg.norm(med-np.array(CEN[best[1]]))*1000/0.677
    print(f"{os.path.basename(p)[7:-25]:22s} median {np.round(med,4)}"
          f"  -> nearest {best[0]} at {d0:7.1f} kpc"
          f"   (next {best[1]} at {d1:8.1f} kpc)")
