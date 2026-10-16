"""Reduce one HESTIA Arrow file to what the J-factor needs.

Memory matters: these files hold ~17M particles, so columns are read one at a
time as float32 and masked immediately rather than stacked as float64.

The halo is identified by matching the dark-matter median position to the AHF
centres in the readme, NOT by the filename: the files named G11/G12 and
G31/G32 each carry the other member of their pair, 850+ kpc away.

Disc orientation comes from the stellar inertia tensor rather than angular
momentum, so velocities never need loading.
"""
import os, sys, gc, numpy as np, pyarrow as pa, pyarrow.ipc as ipc
H, R_KEEP, R_STARS = 0.677, 60.0, 20.0
CEN = {"G1.1":(47.3092422,48.8026602,50.0027344),
       "G2.1":(48.8151992,46.7045156,53.6050664),
       "G3.1":(46.3911250,50.7461367,47.9353516),
       "G1.2":(46.7936523,49.0552930,49.8811797),
       "G2.2":(48.7157578,47.0623281,53.3371875),
       "G3.2":(46.7533242,50.3214297,47.7923984)}
class JS(pa.ExtensionType):
    def __init__(self): super().__init__(pa.string(), "JuliaLang.Symbol")
    def __arrow_ext_serialize__(self): return b""
    @classmethod
    def __arrow_ext_deserialize__(cls, st, s): return cls()
try: pa.register_extension_type(JS())
except pa.lib.ArrowKeyError: pass

p = sys.argv[1]
with pa.memory_map(p, "r") as s:
    t = ipc.open_file(s).read_all()
pt = t["ptype"].combine_chunks().storage.to_numpy(zero_copy_only=False)
is_dm = pt == "dm"
is_st = pt == "stars"
del pt; gc.collect()

cols = []
for i in (1, 2, 3):
    c = t[f"Coordinates{i}"].to_numpy(zero_copy_only=False).astype(np.float32)
    cols.append(c)
med = np.array([np.median(c[is_dm]) for c in cols], dtype=np.float64)
tag = min(CEN, key=lambda k: np.linalg.norm(med - np.array(CEN[k])))
off = float(np.linalg.norm(med - np.array(CEN[tag])) * 1000 / H)
cen = np.array(CEN[tag], dtype=np.float32)

for i in range(3):
    cols[i] = (cols[i] - cen[i]) * np.float32(1000.0 / H)
r2 = cols[0] ** 2 + cols[1] ** 2 + cols[2] ** 2
dm = is_dm & (r2 < R_KEEP ** 2)
st = is_st & (r2 < R_STARS ** 2)
del r2, is_dm, is_st; gc.collect()

mass = t["Masses"].to_numpy(zero_copy_only=False).astype(np.float32)
o = "reduced/" + tag.replace(".", "")
np.save(o + "_dm_pos.npy", np.stack([c[dm] for c in cols], axis=1))
np.save(o + "_dm_mass.npy", mass[dm])
np.save(o + "_star_pos.npy", np.stack([c[st] for c in cols], axis=1))
np.save(o + "_star_mass.npy", mass[st])
print(f"{os.path.basename(p)[14:-25]:12s} -> {tag}  (centre off {off:5.1f} kpc)"
      f"   dm {int(dm.sum()):>9,}  stars {int(st.sum()):>9,}"
      f"   Mdm<60kpc {mass[dm].sum()*1e10:.3g} Msol")
