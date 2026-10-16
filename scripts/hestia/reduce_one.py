"""Reduce one HESTIA Arrow file to what the J-factor needs.

Keeps dark matter inside R_KEEP of the AHF centre, and stars inside 20 kpc so
the disc plane can be defined from their angular momentum. Coordinates go from
absolute box Mpc/h to kpc centred on the halo.
"""
import os, sys
import numpy as np, pyarrow as pa, pyarrow.ipc as ipc

H = 0.677
R_KEEP = 60.0        # kpc; the ROI needs 15, this leaves room for the LOS
R_STARS = 20.0

CENTRES = {   # Mpc/h, from the readme
 "G11": (4.73092422e+01, 4.88026602e+01, 5.00027344e+01),
 "G21": (4.88151992e+01, 4.67045156e+01, 5.36050664e+01),
 "G31": (4.63911250e+01, 5.07461367e+01, 4.79353516e+01),
 "G12": (4.67936523e+01, 4.90552930e+01, 4.98811797e+01),
 "G22": (4.87157578e+01, 4.70623281e+01, 5.33371875e+01),
 "G32": (4.67533242e+01, 5.03214297e+01, 4.77923984e+01),
}

class JuliaSymbol(pa.ExtensionType):
    def __init__(self): super().__init__(pa.string(), "JuliaLang.Symbol")
    def __arrow_ext_serialize__(self): return b""
    @classmethod
    def __arrow_ext_deserialize__(cls, st, s): return cls()
try: pa.register_extension_type(JuliaSymbol())
except pa.lib.ArrowKeyError: pass

path = sys.argv[1]
tag = [k for k in CENTRES if f"_{k}_" in os.path.basename(path)][0]
with pa.memory_map(path, "r") as src:
    t = ipc.open_file(src).read_all()
pt = t["ptype"].combine_chunks().storage.to_numpy(zero_copy_only=False)
xyz = np.stack([t[f"Coordinates{i}"].to_numpy(zero_copy_only=False)
                for i in (1, 2, 3)], axis=1).astype(np.float64)
mass = t["Masses"].to_numpy(zero_copy_only=False).astype(np.float64)

c = np.array(CENTRES[tag])
xyz = (xyz - c) * 1000.0 / H          # Mpc/h -> kpc, centred
r = np.linalg.norm(xyz, axis=1)

dm = (pt == "dm") & (r < R_KEEP)
st = (pt == "stars") & (r < R_STARS)
vel = np.stack([t[f"Velocities{i}"].to_numpy(zero_copy_only=False)
                for i in (1, 2, 3)], axis=1).astype(np.float64)[st]

np.save(f"reduced/{tag}_dm_pos.npy", xyz[dm].astype(np.float32))
np.save(f"reduced/{tag}_dm_mass.npy", mass[dm].astype(np.float32))
np.save(f"reduced/{tag}_star_pos.npy", xyz[st].astype(np.float32))
np.save(f"reduced/{tag}_star_vel.npy", vel.astype(np.float32))
np.save(f"reduced/{tag}_star_mass.npy", mass[st].astype(np.float32))
print(f"{tag}: {dm.sum():,} dm within {R_KEEP:g} kpc, "
      f"{st.sum():,} stars within {R_STARS:g} kpc, "
      f"dm mass {mass[dm].sum():.4g} (code units), "
      f"r_max kept {r[dm].max():.1f} kpc")
