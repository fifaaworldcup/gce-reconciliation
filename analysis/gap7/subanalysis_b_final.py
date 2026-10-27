from pathlib import Path
import numpy as np
import pandas as pd

INPUT = Path("outputs/gap7/covariance_input.csv")
OUTPUT = Path("outputs/gap7")
OUTPUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT)
df = df[
    df["admissible"].astype(str).str.strip().str.lower().eq("true")
    & df["status"].astype(str).str.startswith("OK")
    & df["gap"].eq(5)
    & df["mask"].eq("B")
    & df["arm"].eq("neutral")
    & df["energy_bin"].eq("1-10GeV")
].copy()

df["S"] = pd.to_numeric(df["dAIC_vs_best_gnfw"], errors="coerce")
df = df.dropna(subset=["S"])

def decompose(data, name):
    counts = data.groupby(["diffuse_model", "template"]).size()
    nb = data["diffuse_model"].nunique()
    nt = data["template"].nunique()

    assert len(counts) == nb * nt and counts.eq(1).all()

    grand = data["S"].mean()
    bg = data.groupby("diffuse_model")["S"].mean()
    template = data.groupby("template")["S"].mean()

    ss_bg = nt * ((bg - grand) ** 2).sum()
    ss_template = nb * ((template - grand) ** 2).sum()
    ss_total = ((data["S"] - grand) ** 2).sum()

    components = {
        "diffuse_background": ss_bg,
        "template": ss_template,
        "remainder": ss_total - ss_bg - ss_template,
    }

    return [
        {
            "analysis": name,
            "component": component,
            "sum_of_squares": value,
            "fraction_total": value / ss_total,
        }
        for component, value in components.items()
    ]

assert len(df) == 960
assert df["diffuse_model"].nunique() == 24

results = decompose(df, "all_gap5_neutral")

for name, group in df.groupby("template_class"):
    results.extend(decompose(group, name))

results = pd.DataFrame(results)
results.to_csv(
    OUTPUT / "subanalysis_b_all_classes_decomposition.csv",
    index=False,
)

print("\nVARIANCE DECOMPOSITION")
print(results[["analysis", "component", "fraction_total"]].to_string(index=False))

# Exploratory halo-level bootstrap for HESTIA
hestia = df[df["template_class"] == "dark_matter_simulated"].copy()
halos = hestia["halo"].dropna().unique()
assert len(halos) == 6
assert len(hestia) == 576

def background_fraction(data):
    grand = data["S"].mean()
    means = data.groupby("diffuse_model")["S"].mean()
    counts = data.groupby("diffuse_model")["S"].size()
    ss_bg = (counts * (means - grand) ** 2).sum()
    ss_total = ((data["S"] - grand) ** 2).sum()
    return ss_bg / ss_total

rng = np.random.default_rng(42)
samples = []

for _ in range(200):
    selected = rng.choice(halos, size=len(halos), replace=True)
    sample = pd.concat(
        [hestia[hestia["halo"] == halo] for halo in selected],
        ignore_index=True,
    )
    samples.append(background_fraction(sample))

samples = np.array(samples)
pd.DataFrame({"background_fraction": samples}).to_csv(
    OUTPUT / "subanalysis_b_bootstrap.csv",
    index=False,
)

print("\nHESTIA HALO BOOTSTRAP")
print("Original:", round(background_fraction(hestia), 4))
print("Median:", round(np.median(samples), 4))
print("95% percentile interval:",
      np.round(np.percentile(samples, [2.5, 97.5]), 4))
print("\nDone. Results saved in outputs/gap7/")
