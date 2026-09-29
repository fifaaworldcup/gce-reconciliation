# audit of the 12 extreme IEMs picked for Family 2
#
# their stated procedure:
#   1. pick 6 IEMs from Table 2 by hand
#   2. min-max normalise every column of all 80
#   3. add the model with the LARGEST SUM OF DISTANCES to the current set
#   4. repeat to 12
#
# step 3 is the problem. "largest sum of distances" is not a spread objective.
# it rewards being far from the average of the set, which outliers achieve
# trivially - and two outliers sitting next to each other BOTH score highly,
# so the method can and does pick near-duplicates. the objective that
# guarantees coverage is maximin (farthest-point): add the model whose
# distance to its NEAREST already-chosen model is largest.
#
# their own pitfall note says "selecting extreme IEMs that are too similar to
# one another could create bias". this checks whether that happened.

import itertools
import numpy as np

# transcribed from the Gap 5 table image
COLS = ["zL", "D0", "delta", "vA", "dvc_dz", "ap1", "ap2", "ae1", "ae2", "Np", "Ne"]
CATS = ["ISRF", "gas", "SNSe", "Bfield"]

ROWS = [
    # num  zL   D0   delta vA    dvc  ap1   ap2   ae1  ae2   Np   Ne    ISRF   gas      SNSe      Bfield
    (5,  10.0, 10.0, 0.33, 32.2,   0, 1.70, 2.39, 1.6, 2.33, 1.0, 2.66, "1.4", "0,5,2", "Pul/SNR", "200040050"),
    (7,  10.0,  8.0, 0.33,  0.0,   0, 1.40, 1.80, 1.4, 2.30, 1.3, 3.33, "1.4", "0,5,2", "Pul/SNR", "200040050"),
    (17,  5.6,  4.85,0.40, 24.0,   1, 2.00, 2.25, 1.6, 2.30, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (19,  6.0,  6.5, 0.33, 30.0,   0, 2.04, 2.41, 1.6, 2.43, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (22,  5.5,  5.5, 0.37, 30.0, 2.5, 2.00, 2.38, 1.6, 2.43, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (25,  5.7,  3.9, 0.45, 25.7,   6, 1.99, 2.36, 1.6, 2.43, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (27,  6.0,  3.1, 0.50, 23.0,   9, 2.02, 2.38, 1.6, 2.43, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (30,  3.0,  2.67,0.40, 22.0,   3, 2.08, 2.41, 1.6, 2.43, 5.8, 2.00, "0.7", "9,4,1", "Pul/Pul", "025200010"),
    (51,  6.0,  2.0, 0.33,  0.0,  60, 1.60, 2.30, 1.6, 2.26, 1.5, 5.90, "1.4", "9,5,1", "Pul/SNR", "200030050"),
    (52,  6.0,  2.0, 0.33,  0.0, 100, 1.60, 2.30, 1.6, 2.26, 1.5, 5.90, "1.4", "9,5,1", "Pul/SNR", "200030050"),
    (53, 10.0,  8.0, 0.33, 32.2, 100, 1.40, 1.80, 1.4, 2.30, 1.3, 3.33, "1.4", "0,5,2", "Pul/SNR", "200040050"),
    (56,  6.0, 40.0, 0.33,  0.0,   0, 1.89, 2.10, 1.4, 2.10, 2.0, 1.10, "1",   "0,5,2", "Pul/SNR", "050100020"),
]

nums = np.array([r[0] for r in ROWS])
X = np.array([r[1:1 + len(COLS)] for r in ROWS], float)
cat = [r[1 + len(COLS):] for r in ROWS]


def minmax(a):
    lo, hi = a.min(0), a.max(0)
    rng = np.where(hi > lo, hi - lo, 1.0)
    return (a - lo) / rng


def report():
    Z = minmax(X)
    D = np.sqrt(((Z[:, None, :] - Z[None, :, :]) ** 2).sum(-1))
    np.fill_diagonal(D, np.inf)

    print("=" * 74)
    print("1. NEAR-DUPLICATES  (distance to nearest other model, normalised)")
    print("=" * 74)
    nn = D.min(1)
    order = np.argsort(nn)
    for i in order:
        j = int(np.argmin(D[i]))
        bar = "#" * max(int(nn[i] * 40), 1)
        print(f"  #{nums[i]:<3} nearest is #{nums[j]:<3} at {nn[i]:.3f}  {bar}")
    print(f"\n  min {nn.min():.3f}   median {np.median(nn):.3f}   max {nn.max():.3f}")
    print(f"  ratio max/min = {nn.max()/nn.min():.1f}x  (1.0 would be perfectly even)")

    print()
    print("=" * 74)
    print("2. WHAT ACTUALLY DIFFERS between the closest pair")
    print("=" * 74)
    i, j = np.unravel_index(np.argmin(D), D.shape)
    for k, c in enumerate(COLS):
        a, b = X[i, k], X[j, k]
        mark = "   <-- differs" if not np.isclose(a, b) else ""
        print(f"  {c:<7} #{nums[i]}={a:<8g} #{nums[j]}={b:<8g}{mark}")
    same_cat = [c for c, a, b in zip(CATS, cat[i], cat[j]) if a == b]
    print(f"  categorical: identical in {len(same_cat)}/{len(CATS)} ({', '.join(same_cat)})")

    print()
    print("=" * 74)
    print("3. FAMILY CLUSTERING  (the categorical axes the distance can't see)")
    print("=" * 74)
    from collections import Counter
    key = [f"{c[2]:<8} ISRF={c[0]:<4} gas={c[1]:<6} B={c[3]}" for c in cat]
    for k, n in Counter(key).most_common():
        members = [str(nums[t]) for t in range(len(ROWS)) if key[t] == k]
        print(f"  {n:2d}/12  {k}   models {', '.join(members)}")
    print(f"\n  {len(set(key))} distinct configurations among 12 models")


def why_maxsum_fails():
    print()
    print("=" * 74)
    print("4. WHY 'LARGEST SUM OF DISTANCES' PICKS DUPLICATES")
    print("=" * 74)
    rng = np.random.default_rng(0)
    # a cloud plus two outliers sitting right next to each other
    pool = np.vstack([rng.normal(0.5, 0.08, (60, 2)),
                      np.array([[3.0, 3.0], [3.02, 2.98]])])
    seed = [0, 1, 2]

    def greedy(obj, k=4):
        sel = list(seed)
        for _ in range(k):
            best, bi = -np.inf, None
            for c in range(len(pool)):
                if c in sel:
                    continue
                d = np.linalg.norm(pool[sel] - pool[c], axis=1)
                v = d.sum() if obj == "sum" else d.min()
                if v > best:
                    best, bi = v, c
            sel.append(bi)
        return sel[len(seed):]

    for obj in ("sum", "min"):
        picked = greedy(obj)
        got_both = 60 in picked and 61 in picked
        name = "max-SUM (their method)" if obj == "sum" else "max-MIN (farthest-point)"
        print(f"  {name:26} picked {picked}"
              f"   {'<-- took BOTH twin outliers' if got_both else '<-- took at most one twin'}")
    print("\n  the two twins are 0.03 apart in a space spanning ~3.5.")
    print("  max-sum takes both because each is far from the CLOUD.")
    print("  max-min takes one, then looks elsewhere.")


if __name__ == "__main__":
    report()
    why_maxsum_fails()
