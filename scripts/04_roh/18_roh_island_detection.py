import sys
import pandas as pd

pop = sys.argv[1]
npop = int(sys.argv[2])
thr = 0.5 * npop
gap_bp = int(sys.argv[3])  # jarak maksimum untuk menggabungkan SNP jadi 1 region

df = pd.read_csv(f"roh_occ_{pop}.tsv", sep="\t")
df = df[df["ROH_OCC"] >= thr].sort_values(["CHR","POS"]).reset_index(drop=True)

out = []
for chr_id, sub in df.groupby("CHR"):
    positions = sub["POS"].tolist()
    if not positions:
        continue
    start = prev = positions[0]
    for p in positions[1:]:
        if p - prev <= gap_bp:
            prev = p
        else:
            out.append((chr_id, start, prev))
            start = prev = p
    out.append((chr_id, start, prev))

bed = pd.DataFrame(out, columns=["CHR","START","END"])
bed.to_csv(f"ROH_islands_{pop}.bed", sep="\t", index=False, header=True)
print(f"{pop}: threshold={thr} (N={npop}), islands={len(bed)} saved to ROH_islands_{pop}.bed")
