import sys, numpy as np, pandas as pd

snps_file, seg_file, out_file = sys.argv[1], sys.argv[2], sys.argv[3]

snps = pd.read_csv(snps_file, sep=r"\s+", header=None, names=["CHR","POS","SNPID"])
segs = pd.read_csv(seg_file,  sep=r"\s+", header=None, names=["CHR","START","END"])

snps["CHR"] = snps["CHR"].astype(int)
segs["CHR"] = segs["CHR"].astype(int)

snps = snps.sort_values(["CHR","POS"]).reset_index(drop=True)
segs = segs.sort_values(["CHR","START","END"]).reset_index(drop=True)

occ = np.zeros(len(snps), dtype=np.int32)

for chr_id, idx in snps.groupby("CHR").groups.items():
    pos = snps.loc[idx, "POS"].to_numpy()
    chr_segs = segs[segs["CHR"]==chr_id]
    if chr_segs.empty:
        continue
    starts = chr_segs["START"].to_numpy()
    ends   = chr_segs["END"].to_numpy()
    starts.sort()
    ends.sort()

    i = j = 0
    active = 0
    for k, p in zip(idx, pos):
        while i < len(starts) and starts[i] <= p:
            active += 1; i += 1
        while j < len(ends) and ends[j] < p:
            active -= 1; j += 1
        occ[k] = active

snps["ROH_OCC"] = occ
snps.to_csv(out_file, sep="\t", index=False)
