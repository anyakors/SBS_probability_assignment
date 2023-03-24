from tqdm import tqdm
import os
import argparse
import pandas as pd
import numpy as np
import scipy.sparse
from scipy.stats import fisher_exact, mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
import math


sns.set_style("ticks")
sns.set_context("paper")

plt.rcParams["figure.figsize"] = [1.1, 5] 
plt.rcParams["font.family"] = ["arial"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 8
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["ytick.labelsize"] = 8
plt.rcParams["axes.linewidth"] = 0.5
plt.rcParams["xtick.minor.width"] = 0.5
plt.rcParams["ytick.minor.width"] = 0.5
plt.rcParams["xtick.major.width"] = 0.5
plt.rcParams["ytick.major.width"] = 0.5
plt.rcParams["xtick.major.size"] = 3.5
plt.rcParams["ytick.major.size"] = 3.5
COLOR = "black"
plt.rcParams["text.color"] = COLOR
plt.rcParams["axes.labelcolor"] = COLOR
plt.rcParams["xtick.color"] = COLOR
plt.rcParams["ytick.color"] = COLOR


def check_argv():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--folder", help="folder with generated probability files")
    parser.add_argument("--activities", help="txt file with COSMIC activities")
    parser.add_argument("--p_cutoff", default=0.05, type=float, help="p-value cutoff")

    return parser.parse_args()


args = check_argv()

print(f"Running tests for mutational occupancy...")

activities = pd.read_csv(
    args.activities, sep="\t", header=0, index_col=0
)

newindex = activities.mean(axis=0)[activities.mean(axis=0) != 0].index

print(activities.mean(axis=0)[activities.mean(axis=0) != 0])

# derive FC from intrasample values

folder = args.folder
dirs = os.listdir(folder)

df_2 = {}

for f in dirs:
    if "probas_cont.csv.npz" in f:
        sample = f.split("_")[0]
        sample_no = f.split("_")[1]
        s = f"{sample}_{sample_no}"
        # cancer = f.split('_')[1]
        npz = scipy.sparse.load_npz(os.path.join(folder, f))
        mean_context = np.sum(npz.toarray(), axis=0)

        for f2 in dirs:
            sample1 = f2.split("_")[0]
            sample1_no = f2.split("_")[1]
            s1 = f"{sample1}_{sample1_no}"
            if "probas.csv.npz" in f2 and s1 == s:
                npz2 = scipy.sparse.load_npz(os.path.join(folder, f2))
                mean_feat = np.sum(npz2.toarray(), axis=0)

        if len(mean_feat) == len(newindex) and len(mean_context) == len(newindex):
            df_2[s] = mean_feat / mean_context

df_2 = pd.DataFrame(df_2).T
df_2.columns = newindex

fcs = []
for c in df_2.columns:
    fcs.append(df_2[(df_2[c] != np.inf) & (df_2[c] != np.nan)][c].mean())
df_plot = pd.DataFrame({"signature": df_2.columns, "FC": fcs})
df_plot.set_index("signature", inplace=True)

ss = []
pp = []

for c in df_2.columns:
    test = np.log2(df_2[df_2[c] != np.inf][c].dropna().values)
    control = np.array([np.log2(1)] * len(test))
    # res = stats.ttest_rel(control, test)
    if len(test) > 1:
        res = mannwhitneyu(test, control)
    else:
        res = (None, None)
    print(f"{c}: Mann-Whitney U test statistic: {res[0]}, p: {res[1]}")
    ss.append(res[0])
    pp.append(res[1])

df_p = pd.DataFrame({"MW-U": ss, "p": pp})
df_p["SBS"] = newindex
df_p.set_index("SBS").to_csv(os.path.join(folder, f"all_SBS_intersample.csv"))

p_cutoff = args.p_cutoff

for i in range(len(df_plot.index)):
    if (
        df_p.iloc[i]["p"] == None
        or np.isnan(df_p.iloc[i]["p"])
        or df_p.iloc[i]["p"] > p_cutoff
    ):
        df_plot.iloc[i]["FC"] = np.nan

sns.heatmap(
    df_plot,
    annot=True,
    fmt=".3f",
    cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True),
    vmin=0.1,
    vmax=2,
)
plt.yticks(rotation=0)

plt.savefig(
    os.path.join(folder, f"FC_heatmap.pdf"), format="pdf", bbox_inches="tight"
)
