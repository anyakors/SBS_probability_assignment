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

plt.rcParams["figure.figsize"] = [4, 5]
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

activities = pd.read_csv(args.activities, sep="\t", header=0, index_col=0)

if "cancer" in list(activities.columns):
    activities.drop(["cancer"], axis=1, inplace=True)

activities = activities.divide(activities.sum(axis=1).values, axis="rows")
newindex = activities.mean(axis=0)[activities.mean(axis=0) != 0].index

print(activities.mean(axis=0)[activities.mean(axis=0) != 0])

# derive FC from intrasample values

folder = args.folder
dirs = os.listdir(folder)

fc_b = []
act_b = []
sbs_b = []

for f in dirs:
    if "bootstrap.csv" in f:
        sample_name = f[:-18]
        df_temp = pd.read_csv(os.path.join(folder, f), header=0, index_col=0)
        df_temp = df_temp.dropna(axis=0)
        print(sample_name)
        if sample_name in list(activities.index):
            for i in df_temp.index:
                fc_b.append(float(df_temp.loc[i]["FC"]))
                sbs_b.append(i)
                act_b.append(
                    activities.loc[sample_name][i] / activities.loc[sample_name].sum()
                )

df_b = pd.DataFrame(
    {
        "SBS": sbs_b,
        "FC": fc_b,
        "Activities": act_b,
        "FC_weighted": [x * y for x, y in zip(fc_b, act_b)],
    }
)
# df_b.to_csv(os.path.join(folder, f"all_SBS_allsamples_FC_b.csv"))

print(df_b)

fcs = []
for c in newindex:
    fcs.append(
        np.sum(df_b[df_b["SBS"] == c]["FC_weighted"])
        / np.sum(df_b[df_b["SBS"] == c]["Activities"])
    )
df_plot = pd.DataFrame({"signature": newindex, "FC": fcs})
df_plot.set_index("signature", inplace=True)

ss = []
pp = []

for c in newindex:
    # test = np.log2(df_2[df_2[c] != np.inf][c].dropna().values)
    test = np.log2(df_b[df_b["SBS"] == c]["FC"].dropna().values)
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

p_cutoff = args.p_cutoff

for i in range(len(df_plot.index)):
    if (
        df_p.iloc[i]["p"] == None
        or np.isnan(df_p.iloc[i]["p"])
        or df_p.iloc[i]["p"] > p_cutoff
    ):
        df_plot.iloc[i]["FC"] = np.nan

df_p["FC"] = df_plot["FC"]
df_p.set_index("SBS").to_csv(
    os.path.join(folder, f"all_SBS_intersample_bootstrap_weighted.csv")
)

plt.cla()
plt.clf()
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={"width_ratios": [3, 1]})

df_b["FC"] = [np.log2(x) if not np.isnan(x) else np.nan for x in df_b["FC"]]

ax1 = sns.boxplot(
    data=df_b,
    order=newindex,
    y="SBS",
    x="FC",
    whis=[5, 95],
    # orient="h",
    fliersize=1,
    notch=True,
    showcaps=False,
    flierprops={"marker": "x", "markeredgecolor": "salmon"},
    boxprops={"facecolor": (0.4, 0.6, 0.8, 0.5)},
    medianprops={"color": "gray"},
    ax=ax1,
)

ax1.set_xlabel(r"FC")

ax2 = sns.heatmap(
    df_plot,
    annot=True,
    fmt=".3f",
    cmap=sns.color_palette("ch:s=-.2,r=.6", as_cmap=True),
    vmin=0.1,
    vmax=2,
    ax=ax2,
)

plt.ylabel("")
ax2.set_yticklabels("")
ax2.set_yticks(ax1.get_yticks() + 0.5)

ax1.axvline(
    x=0,
    ymin=0,
    ymax=len(df_b["SBS"].unique()),
    lw=1.5,
    ls="--",
    color="red",
    zorder=-10,
)

print(ax1.get_xlim())
if ax1.get_xlim()[1] > 10 or ax1.get_xlim()[0] < -10:
    ax1.set_xlim(-10, 10)

xticks = ax1.get_xticks()
ax1.set_xticklabels([r"2$^{{{x}}}$".format(x=int(x)) for x in xticks])
ax2.set_xticklabels("")
ax2.set_xticks([])

plt.tight_layout(w_pad=-0.2)

plt.savefig(
    os.path.join(folder, f"FC_heatmap_fromBootstrap_weighted.pdf"),
    format="pdf",
    bbox_inches="tight",
)
