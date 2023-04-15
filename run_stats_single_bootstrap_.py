from tqdm import tqdm
import os
import argparse
import pandas as pd
import numpy as np
import scipy.sparse
from scipy.stats import fisher_exact, mannwhitneyu, ttest_1samp
import seaborn as sns
import matplotlib.pyplot as plt
import random
import seaborn as sns
import matplotlib.pyplot as plt
import math


sns.set_style("ticks")
sns.set_context("paper")

plt.rcParams["figure.figsize"] = [1.1, 5]  # [1.1, 3.5]
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
    parser.add_argument("--folder", help="folder with data")
    parser.add_argument("--activities", help="txt file with COSMIC activities")
    parser.add_argument("--p_cutoff", default=0.05, type=float, help="p-value cutoff")
    parser.add_argument("--N_perm", default=1000, type=int, help="number of bootstrapping experiments")
    parser.add_argument("--N_sampled", default=10000, type=int, help="number of samples in one bootstrap")

    return parser.parse_args()


args = check_argv()

print(f"Running tests for individual samples...")

activities = pd.read_csv(
    args.activities, sep="\t", header=0, index_col=0
)

newindex = activities.mean(axis=0)[activities.mean(axis=0) != 0].index

print(activities.mean(axis=0)[activities.mean(axis=0) != 0])

folder = args.folder
dirs = os.listdir(folder)
samples_ = []

for f in dirs:
    if "probas_cont.csv.npz" in f:
        samples_.append(f)

# how many bootstrapping operations we perform and how many bootstrap samples we will draw at once
N_perm = args.N_perm
N_sampled = args.N_sampled

for m in tqdm(range(len(samples_))):

    f = samples_[m]

    df_2 = {c: [] for c in newindex}
    df_context = {}

    sample = f.split(".")[0]
    s = sample.split('_')[0]
    npz = scipy.sparse.load_npz(os.path.join(folder, f)).toarray()

    if N_sampled>np.shape(npz)[0]:
        N_sampled = int(np.shape(npz)[0]/2)

    found_s1 = False
    for f2 in dirs:
        sample1 = f2.split(".")[0]
        s1 = sample1.split('_')[0]
        if "probas.csv.npz" in f2 and s1 == s:
            npz2 = scipy.sparse.load_npz(os.path.join(folder, f2)).toarray()
            found_s1 = True

    if not found_s1:
        print(f"Have not found the context file for {s}", sample)
        continue

    for j, c in enumerate(newindex):
        df_context[c] = np.mean(npz[:, j])

    for i in range(N_perm):
        rand_ind = random.sample(range(np.shape(npz)[0]), N_sampled)
        feat = npz2[rand_ind, :]

        for j, c in enumerate(newindex):
            fc = np.mean(feat[:, j]) / df_context[c]
            if not np.isnan(fc):
                df_2[c].append(fc)

    fc_bootstrap = []
    pp_bootstrap = []

    for j, c in enumerate(newindex):
        fc = np.mean(df_2[c])
        res = ttest_1samp(df_2[c], df_context[c])
        pp_bootstrap.append(res[1])
        fc_bootstrap.append(fc)

    df_persamp_ = pd.DataFrame({'FC': fc_bootstrap,
                                'p-value': pp_bootstrap}, index=newindex)

    df_persamp_.to_csv(os.path.join(folder, f"{s}_SBS_bootstrap.csv"))
    df_persamp_.head()
