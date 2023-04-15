import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess
from tqdm import tqdm
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
import pandas
import time


sns.set_style("ticks")
sns.set_context("paper")

plt.rcParams["figure.figsize"] = [4, 1.2]
plt.rcParams["font.family"] = ["arial"]
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["xtick.labelsize"] = 7
plt.rcParams["ytick.labelsize"] = 7
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


def rev_compl_str(s):
    s_ = "".join([rev_dict[l] for l in s])[::-1]
    return s_


def rev_compl_mut(mut_tuple):
    ini, mut = mut_tuple
    ini_ = "".join([rev_dict[l] for l in ini])[::-1]
    mut_ = "".join([rev_dict[l] for l in mut])[::-1]
    return (ini_, mut_)


def check_argv():
    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument("--input", help="input folder with .vcf files")
    parser.add_argument("--feature", help=".bed feature file path")
    parser.add_argument("--out", help="output folder name")
    parser.add_argument("--hg", help="path to reference genome .fa")
    parser.add_argument("--hgver", help="hg version, either hg19 or hg38")
    parser.add_argument(
        "--cosmicver", help="COSMIC version, either 3.3, 3.0, 2.0 or PCAWG"
    )
    parser.add_argument(
        "--activities", help="file with COSMIC activities for given samples"
    )
    parser.add_argument("--perm", help="number of permutations")
    parser.add_argument("--chrom_lengths", help="chromosome length file")
    parser.add_argument("--temp_folder", help="temporary folder name")
    parser.add_argument(
        "--correct_trinucl",
        action="store_true",
        help="correct for trinucleotide context",
    )

    return parser.parse_args()


args = check_argv()

print(
    f"Deriving probabilities of SBS occurence at {args.feature.split('.')[0]} regions"
)

if not os.path.exists(args.out):
    os.makedirs(args.out)

if not os.path.exists(args.temp_folder):
    os.makedirs(args.temp_folder)

dict_hgver = {"hg19": "GRCh37", "hg38": "GRCh38"}

if args.cosmicver == "3.3":
    df = pd.read_csv(f"data/COSMIC_v3.3.1_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "3.0":
    df = pd.read_csv(f"data/COSMIC_v3_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "2.0":
    df = pd.read_csv(f"data/COSMIC_v2_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "PCAWG":
    df = pd.read_csv(f"data/pcawg_published_reference.txt", sep="\t")

activities = pd.read_csv(args.activities, sep="\t", header=0, index_col=0)

if "cancer" in list(activities.columns):
    activities.drop(["cancer"], axis=1, inplace=True)

activities = activities.divide(activities.sum(axis=1).values, axis="rows")

sbs_active = activities.mean(axis=0)[activities.mean(axis=0) != 0].index

tuples_mut = [
    (f"{b5}{x}{b3}", f"{b5}{y}{b3}")
    for b5 in "AGTC"
    for x, y in zip(["C", "C", "C", "T", "T", "T"], ["A", "G", "T", "A", "C", "G"])
    for b3 in "AGTC"
]

region = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
]

proc = subprocess.Popen(["which", "shuffleBed"], stdout=subprocess.PIPE)
out = proc.stdout.read().decode("utf-8")
bedtools_exec = "/".join(out.strip("\n").split("/")[:-1])
print("bedtools executable path to be used:", bedtools_exec)

# convert into dict and assign numeric mutation codes
dict_mut = {k: v for k, v in zip(tuples_mut, range(len(tuples_mut)))}
df_mut = pd.DataFrame({k: [dict_mut[k]] for k in dict_mut.keys()}).T
dict_mut_inv = {k: v for k, v in zip(range(len(tuples_mut)), tuples_mut)}

df["Type1"] = [
    (
        s.split("[")[0] + s.split("[")[1][0] + s.split("]")[-1],
        s.split("[")[0] + s.split("]")[0][-1] + s.split("]")[-1],
    )
    for s in df["Type"]
]

df.set_index("Type1", inplace=True)
df.drop("Type", axis=1, inplace=True)
S = df.sum(axis=1)

df_ = df.div(S, axis=0)

rev_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

# dictionary to store hg sequences
seq_dict = {}

print("Importing the genome...")
with open(args.hg, mode="r") as handle:
    # process each record in .fa file if there's more than one
    for record in SeqIO.parse(handle, "fasta"):
        identifier = record.id
        description = record.description
        # ignore alternative contigs
        if identifier in region:
            # print(description)
            sequence = record.seq
            seq_dict[identifier] = str(sequence)

# the input folder is hardcoded here
samples = os.listdir(args.input)
samples_ = []
for s in samples:
    if ".vcf" in s:
        samples_.append(s)

sample_list_clean = []
dict_sign_contrib = {}

enh = pd.read_csv(args.feature, sep="\t", header=None)
enh_ = enh

n_enh = len(enh_[1])
enh_["length"] = enh[2] - enh[1]

feature_coords = {i: [x, y] for i, x, y in zip(enh_[3], enh_[1], enh_[2])}

enh_["newstart"] = enh_[1]
enh_["newend"] = enh_[2]
enh_["id"] = [f"id{i}" for i in range(len(enh_))]
enh_pad = enh[[0, "newstart", "newend", "id"]]
enh_pad.to_csv(
    f"{args.temp_folder}/feature_padded.bed", sep="\t", index=None, header=None
)

if args.correct_trinucl:
    print("Calculating the odds for the tri-nucleotide frequency corection...")
    df_tri_freq = pd.read_csv("data/trinucleotide_freq.csv", index_col=0, header=0)

    feature_seqs = []
    for c, s, e in zip(enh_[0], enh[1], enh[2]):
        if c in seq_dict.keys():
            feature_seqs.append(seq_dict[c][s:e].upper())

    feature_freqs = []
    for tri in df_tri_freq["trinucl"]:
        n_occur = 0
        len_all = 0
        for fs in feature_seqs:
            n_occur += fs.count(tri)
            len_all += len(fs)
        feature_freqs.append(n_occur / len_all)

    df_tri_freq["freq_feature"] = feature_freqs
    df_tri_freq["OR_feature"] = [
        x / y if y != 0 else 1.0
        for x, y in zip(df_tri_freq["freq_feature"], df_tri_freq["frequency"])
    ]

    dict_tri_odds = dict(zip(df_tri_freq["trinucl"], df_tri_freq["OR_feature"]))

for sample in tqdm(samples_):
    dict_residuals = {}
    dict_residuals_cont = {}
    st = time.time()

    sample_list_clean.append(sample)

    mut1 = pd.read_csv(
        os.path.join(args.input, sample),
        sep="\t",
        comment="#",
        header=None,
        low_memory=False,
    )
    mut1["chr"] = ["chr" + str(x) for x in mut1[0]]
    mut1["start"] = [x - 1 for x in mut1[1]]
    mut1["end"] = mut1[1]

    mut_numbers = []
    mut_types = []
    mut_contexts = []
    strands = []

    for c, s, r in zip(mut1["chr"], mut1[1], mut1[3]):
        if seq_dict[c][s - 1].upper() == r:
            strands.append("+")
        elif rev_dict[seq_dict[c][s - 1].upper()] == r:
            strands.append("-")
        else:
            strands.append(np.nan)

    mut1["strand"] = strands

    for c, ss, st, r, a in zip(mut1["chr"], mut1["strand"], mut1[1], mut1[3], mut1[4]):
        context = seq_dict[c][st - 2 : st + 1].upper()
        ty = (context, context[0] + a + context[-1])
        if ss == "-":
            context = rev_compl_str(context)
            ty = rev_compl_mut(ty)
        mut_contexts.append(context)
        mut_types.append(ty)

    mut1["context"] = mut_contexts
    mut1["type"] = mut_types
    mut_numbers = []
    strands = []

    for name in mut1["type"]:
        if name in dict_mut.keys():
            mut_numbers.append(dict_mut[name])
        elif rev_compl_mut(name) in dict_mut.keys():
            mut_numbers.append(dict_mut[rev_compl_mut(name)])
        else:
            mut_numbers.append(None)

    mut1["number"] = mut_numbers

    mut2 = mut1[["chr", "start", "end", "number", "strand"]].dropna()
    mut2.to_csv(
        os.path.join(args.temp_folder, sample.split(".")[0] + ".bed"),
        sep="\t",
        index=None,
        header=None,
    )

    all_mut = pd.read_csv(
        os.path.join(args.temp_folder, sample.split(".")[0] + ".bed"),
        sep="\t",
        header=None,
    )

    f = open(f"{args.temp_folder}/feature_vcf_intersect.bed", "w")
    subprocess.call(
        [
            f"{bedtools_exec}/intersectBed",
            "-b",
            f"{args.temp_folder}/{sample.split('.')[0]}.bed",
            "-a",
            f"{args.temp_folder}/feature_padded.bed",
            "-wa",
            "-wb",
        ],
        stdout=f,
    )
    f = open(f"{args.temp_folder}/feature_vcf_nointersect.bed", "w")
    subprocess.call(
        [
            f"{bedtools_exec}/intersectBed",
            "-b",
            f"{args.temp_folder}/{sample.split('.')[0]}.bed",
            "-a",
            f"{args.temp_folder}/feature_padded.bed",
            "-v",
        ],
        stdout=f,
    )

    for k in range(int(args.perm)):
        print(f"Running random permutation {k} for context generation...")
        f = open(f"{args.temp_folder}/feature_context_{k}.bed", "w")
        subprocess.call(
            [
                f"{bedtools_exec}/shuffleBed",
                "-i",
                f"{args.temp_folder}/feature_padded.bed",
                "-g",
                f"data/{args.hgver}.genome",
            ],
            stdout=f,
        )
        f = open(f"{args.temp_folder}/context_vcf_intersect_{k}.bed", "w")
        subprocess.call(
            [
                f"{bedtools_exec}/intersectBed",
                "-b",
                f"{args.temp_folder}/{sample.split('.')[0]}.bed",
                "-a",
                f"{args.temp_folder}/feature_context_{k}.bed",
                "-wa",
                "-wb",
            ],
            stdout=f,
        )
        f = open(f"{args.temp_folder}/context_vcf_nointersect_{k}.bed", "w")
        subprocess.call(
            [
                f"{bedtools_exec}/intersectBed",
                "-b",
                f"{args.temp_folder}/{sample.split('.')[0]}.bed",
                "-a",
                f"{args.temp_folder}/feature_context_{k}.bed",
                "-v",
            ],
            stdout=f,
        )

    try:
        enh_int = pd.read_csv(
            f"{args.temp_folder}/feature_vcf_intersect.bed", sep="\t", header=None
        )
        enh_int_none = pd.read_csv(
            f"{args.temp_folder}/feature_vcf_nointersect.bed", sep="\t", header=None
        )
        for k in range(int(args.perm)):
            enh_int_c = pd.read_csv(
                f"{args.temp_folder}/context_vcf_intersect_{k}.bed",
                sep="\t",
                header=None,
            )
            enh_int_cn = pd.read_csv(
                f"{args.temp_folder}/context_vcf_nointersect_{k}.bed",
                sep="\t",
                header=None,
            )
            enh_int_c["id"] = [x + f"_{k}" for x in enh_int_c[3]]
            if k == 0:
                enh_int_cont = enh_int_c
                enh_int_cont_none = enh_int_cn
            else:
                enh_int_cont = pd.concat([enh_int_cont, enh_int_c])
                enh_int_cont_none = pd.concat([enh_int_cont_none, enh_int_cn])
    except pandas.errors.EmptyDataError:
        continue

    enh_int["mut_name"] = [dict_mut_inv[mut_number] for mut_number in enh_int[7]]
    enh_int_cont["mut_name"] = [
        dict_mut_inv[mut_number] for mut_number in enh_int_cont[7]
    ]
    enh_int["length"] = enh_int[2] - enh_int[1]
    enh_int_cont["length"] = enh_int_cont[2] - enh_int_cont[1]

    for col in sbs_active:
        # print(f"Processing {col}...")
        enh_int[f"weight_{col}"] = [
            float(df_.loc[[mut_name]][col]) for mut_name in enh_int["mut_name"]
        ]
        enh_int_cont[f"weight_{col}"] = [
            float(df_.loc[[mut_name]][col]) for mut_name in enh_int_cont["mut_name"]
        ]

        if args.correct_trinucl:
            enh_int[f"weight_{col}"] = [
                w / dict_tri_odds[mut_name[0]]
                for mut_name, w in zip(enh_int["mut_name"], enh_int[f"weight_{col}"])
            ]
            enh_int_cont[f"weight_{col}"] = [
                w / dict_tri_odds[mut_name[0]]
                for mut_name, w in zip(
                    enh_int_cont["mut_name"], enh_int_cont[f"weight_{col}"]
                )
            ]

        residuals = []
        residuals_cont = []
        activity = activities.loc[sample.split(".")[0]][col]
        # tot_length = enh_int_cont["length"].sum()
        # R_fid_cont = enh_int_cont[f"weight_{col}"].sum() * activity / tot_length
        # collated_feature_bed = {"chr": [], "start": [], "end": [], "activity": []}
        for fid in enh_["id"].unique():
            enh_int_fid = enh_int[enh_int[3] == fid]
            if len(enh_int_fid) > 0:
                # collated_feature_bed["chr"].append(enh_int_fid.iloc[0][0])
                # collated_feature_bed["start"].append(enh_int_fid.iloc[0][1])
                # collated_feature_bed["end"].append(enh_int_fid.iloc[0][2])
                length = int(enh_int_fid.iloc[0]["length"])
                R_fid = enh_int_fid[f"weight_{col}"].sum() * activity / length
                # collated_feature_bed["activity"].append(R_fid)
            else:
                length = 1
                R_fid = 0

            residuals.append(R_fid)
        residuals.extend([0] * len(enh_int_none))

        for fid in enh_int_cont["id"].unique():
            enh_int_cont_fid = enh_int_cont[enh_int_cont["id"] == fid]
            if len(enh_int_cont_fid) > 0:
                # collated_feature_bed["chr"].append(enh_int_fid.iloc[0][0])
                # collated_feature_bed["start"].append(enh_int_fid.iloc[0][1])
                # collated_feature_bed["end"].append(enh_int_fid.iloc[0][2])
                length = int(enh_int_cont_fid.iloc[0]["length"])
                R_fid_cont = enh_int_cont_fid[f"weight_{col}"].sum() * activity / length
                # collated_feature_bed["activity"].append(R_fid)
            else:
                length = 1
                R_fid_cont = 0

            residuals_cont.append(R_fid_cont)
        residuals_cont.extend([0] * len(enh_int_cont_none))

        dict_residuals[f"{col}"] = residuals
        dict_residuals_cont[f"{col}"] = residuals_cont
        # collated_feature_bed = pd.DataFrame(collated_feature_bed)
        # collated_feature_bed.to_csv(
        #    os.path.join("temp", f"{sample.split('.')[0]}_{col}.bed"),
        #    sep="\t",
        #    index=None,
        # )

    df_residuals = pd.DataFrame(dict_residuals)
    # df_residuals.to_csv(os.path.join(args.out, f"{sample.split('.')[0]}_probas.csv"))
    df_residuals_sp = scipy.sparse.csr_matrix(df_residuals.values)
    scipy.sparse.save_npz(
        os.path.join(args.out, f"{sample.split('.')[0]}_probas.csv"),
        df_residuals_sp,
    )

    df_residuals_cont = pd.DataFrame(dict_residuals_cont)
    df_residuals_cont_sp = scipy.sparse.csr_matrix(df_residuals_cont.values)
    scipy.sparse.save_npz(
        os.path.join(args.out, f"{sample.split('.')[0]}_probas_cont.csv"),
        df_residuals_cont_sp,
    )

et = time.time()
elapsed_time = et - st
print("Execution time:", elapsed_time, "seconds")
