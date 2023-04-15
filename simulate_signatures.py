from Bio import SeqIO
from tqdm import tqdm

import pysam
import os
import argparse
import sys
import re
import pandas as pd
import numpy as np
import random


# check arguments
def check_argv():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--signatures", nargs="+", type=str, help="signatures to simulate"
    )
    parser.add_argument(
        "--activities",
        nargs="+",
        type=float,
        help="activities to simulate; must have the same length as signatures",
    )
    parser.add_argument(
        "--FCs",
        nargs="+",
        type=float,
        help="FCs to simulate; must have the same length as signatures",
    )
    parser.add_argument(
        "--feature",
        required=False,
        type=str,
        help="genomic feature to populate",
    )
    parser.add_argument(
        "--N_samples",
        default=1000,
        required=True,
        type=int,
        help="number of sample to simulate",
    )
    parser.add_argument(
        "--N_mut",
        default=5000,
        required=True,
        type=int,
        help="approximate number of mutations to populate per sample",
    )
    parser.add_argument(
        "--N_noise",
        default=1000,
        required=True,
        type=int,
        help="number of noise mutations",
    )
    parser.add_argument(
        "--hg19", default="hg19.fa", required=True, type=str, help="path to hg19.fa"
    )
    parser.add_argument("--hgver", help="hg version, either hg19 or hg38")
    parser.add_argument(
        "--cosmicver",
        required=True,
        help="COSMIC version, either 3.3, 3.0, 2.0 or PCAWG",
    )
    parser.add_argument(
        "--output",
        default="out",
        required=True,
        type=str,
        help="output folder name",
    )

    return parser.parse_args()


def rev_compl_mut(mut_tuple):
    rev_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    ini, mut = mut_tuple
    ini_ = "".join([rev_dict[l] for l in ini])[::-1]
    mut_ = "".join([rev_dict[l] for l in mut])[::-1]

    return (ini_, mut_)


def rev_compl(s):
    rev_dict = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    s_ = "".join([rev_dict[l] for l in ini])[::-1]
    return s_


dict_hgver = {"hg19": "GRCh37", "hg38": "GRCh38"}


args = check_argv()

assert len(args.signatures) == len(
    args.activities
), f"expected the same number of signatures and activities, got {len(args.signatures)} and {len(args.activities)}"

assert len(args.signatures) == len(
    args.FCs
), f"expected the same number of signatures and activities, got {len(args.signatures)} and {len(args.FCs)}"

activities = [a / np.sum(np.array(args.activities)) for a in args.activities]

if args.cosmicver == "3.3":
    df = pd.read_csv(f"data/COSMIC_v3.3.1_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "3.0":
    df = pd.read_csv(f"data/COSMIC_v3_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "2.0":
    df = pd.read_csv(f"data/COSMIC_v2_SBS_{dict_hgver[args.hgver]}.txt", sep="\t")
elif args.cosmicver == "PCAWG":
    df = pd.read_csv(f"data/pcawg_published_reference.txt", sep="\t")

df["Type1"] = [
    (
        s.split("[")[0] + s.split("[")[1][0] + s.split("]")[-1],
        s.split("[")[0] + s.split("]")[0][-1] + s.split("]")[-1],
    )
    for s in df["Type"]
]

df.set_index("Type1", inplace=True)
df.drop("Type", axis=1, inplace=True)

# uncomment if this is not a test run and you want to process all chromosomes
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

# dictionary to store hg19 sequences
seq_dict = {}

with open(args.hg19, mode="r") as handle:
    # process each record in .fa file if there's more than one
    for record in SeqIO.parse(handle, "fasta"):
        identifier = record.id
        description = record.description
        # ignore alternative contigs
        if identifier in region:
            print(description)
            sequence = record.seq
            seq_dict[identifier] = str(sequence)

# create output directory if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

# list of mutation types in tuples; ('CTA', 'CAA') means from 'CTA' > to 'CAA'
tuples_mut = [
    (f"{b5}{x}{b3}", f"{b5}{y}{b3}")
    for b5 in "AGTC"
    for x, y in zip(["C", "C", "C", "T", "T", "T"], ["A", "G", "T", "A", "C", "G"])
    for b3 in "AGTC"
]

channels = [x[0] for x in tuples_mut]

# convert into dict and assign numeric mutation codes
dict_mut = {k: v for k, v in zip(tuples_mut, range(len(tuples_mut)))}
df_mut = pd.DataFrame({k: [dict_mut[k]] for k in dict_mut.keys()})
# save the numeric codes of single base mismatches
df_mut.to_csv(os.path.join(args.output, "sbs_codes.csv"))

extended_dict = {}
total_length = 0

# loop over chromosomes: whatever we saved into seq_dict according to chromosomal region of interest
for chr_id in tqdm(seq_dict.keys()):
    # remove "chr" from chromosome name if needed
    print(f"......Processing chromosome {chr_id}.......")
    ref = seq_dict[chr_id].upper()
    extended_dict[chr_id] = {}
    total_length += len(ref)

    for x in tqdm(range(1, len(ref) - 2)):
        extended_dict[chr_id][x] = ref[x - 1 : x + 2]

inverse_dict = {x: [] for x in channels}

for chr_id in tqdm(extended_dict.keys()):
    print(f"......Inverting the genomic dictionary {chr_id}.......")
    for loc in extended_dict[chr_id].keys():
        if extended_dict[chr_id][loc] in inverse_dict.keys():
            inverse_dict[extended_dict[chr_id][loc]].append((chr_id, loc))
        if rev_compl(extended_dict[chr_id][loc]) in inverse_dict.keys():
            inverse_dict[extended_dict[chr_id][loc]].append((chr_id, loc))

# construct a dict with feature
feature = pd.read_csv(args.feature, sep="\t", header=None)

feature_dict = {}
feature_length = 0
print(f"......Processing genomic feature.......")
feature_dict = {chr_id: {} for chr_id in feature[0].unique()}
# loop over chromosomes: whatever we saved into seq_dict according to chromosomal region of interest
for chr_id, start, end in tqdm(zip(feature[0], feature[1], feature[2])):
    ref = seq_dict[chr_id][start:end].upper()
    feature_length += len(ref)
    for x in tqdm(range(1, len(ref) - 2)):
        feature_dict[chr_id][start + x] = ref[x - 1 : x + 2]


inverse_dict_feature = {x: [] for x in channels}

for chr_id in tqdm(feature_dict.keys()):
    print(f"......Inverting the feature dictionary {chr_id}.......")
    for loc in feature_dict[chr_id].keys():
        if feature_dict[chr_id][loc] in inverse_dict_feature.keys():
            inverse_dict_feature[feature_dict[chr_id][loc]].append((chr_id, loc))
        if rev_compl(feature_dict[chr_id][loc]) in inverse_dict_feature.keys():
            inverse_dict_feature[feature_dict[chr_id][loc]].append((chr_id, loc))


N_bck = {x: [] for x in args.signatures}
N_feat = {x: [] for x in args.signatures}

for i, sign in enumerate(args.signatures):
    N_mut_signature = args.N_mut * activity[i]
    N_bck[sign] = int(
        N_mut_signature / (1 + args.FCs[i] * feature_length / total_length)
    )
    N_feat[sign] = N_mut_signature - N_bck[sign]

print("Background mutations:", N_bck, ", feature mutations:", N_feat)

for k in args.N_samples:
    # now we need to count mutation types by channel from the COSMIC table and sample the respective trinucl,
    # populate a new dictionary with the trinucl location and (from)>(to) channels for the final vcf
    sampled_bck = {x: 0 for x in df.index}
    sampled_feat = {x: 0 for x in df.index}

    for i, sign in enumerate(args.signatures):
        signature_dict = df[sign]
        for ch, v in zip(signature_dict.index, signature_dict.values):
            N_b = int(N_bck[sign] * v)
            if N_b > 0:
                sampled_bck[ch] += 1
            N_f = int(N_feat[sign] * v)
            if N_f > 0:
                sampled_feat[ch] += 1

    # sampling
    loci, chs = [], [], []
    for ch in sampled_bck.keys():
        chs.extend([ch] * sampled_bck[ch])
        loci.extend(random.sample(inverse_dict[ch[0]], sampled_bck[ch]))

    for ch in sampled_feat.keys():
        chs.extend([ch] * sampled_feat[ch])
        loci.extend(random.sample(inverse_dict_feature[ch[0]], sampled_feat[ch]))

    # VCF: CHROM POS ID REF ALT QUAL FILTER INFO
    # 1 2261023 .   C   G   .   .   Callers=dkfz
    df_vcf = pd.DataFrame(
        {
            "chr": [x[0][3:] for x in loci],
            "pos": [x[1] + 1 for x in loci],
            "ID": ["."] * len(loci),
            "ref": [x[0][1] for x in chs],
            "alt": [x[1][1] for x in chs],
            "qual": ["."] * len(loci),
            "filter": ["."] * len(loci),
            "info": ["None"] * len(loci),
        }
    )
    df_vcf.to_csv(
        os.path.join(args.output, f"simulated_{k}.vcf"),
        sep="\t",
        header=None,
        index=None,
    )
