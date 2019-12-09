import os
import sys
import pandas as pd

cosmic2_dir = "cosmic_v2"
cosmic3_dir = "cosmic_v3"

# -------------------------------
# COSMIC2
# -------------------------------
cosmic2 = pd.read_csv(os.path.join(cosmic2_dir, "cosmic_v2.txt"),sep="\t").dropna(1)
cosmic2.to_csv(os.path.join(cosmic2_dir, "sa_cosmic2.tsv"),sep="\t")

# -------------------------------
# COSMIC3 SBS Exome
# -------------------------------
cosmic3_dir = "/home/sanand/getzlab-SignatureAnalyzer/siganalyzer/ref/cosmic_v3/"

sbs_ex = pd.read_csv(os.path.join(cosmic3_dir, "sigProfiler_exome_SBS_signatures.csv"))
sbs_ex = sbs_ex.rename(columns={'Type':'Substitution Type', 'SubType':'Trinucleotide'})

def func(row):
    tnc = row['Trinucleotide']
    return "{}[{}]{}".format(tnc[0],row['Substitution Type'],tnc[-1])

sbs_ex['Somatic Mutation Type'] = sbs_ex.apply(func,1)

_reorder = ['Substitution Type','Trinucleotide','Somatic Mutation Type']
_reorder += [x for x in list(sbs_ex) if x.startswith("SBS")]

sbs_ex = sbs_ex.loc[:,_reorder]
sbs_ex.to_csv(os.path.join(cosmic3_dir, "sa_cosmic3_sbs_exome.tsv"),sep="\t")

# -------------------------------
# COSMIC3 SBS
# -------------------------------
sbs = pd.read_csv(os.path.join(cosmic3_dir, "sigProfiler_SBS_signatures_2019_05_22.csv"))
sbs = sbs.rename(columns={'Type':'Substitution Type', 'SubType':'Trinucleotide'})

sbs['Somatic Mutation Type'] = sbs.apply(func,1)

_reorder = ['Substitution Type','Trinucleotide','Somatic Mutation Type']
_reorder += [x for x in list(sbs) if x.startswith("SBS")]

sbs = sbs.loc[:,_reorder]
sbs.to_csv(os.path.join(cosmic3_dir, "sa_cosmic3_sbs.tsv"),sep="\t")

# -------------------------------
# COSMIC3 SBS-T
# -------------------------------
sbs_t = pd.read_csv(os.path.join(cosmic3_dir, "sigProfiler_TSB_signatures.csv"))
sbs_t = sbs_t.rename(columns={'Type':'Substitution Type', 'Subtype':'Trinucleotide'})

sbs_t['Somatic Mutation Type'] = sbs_t.apply(func,1)

_reorder = ['Strand','Substitution Type','Trinucleotide','Somatic Mutation Type']
_reorder += [x for x in list(sbs_t) if x.startswith("SBS")]

sbs_t = sbs_t.loc[:,_reorder]

sbs_t.head()
sbs_t.to_csv(os.path.join(cosmic3_dir, "sa_cosmic3_sbs_t.tsv"),sep="\t")

# -------------------------------
# COSMIC 3 DBS
# -------------------------------
dbs = pd.read_csv(os.path.join(cosmic3_dir, "sigProfiler_DBS_signatures.csv"))
dbs = dbs.rename(columns={'Mutation Type':'Somatic Mutation Type'})

dbs.to_csv(os.path.join(cosmic3_dir, "sa_cosmic3_dbs.tsv"),sep="\t")
