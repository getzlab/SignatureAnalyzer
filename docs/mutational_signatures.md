## Mutational Signatures

Decomposition of mutational signatures using  `signatureanalyzer`. For a comprehensive  description of mutational signatures, their relevance, and references, please see the Catalogue of Somatic Mutations in Cancer, or COSMIC, [here](https://cancer.sanger.ac.uk/cosmic/signatures). The following document is a reference for important considerations when running this method.

---

#### Objective Function
For mutational signatures, we assume a poisson distribution of counts and use Fevotte & Tan's derivation of a poisson objective function for ARD-NMF. Thus, it is important to use the default value for the objective function (`poisson`).

Use:
```{bash}
--objective poisson
```
---

#### Human Genome (Hg) Build
Select which human genome build to use for mapping. We build base contexts using a 2-bit representation of the genome build. These may be downloaded here:
* hg38: `wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit`
* hg19: `wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit`

Use:
```{bash}
--hg_build <PATH>/hg19.2bit
```

---

#### COSMIC Signatures
Signature Analyzer supports encoding of:
* Single Base Substitution (SBS) Signatures (WGS: `cosmic3`, WES: `cosmic3_exome`)
* Doublet Base Substitution (DBS) Signatures (DBS: `cosmic3_DBS`)
* Small Insertion & Deletion (ID) Signatures (ID: `cosmic3_ID`)

Use:
```{bash}
--cosmic {cosmic2,cosmic3,cosmic3_exome,cosmic3_DBS,cosmic3_ID,cosmic3_TSB}
```

---

#### Prior on H & W
We generally impose an exponential (`L1`) prior on the W & H matrices for non-negative matrix factorization.

Use:
```{bash}
--prior_on_H L1 --prior_on_W L1
```

---

#### Running the method
This method may be run in two ways, from a `.maf` file or a spectra file (`.txt`, `.parquet`, `.txt.gz`, `.csv`).
* _Mutation Annotation Format_: for details on this format (`.maf`), please see [this](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) reference from NCI's Genomic Data Commons website
  * If this option is used, `signatureanalyzer` will generate a spectra using the `.maf` based on what `--cosmic` option is selected
* _Spectra_: this option is provided if the user wants to provide a pre-computed mutational spectra (ex. 96-base context; see COSMIC site or [**generating_mutational_spectra.md**](https://github.com/broadinstitute/getzlab-SignatureAnalyzer/blob/master/docs/generating_mutational_spectra.md)

Use:
```
signatureanalyzer -n 10 \
                  --cosmic cosmic3_exome \
                  --hg_build hg38.2bit \
                  --objective poisson \
                  --max_iter 30000 \
                  --prior_on_H L1 \
                  --prior_on_W L1 \
                  input.maf
```
