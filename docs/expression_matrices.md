## Expression Matrices

Decomposition of expression signatures using  `signatureanalyzer`. For identifying _de novo_ signatures in expression matrices (ex. single-cell RNA-Seq, bulk RNA-Seq, etc.). The following document is a reference of important considerations when running this method for these data-types.

---
#### Feature Selection
With most clustering methods on transcriptional data, we recommend a highly variable gene selection step. These are well documented for single-cell and bulk RNA-Seq data and may reduce the input feature space to a few thousand genes of interest.

#### Objective Function
For mutational signatures, we assume a gaussian distribution of normalized counts and use Fevotte & Tan's derivation of a gaussian objective function for ARD-NMF. Thus, it is important to use the `gaussian` objective function. We recommend the following normalization methods:
* Bulk RNA-seq: `log2(TPM+1)`; normalize TPMs with DESeq2 size factors
* Single-cell RNA-seq: `ln(CP10K+1)`; normalize CP10K w/ scran size factors

Use:
```{bash}
--objective gaussian
```

---

#### Type of Run
This specifies whether or not to do Cosmic mapping for your dataset. For expression matrices, this is not relevent. 

Use:
```{bash}
-t matrix
```

#### Prior on H & W
We generally impose an exponential (`L1`) or half-normal (`L2`) prior on the W & H matrices for non-negative matrix factorization.

Use:
```{bash}
--prior_on_H L1 --prior_on_W L1

# or

--prior_on_H L2 --prior_on_W L2
```

---

#### Running the method
This method may be run using an input of (n samples x m features).

Use:
```
signatureanalyzer -n 10 \
                  -t matrix \
                  --objective gaussian \
                  --max_iter 30000 \
                  --prior_on_H L1 \
                  --prior_on_W L1 \
                  input.tsv
```
