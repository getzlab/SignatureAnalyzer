## Expression Matrices

Decomposition of expression signatures using  `signatureanalyzer`. For identifying _de novo_ signatures in expression matrices (ex. single-cell RNA-seq, bulk RNA-seq, etc.). The following document is a reference for important considerations when running this method for these data-types.

---

#### Objective Function
For mutational signatures, we assume a gaussian distribution of normalized and use Fevotte & Tan's derivation of a gaussian objective function for ARD-NMF. Thus, it is important to use the default value for the objective function (`gaussian`).
* Bulk RNA-seq: `log2(TPM+1)`; normalize TPMs with DESeq2 size factors
* Single-cell RNA-seq: `ln(CP10K+1)`; normalize CP10K w/ scran size factors

Use:
```{bash}
--objective gaussian
```

---

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
This method may be run using an input of (n x m), with n: samples, m: variables.

Use:
```
signatureanalyzer -i input.tsv \
                  -n 10 \
                  --objective gaussian \
                  --max_iter 30000 \
                  --prior_on_H L1 \
                  --prior_on_W L1
```
