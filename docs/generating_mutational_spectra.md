## Generating Mutational Spectra

Generating mutational signatures using  `signatureanalyzer`. For a comprehensive  description of mutational signatures, their relevance, and references, please see the Catalogue of Somatic Mutations in Cancer, or COSMIC, [here](https://cancer.sanger.ac.uk/cosmic/signatures). The following document is a reference for easy creation of spectra using `.mafs`.

---

Generating mutational spectra from an `hg19` `.maf` using `signatureanalyzer`:

```{python}
import signatureanalyzer as sa
import pandas as pd

maf_df = pd.read_csv(<PATH_TO_MAF>, sep='\t')

# ------------------------------------
# Single-base-substitution (SBS) spectra
# -
# This encodes the 96-base context & REQUIRES a 2-bit human genome build
# ------------------------------------
_,spectra_sbs = sa.spectra.get_spectra_from_maf(maf_df, cosmic='cosmic3_exome', hgfile='hg19.2bit')

# ------------------------------------
# Doublet-base-substitution (DBS) spectra
# -
# This encodes the 96-base context & REQUIRES a 2-bit human genome build
# ------------------------------------
_,spectra_dbs = sa.spectra.get_spectra_from_maf(maf_df, cosmic='cosmic3_DBS')
```
