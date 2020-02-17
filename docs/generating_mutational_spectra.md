## Generating Mutational Spectra

Generating mutational signatures using  `signatureanalyzer`. For a comprehensive  description of mutational signatures, their relevance, and references, please see the Catalogue of Somatic Mutations in Cancer, or COSMIC, [here](https://cancer.sanger.ac.uk/cosmic/signatures). The following document is a reference for easy creation of spectra using `.mafs`.

---

#### Load Mutational Annotation File (`maf`)

For more information about the `.maf` format, see [here](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/). The following examples show generating mutational spectra from an `hg19` `.maf` using `signatureanalyzer`. First, load your `maf` as a pandas dataframe.

The following columns are **required**:
* `Hugo_Symbol`
* `Tumor_Sample_Barcode`
* `Chromosome`
* `Start_Position`
* `Reference_Allele`
* `Tumor_Seq_Allele2`
* `Variant_Type`

The following columns are **optional**:
* `ref_context`: will map context if not provided

```{python}
import signatureanalyzer as sa
import pandas as pd

maf_df = pd.read_csv(<PATH_TO_MAF>, sep='\t').loc[:,[
  'Hugo_Symbol',
  'Tumor_Sample_Barcode',
  'Chromosome',
  'Start_Position',
  'Reference_Allele',
  'Tumor_Seq_Allele2'
  ]]

print(maf_df.head())
```
|    | Hugo_Symbol   | Tumor_Sample_Barcode   |   Chromosome |   Start_Position | Reference_Allele   | Tumor_Seq_Allele2   |
|---:|:--------------|:-----------------------|-------------:|-----------------:|:-------------------|:--------------------|
|  0 | CPN1          | sample_0               |           10 |        101814119 | G                  | C                   |
|  1 | MKI67         | sample_1               |           10 |        129902901 | G                  | A                   |
|  2 | NEBL          | sample_2               |           10 |         21104601 | TTACAC             | -                   |
|  3 | RP11-445N18.7 | sample_3               |           10 |         45652518 | G                  | A                   |
|  4 | ...         | sample_4               |           10 |         50667200 | G                  | A                   |

---

#### Single-base-substitution (SBS) spectra
* This encodes the 96-base context
* **note**: two forms of this exist - either input should work
  * word: ACAG --> (REF)(MUT)(LEFT NT)(RIGHT NT)
  * arrow: A[A>C]G --> (LEFT NT)(REF)>(MUT)(RIGHT NT)
* REQUIRES a 2-bit human genome build

```
_,spectra_sbs = sa.spectra.get_spectra_from_maf(maf_df, cosmic='cosmic3_exome', hgfile='hg19.2bit')

print(spectra_sbs.head().iloc[:,:5])
```
| context96.word   |   sample_0 |   sample_1 |   sample_2 |   sample_3 |   sample_4 |
|:-----------------|-----------:|-----------:|-----------:|-----------:|-----------:|
| ACAA             |          0 |          2 |          1 |          2 |          0 |
| ACAC             |          0 |          0 |          0 |          0 |          0 |
| ACAG             |          0 |          0 |          2 |          6 |          0 |
| ACAT             |          0 |          0 |          0 |          0 |          2 |
| ...             |          2 |          1 |          2 |          5 |          1 |


---

#### Doublet-base-substitution (DBS) spectra
* This encodes the 78-base context

```
_,spectra_dbs = sa.spectra.get_spectra_from_maf(maf_df, cosmic='cosmic3_DBS')

print(spectra_dbs.head())
```

| context78.word   |   sample_0 |   sample_1 |   sample_2 |   sample_3 |   sample_4 |
|:-----------------|-----------:|-----------:|-----------:|-----------:|-----------:|
| AC>CA            |          0 |          0 |          0 |          0 |          0 |
| AC>CG            |          0 |          0 |          0 |          0 |          0 |
| AC>CT            |          0 |          0 |          0 |          0 |          0 |
| AC>GA            |          0 |          0 |          0 |          0 |          0 |
| ...            |          0 |          0 |          0 |          0 |          0 |

---

#### Insertions & Deletions (ID) spectra
* This encodes the 83-base context
* REQUIRES a 2-bit human genome build

```
_,spectra_id = sa.spectra.get_spectra_from_maf(maf_df, cosmic='cosmic3_ID', hgfile='hg19.2bit')

print(spectra_id.head().iloc[:,:5])
```

| context83.word   |   sample_0 |   sample_1 |   sample_2 |   sample_3 |   sample_4 |
|:-----------------|-----------:|-----------:|-----------:|-----------:|-----------:|
| Cdel1            |          0 |          1 |          3 |          1 |          0 |
| Cdel2            |          0 |          4 |          3 |          6 |          0 |
| Cdel3            |          2 |          2 |          2 |          9 |          0 |
| Cdel4            |          0 |          0 |          0 |          3 |          0 |
| ...            |          0 |          0 |          1 |          3 |          0 |
