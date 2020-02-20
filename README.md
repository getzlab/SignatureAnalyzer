# SignatureAnalyzer

Automatic Relevance Determination (ARD) - NMF of mutational signature &amp; expression data. Designed for scalability using Pytorch to run using GPUs if available.
* See `docs` for a more in-depth description of how to use method.

_Requires Python 3.6.0 or higher._

## Installation

##### PIP

`pip3 install signatureanalyzer`

or

##### Git Clone

* `git clone --recursive https://github.com/broadinstitute/getzlab-SignatureAnalyzer.git`
* `cd getzlab-SignatureAnalyzer`
* `pip3 install -e .`

Note `--recurisve` flag is required to clone submodules.

##### Docker

Link: `http://gcr.io/broad-cga-sanand-gtex/signatureanalyzer`

* `docker pull gcr.io/broad-cga-sanand-gtex/signatureanalyzer:latest`
* `docker run -it --rm gcr.io/broad-cga-sanand-gtex/signatureanalyzer`

---

## Source Publications

**PCAWG Mutational Signatures**
* Alexandrov, L. B., Kim, J., Haradhvala, N. J., Huang, M. N., Ng, A. W. T., Wu, Y., ... & Islam, S. A. (2020). The repertoire of mutational signatures in human cancer. Nature, 578(7793), 94-101.
 * see: https://www.nature.com/articles/s41586-020-1943-3
 * see `./PCAWG/`

**SignatureAnalyzer-GPU source publication**
* Taylor-Weiner, A., Aguet, F., Haradhvala, N.J. et al. Scaling computational genomics to millions of individuals with GPUs. Genome Biol 20, 228 (2019) doi:10.1186/s13059-019-1836-7
(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1836-7)
  * see: https://github.com/broadinstitute/SignatureAnalyzer-GPU

**SignatureAnalyzer-CPU source publications**
* Kim, J. et al. Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors. Nat. Genet. 48, 600–606 (2016). (https://www.nature.com/articles/ng.3557)

* Kasar, S. et al. Whole-genome sequencing reveals activation-induced cytidine deaminase signatures during indolent chronic lymphocytic leukaemia evolution. Nat. Commun. 6, 8866 (2015). (https://www.nature.com/articles/ncomms9866)

**Mathematical details**
* Tan, V. Y. F., Edric, C.  & Evotte, F. Automatic Relevance Determination in Nonnegative Matrix Factorization with the β-Divergence. (2012). (https://arxiv.org/pdf/1111.6085.pdf)


---
## Command Line Interface

```
usage: signatureanalyzer [-h] [-t {maf,spectra,matrix}] [-n NRUNS] [-o OUTDIR]
                         [--cosmic {cosmic2,cosmic3,cosmic3_exome,cosmic3_DBS,cosmic3_ID,cosmic3_TSB}]
                         [--hg_build HG_BUILD] [--cuda_int CUDA_INT]
                         [--verbose] [--K0 K0] [--max_iter MAX_ITER]
                         [--del_ DEL_] [--tolerance TOLERANCE] [--phi PHI]
                         [--a A] [--b B] [--objective {poisson,gaussian}]
                         [--prior_on_W {L1,L2}] [--prior_on_H {L1,L2}]
                         [--report_freq REPORT_FREQ]
                         [--active_thresh ACTIVE_THRESH] [--cut_norm CUT_NORM]
                         [--cut_diff CUT_DIFF]
                         input
```

#### Example:

```
signatureanalyzer input.maf -n 10 --cosmic cosmic2 --objective poisson
```


## Python API

```python
import signatureanalyzer as sa

# ---------------------
# RUN SIGNATURE ANALYZER
# ---------------------

# Run array of decompositions with mutational signature processing
sa.run_maf(input.maf, outdir='./ardnmf_output/', cosmic='cosmic2', hg_build='./ref/hg19.2bit', nruns=10)

# Run ARD-NMF algorithm standalone
sa.ardnmf(...)

# ---------------------
# LOADING RESULTS
# ---------------------
import pandas as pd

H = pd.read_hdf('nmf_output.h5', 'H')
W = pd.read_hdf('nmf_output.h5', 'W')
Hraw = pd.read_hdf('nmf_output.h5', 'Hraw')
Wraw = pd.read_hdf('nmf_output.h5', 'Wraw')
feature_signatures = pd.read_hdf('nmf_output.h5', 'signatures')
markers = pd.read_hdf('nmf_output.h5', 'markers')
cosine = pd.read_hdf('nmf_output.h5', 'cosine')
log = pd.read_hdf('nmf_output.h5', 'log')

# Output for each run may be found at...
Hrun1 = pd.read_hdf('nmf_output.h5', 'run1/H')
Wrun1 = pd.read_hdf('nmf_output.h5', 'run1/W')
# etc...

# Aggregate output information for each run
aggr = pd.read_hdf('nmf_output.h5', 'aggr')

# ---------------------
# PLOTTING
# ---------------------
sa.pl.marker_heatmap(...)
sa.pl.signature_barplot(...)
sa.pl.stacked_bar(...)
sa.pl.k_dist(...)
sa.pl.consensus_matrix(...)

```
