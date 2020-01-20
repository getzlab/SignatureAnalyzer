### README

This file contains a brief description of SignatureAnalzyer as used in PCAWG7 analyses.

### Installation and execution
1. All scripts and codes are self-contained and standalone, and don't have any external dependencies.
2. All scripts and codes are written in R and executed in R environments (R-3.3.3 in Mac OS X).
3. Several R libraries are required: gridExtra, ggplot2, gplots, reshape2, grid.

### Directory structure
1. `INPUT_SignatureAnalzyer`: contains all inputs (lego matrix for 96 and 1536 SNV contexts, DBSs, and INDELs) and the information on putative POLE, MSI, and single TMZ sample.
2. `OUTPUT_SignatureAnalzyer`: all final outputs will be saved in this directory.
3. `TEMPORARY_SignatureAnalzyer`: all intermediate files will be saved in this directory.
4. `OUTPUT_DEMO`: all outputs from "SignatureAnalzyer.demo.R" will be saved here.

### Main scripts and code
1. `SignatureAnalzyer.demo.R`
    * Demonstrates how SignatureAnalzyer extracts signatures for 35 PCAWG Biliary samples using 96 contexts.
    * Runtime: about 20 minutes for 10 independent BayesNMF runs with tol = 1.e-07 and K = 25 (maximum signatures) on a 2.6 GHz Intel Core i7 on a MacBook (OS X El Capitan).
2. `SignatureAnalyzer.PCAWG.COMPOSITE.R`: COMPOSITE signature extraction and activity assignment for 2780 PCAWG samples.
3. `SignatureAnalyzer.PCAWG.DNP.R`: DBS (double-base substitution) signature extraction and activity assignment for 2780 PCAWG samples.
4. `SignatureAnalyzer.PCAWG.INDEL.R`: INDEL signature extraction and activity assignment for 2780 PCAWG samples.
5. `SignatureAnalyzer.PCAWG.function.R`: contains all necessary functions.

### Description of key functions
1. `BayesNMF.L1W.L2H`
    * Bayesian non-negative matrix factorization algorithm with an exponential prior for W and a half-normal prior for H.
    * This function is used in all signature extraction in PCAWG activity.
2. `BayesNMF.L1.KL.fixed_W.Z`
    * Adopted from the Bayesian non-negative matrix factorization algorithm with an exponential prior for both W and H.
    * This function is used in the activity attribution step to select an optimal set of signatures.
3. `BayesNMF.L1.KL.fixed_W.Z.sample`
    * Adopted from the Bayesian non-negative matrix factorization algorithm with an exponential prior for both W and H.
    * This function is used for the activity attribution with a set of selected signatures.
4. Detailed descriptions for all other functions are contained in each script or code.

### References
1. Tan, VY, Fevotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence. IEEE Trans Pattern Anal Mach Intell 2013;35:1592-605.
2. Kim J, et al. Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors (2016). Nat Genet 48, 600-606.
3. Kasar, S, Kim J et al. Whole-genome sequencing reveals activation-induced cytidine deaminase signatures during indolent chronic lymphocytic leukaemia evolution. Nat Commun. 6:8866 doi: 10.1038/ncomms9866 (2015).
4. P. Polak, Kim J, L. Brounstein  et al, A mutational signature reveals alterations underlying deficient homologous recombination repair in breast cancer. Nature Genetics (2017), doi:10.1038/ng.3934
5. Haradhvala NJ, Kim J, Maruvka YE et al, Distinct mutational signatures characterize concurrent loss of polymerase proofreading and mismatch repair, Nature Commun. (2018), PMID: 29717118
