import itertools, pandas as pd

_acontext = itertools.product('A', 'CGT', 'ACGT', 'ACGT')
_ccontext = itertools.product('C', 'AGT', 'ACGT', 'ACGT')
context96 = dict(zip(map(''.join, itertools.chain(_acontext, _ccontext)), range(1, 97)))


def get_spectra_from_maf(maf):
    """
    Attaches context categories to maf and gets counts of contexts for each sample
    Args:
        maf: Pandas DataFrame of maf

    Returns:
        Pandas DataFrame of maf with context category attached
        Pandas DataFrame of counts with samples as columns and context as rows
    """
    maf = maf.copy()
    maf['sample'] = maf['Tumor_Sample_Barcode']
    if 'Variant_Type' in maf.columns:
        maf = maf.loc[maf['Variant_Type'] == 'SNP']
    else:
        maf = maf.loc[maf['Reference_Allele'].apply(lambda k: len(k) == 1 and k != '-') &
                      maf['Tumor_Seq_Allele2'].apply(lambda k: len(k) == 1 and k != '-')]
    ref = maf['Reference_Allele'].str.upper()
    alt = maf['Tumor_Seq_Allele2'].str.upper()
    context = maf['ref_context'].str.upper()
    n_context = context.str.len()
    mid = n_context // 2
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    contig = pd.Series([r + a + c[m - 1] + c[m + 1] if r in 'AC'
                        else complements[r] + complements[a] + complements[c[m - 1]] + complements[c[m + 1]]
                        for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)
    try:
        maf['context96.num'] = contig.apply(context96.__getitem__)
    except KeyError as e:
        raise KeyError('Unusual context: ' + str(e))
    maf['context96.word'] = contig
    spectra = maf.groupby(['context96.num', 'sample']).size().unstack().fillna(0).astype(int)
    return maf, spectra
