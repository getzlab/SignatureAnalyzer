import itertools
import numpy as np
import pandas as pd
from twobitreader import TwoBitFile
from typing import Union
from tqdm import tqdm
from sys import stdout
from .utils import compl, get_true_snps_from_maf, get_dnps_from_maf

acontext = itertools.product('A', 'CGT', 'ACGT', 'ACGT')
ccontext = itertools.product('C', 'AGT', 'ACGT', 'ACGT')

context96 = dict(zip(map(''.join, itertools.chain(acontext, ccontext)), range(1, 97)))
context78 = dict(zip(['AC>CA', 'AC>CG', 'AC>CT', 'AC>GA', 'AC>GG', 'AC>GT', 'AC>TA', 'AC>TG', 'AC>TT', 'AT>CA',
                      'AT>CC', 'AT>CG', 'AT>GA', 'AT>GC', 'AT>TA', 'CC>AA', 'CC>AG', 'CC>AT', 'CC>GA', 'CC>GG',
                      'CC>GT', 'CC>TA', 'CC>TG', 'CC>TT', 'CG>AT', 'CG>GC', 'CG>GT', 'CG>TA', 'CG>TC', 'CG>TT',
                      'CT>AA', 'CT>AC', 'CT>AG', 'CT>GA', 'CT>GC', 'CT>GG', 'CT>TA', 'CT>TC', 'CT>TG', 'GC>AA',
                      'GC>AG', 'GC>AT', 'GC>CA', 'GC>CG', 'GC>TA', 'TA>AT', 'TA>CG', 'TA>CT', 'TA>GC', 'TA>GG',
                      'TA>GT', 'TC>AA', 'TC>AG', 'TC>AT', 'TC>CA', 'TC>CG', 'TC>CT', 'TC>GA', 'TC>GG', 'TC>GT',
                      'TG>AA', 'TG>AC', 'TG>AT', 'TG>CA', 'TG>CC', 'TG>CT', 'TG>GA', 'TG>GC', 'TG>GT', 'TT>AA',
                      'TT>AC', 'TT>AG', 'TT>CA', 'TT>CC', 'TT>CG', 'TT>GA', 'TT>GC', 'TT>GG'], range(1, 79)))

def get_spectra_from_maf(
    maf: pd.DataFrame,
    hgfile: Union[str,None] = None,
    cosmic: str = 'cosmic2',
    real_snps: bool = False
    ):
    """
    Attaches context categories to maf and gets counts of contexts for each sample
    ---------------------------
    Args:
        * maf: Pandas DataFrame of maf
        * hgfile: path to 2bit genome build file for computing reference context
        * cosmic: cosmic signatures to decompose to

    Returns:
        * Pandas DataFrame of maf with context category attached
        * Pandas DataFrame of counts with samples as columns and context as rows
    """
    maf = maf.copy()

    if 'Start_Position' in list(maf):
        maf = maf.rename(columns={'Start_Position':'Start_position'})

    maf['sample'] = maf['Tumor_Sample_Barcode']

    if cosmic in ['cosmic2', 'cosmic3', 'cosmic3_exome']:
        # Subset to SNPs
        if 'Variant_Type' in maf.columns:
            maf = maf.loc[maf['Variant_Type'] == 'SNP']
        else:
            maf = maf.loc[maf['Reference_Allele'].apply(lambda k: len(k) == 1 and k != '-') & \
            maf['Tumor_Seq_Allele2'].apply(lambda k: len(k) == 1 and k != '-')]
        if not real_snps:
            maf = get_true_snps_from_maf(maf)

        ref = maf['Reference_Allele'].str.upper()
        alt = maf['Tumor_Seq_Allele2'].str.upper()

        if 'ref_context' in list(maf):
            context = maf['ref_context'].str.upper()
        else:
            assert hgfile is not None, 'Please provide genome build file.'
            hg = TwoBitFile(hgfile)

            # Map contexts
            _contexts = list()
            maf_size = maf.shape[0]
            for idx,(pos,chromosome) in enumerate(zip(maf["Start_position"].astype(int), maf["Chromosome"].astype(str))):
                stdout.write("\r      * Mapping contexts: {} / {}".format(idx, maf_size))

                # Double check version
                if chromosome == '23':
                    chromosome = 'X'
                elif chromosome == '24':
                    chromosome = 'Y'
                elif chromosome == 'MT':
                    chromosome = 'M'

                _contexts.append(hg['chr'+chromosome][pos-2:pos+1].lower())

            maf['ref_context'] = _contexts
            stdout.write("\n")
            context = maf['ref_context'].str.upper()

        n_context = context.str.len()
        mid = n_context // 2

        contig = pd.Series([r + a + c[m - 1] + c[m + 1] if r in 'AC' \
                            else compl(r + a + c[m + 1] + c[m - 1]) \
                            for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)

        try:
            maf['context96.num'] = contig.apply(context96.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        maf['context96.word'] = contig
        spectra = maf.groupby(['context96.word', 'sample']).size().unstack().fillna(0).astype(int)
        for c in context96:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra.sort_index(inplace=True)

    elif cosmic == 'cosmic3_DBS':
        # Subset to DNPs
        if 'Variant_Type' not in maf.columns:
            ref_alt = maf['Reference_Allele'] + '>' + maf['Tumor_Seq_Allele2']

            def get_variant_type(ra):
                r, a = ra.split('>')
                if len(r) == 1 and r != '-' and len(a) == 1 and a != '-':
                    return 'SNP'
                if len(r) == 2 and len(a) == 2:
                    return 'DNP'
            maf['Variant_Type'] = ref_alt.apply(get_variant_type)
        if 'DNP' in maf['Variant_Type']:
            maf = maf.loc[maf['Variant_Type'] == 'DNP']
        else:
            maf = get_dnps_from_maf(maf)

        ref = maf['Reference_Allele']
        alt = maf['Tumor_Seq_Allele2']

        contig = pd.Series([r + '>' + a if r + '>' + a in context78
                            else compl(r, reverse=True) + '>' + compl(a, reverse=True)
                            for r, a in zip(ref, alt)], index=maf.index)

        try:
            maf['context78.num'] = contig.apply(context78.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        maf['context78.word'] = contig
        spectra = maf.groupby(['context78.word', 'sample']).size().unstack().fillna(0).astype(int)
        for c in context78:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra.sort_index(inplace=True)

    else:

        raise NotImplementedError()

    return maf, spectra
