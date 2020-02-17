import itertools
import pandas as pd
from twobitreader import TwoBitFile
from typing import Union
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

context83 = dict(zip(['Cdel1', 'Cdel2', 'Cdel3', 'Cdel4', 'Cdel5', 'Cdel6+',
                       'Tdel1', 'Tdel2', 'Tdel3', 'Tdel4', 'Tdel5', 'Tdel6+',
                       'Cins0', 'Cins1', 'Cins2', 'Cins3', 'Cins4', 'Cins5+',
                       'Tins0', 'Tins1', 'Tins2', 'Tins3', 'Tins4', 'Tins5+',
                       '2del1', '2del2', '2del3', '2del4', '2del5', '2del6+',
                       '3del1', '3del2', '3del3', '3del4', '3del5', '3del6+',
                       '4del1', '4del2', '4del3', '4del4', '4del5', '4del6+',
                       '5+del1', '5+del2', '5+del3', '5+del4', '5+del5', '5+del6+',
                       '2ins0', '2ins1', '2ins2', '2ins3', '2ins4', '2ins5+',
                       '3ins0', '3ins1', '3ins2', '3ins3', '3ins4', '3ins5+',
                       '4ins0', '4ins1', '4ins2', '4ins3', '4ins4', '4ins5+',
                       '5+ins0', '5+ins1', '5+ins2', '5+ins3', '5+ins4', '5+ins5+',
                       '2delm1', '3delm1', '3delm2', '4delm1', '4delm2', '4delm3',
                       '5+delm1', '5+delm2', '5+delm3', '5+delm4', '5+delm5+'], range(1, 84)))

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

            try:
                hg = TwoBitFile(hgfile)
            except:
                raise Exception("{} not a valid 2bit file.".format(hgfile))

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
                if not chromosome.startswith('chr'):
                    chromosome = 'chr' + chromosome

                _contexts.append(hg[chromosome][pos-2:pos+1].lower())

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
        spectra = spectra.loc[context96]

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

        ref = maf['Reference_Allele'].str.upper()
        alt = maf['Tumor_Seq_Allele2'].str.upper()

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
        spectra = spectra.loc[context78]

    elif cosmic == 'cosmic3_ID':

        maf = maf.loc[(maf['Reference_Allele'] == '-') ^ (maf['Tumor_Seq_Allele2'] == '-')]

        ref = maf['Reference_Allele'].str.upper()
        alt = maf['Tumor_Seq_Allele2'].str.upper()

        assert hgfile is not None, 'Please provide genome build file.'

        try:
            hg = TwoBitFile(hgfile)
        except:
            raise Exception("{} not a valid 2bit file.".format(hgfile))

        # Map contexts
        contig = list()
        maf_size = maf.shape[0]
        for idx,(pos,chromosome,r,a) in enumerate(zip(maf["Start_position"].astype(int),
            maf["Chromosome"].astype(str), ref, alt)):
            stdout.write("\r      * Mapping contexts: {} / {}".format(idx, maf_size))

            # Double check version
            if chromosome == '23':
                chromosome = 'X'
            elif chromosome == '24':
                chromosome = 'Y'
            elif chromosome == 'MT':
                chromosome = 'M'
            if not chromosome.startswith('chr'):
                chromosome = 'chr' + chromosome

            if a == '-':
                del_len = len(r)
                _context = hg[chromosome][pos - 1 + del_len:pos - 1 + del_len * 6].upper()
                _context_list = [_context[n: n + del_len] for n in range(0, 5 * del_len, del_len)]
                n_repeats = 1
                for c in _context_list:
                    if c == r:
                        n_repeats += 1
                    else:
                        break
                microhomology = 0
                if n_repeats == 1:
                    for b1, b2 in zip(r, _context_list[0]):
                        if b1 == b2:
                            microhomology += 1
                        else:
                            break
                    prev_context = hg[chromosome][pos - 1 - del_len: pos - 1].upper()
                    for b1, b2 in zip(reversed(r), reversed(prev_context)):
                        if b1 == b2:
                            microhomology += 1
                        else:
                            break
                if del_len == 1:
                    pre = 'C' if r in 'CG' else 'T'
                elif del_len >= 5:
                    pre = '5+'
                else:
                    pre = str(del_len)
                if microhomology >= 5:
                    post = 'm5+'
                elif microhomology:
                    post = 'm' + str(microhomology)
                elif n_repeats == 6:
                    post = '6+'
                else:
                    post = str(n_repeats)
                contig.append(pre + 'del' + post)

            elif r == '-':
                ins_len = len(a)
                _context = hg[chromosome][pos:pos + ins_len * 5].upper()
                _context_list = [_context[n: n + ins_len] for n in range(0, 5 * ins_len, ins_len)]
                n_repeats = 0
                for c in _context_list:
                    if c == a:
                        n_repeats += 1
                    else:
                        break
                if ins_len == 1:
                    pre = 'C' if a in 'CG' else 'T'
                elif ins_len >= 5:
                    pre = '5+'
                else:
                    pre = str(ins_len)
                if n_repeats == 5:
                    post = '5+'
                else:
                    post = str(n_repeats)
                contig.append(pre + 'ins' + post)

        maf['context83.word'] = contig
        try:
            maf['context83.num'] = maf['context83.word'].apply(context83.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        spectra = maf.groupby(['context83.word', 'sample']).size().unstack().fillna(0).astype(int)
        for c in context83:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra = spectra.loc[context83]

        stdout.write("\n")
    else:
        raise NotImplementedError()

    return maf, spectra
