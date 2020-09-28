import itertools
import pandas as pd
from twobitreader import TwoBitFile
from typing import Union
from sys import stdout
from .utils import compl, get_true_snps_from_maf, get_dnps_from_maf
from .context import context96, context1536, context78, context83, context_composite

def get_spectra_from_maf(
    maf: pd.DataFrame,
    hgfile: Union[str,None] = None,
    cosmic: str = 'cosmic2',
    real_snps: bool = False,
    composite: bool = False
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

    # Assign values to variables based on context and whether method called from composite 
    if composite: context_num, context_form, context_use, context_sbs = 'context_composite.num', 'context_composite.arrow', context_composite, context1536
    elif cosmic in ['cosmic2', 'cosmic3', 'cosmic3_exome']: context_num, context_form, context_use, context_sbs = 'context96.num', 'context96.word', context96, context96
    elif cosmic == 'cosmic3_1536': context_num, context_form, context_use, context_sbs = 'context1536.num', 'context1536.arrow', context1536, context1536
    elif cosmic == 'cosmic3_DBS':  context_num, context_form, context_use = 'context78.num', 'context78.word', context78
    elif cosmic == 'cosmic3_ID': context_num, context_form, context_use = 'context83.num', 'context83.word', context83
    else: raise NotImplementedError()
    
    if cosmic in ['cosmic2', 'cosmic3', 'cosmic3_exome', 'cosmic3_1536']:
        # Subset to SNPs
        if 'Variant_Type' in maf.columns:
            maf = maf.loc[maf['Variant_Type'] == 'SNP']
        else:
            maf = maf.loc[maf['Reference_Allele'].apply(lambda k: len(k) == 1 and k != '-') & \
            maf['Tumor_Seq_Allele2'].apply(lambda k: len(k) == 1 and k != '-')]
        if not real_snps:
            # Filter out adjacent SNPs
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

                # 96 context, get reference [pos-1, pos, pos+1]
                if cosmic != 'cosmic3_1536':
                    _contexts.append(hg[chromosome][pos-2:pos+1].lower())
                # 1536 context, get refernece [pos-2, pos-1, pos, pos+1, pos+2]
                else:
                    _contexts.append(hg[chromosome][pos-3:pos+2].lower())

            maf['ref_context'] = _contexts
            stdout.write("\n")
            context = maf['ref_context'].str.upper()

        n_context = context.str.len()
        mid = n_context // 2

        if cosmic != 'cosmic3_1536':
            contig = pd.Series([r + a + c[m - 1] + c[m + 1] if r in 'AC' \
                                else compl(r + a + c[m + 1] + c[m - 1]) \
                                for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)
        else:
            contig = pd.Series([c[m-2:m] + "[" + r + ">" + a + "]" + c[m+1:] if r in 'TC' \
                                else compl(c[::-1][m-2:m] + "[" + r + ">" + a + "]" + c[::-1][m+1:]) \
                                for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)
        try:
            maf[context_num] = contig.apply(context_use.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        maf[context_form] = contig
        spectra = maf.groupby([context_form, 'sample']).size().unstack().fillna(0).astype(int)
        for c in context_sbs:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra = spectra.loc[context_sbs]               
            
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
        if maf['Variant_Type'].str.contains('DNP').any():
            maf = maf.loc[maf['Variant_Type'] == 'DNP']
        else:
            maf = get_dnps_from_maf(maf)
        ref = maf['Reference_Allele'].str.upper()
        alt = maf['Tumor_Seq_Allele2'].str.upper()

        contig = pd.Series([r + '>' + a if r + '>' + a in context78
                            else compl(r, reverse=True) + '>' + compl(a, reverse=True)
                            for r, a in zip(ref, alt)], index=maf.index)

        try:
           maf[context_num] = contig.apply(context_use.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        maf[context_form] = contig
        spectra = maf.groupby([context_form, 'sample']).size().unstack().fillna(0).astype(int)
        for c in context_use:
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

        maf[context_form] = contig
        try:
            maf[context_num] = maf[context_form].apply(context_use.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        spectra = maf.groupby([context_form, 'sample']).size().unstack().fillna(0).astype(int)
        for c in context83:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra = spectra.loc[context83]

        stdout.write("\n")
    
    elif cosmic == 'cosmic3_composite':
        """
        Concatenate SBS, DBS, and ID spectra
        """
        # Get spectra for 3 sections
        maf_sbs,sbs_df = get_spectra_from_maf(maf,hgfile,'cosmic3_1536',real_snps,composite=True)
        maf_dbs,dbs_df = get_spectra_from_maf(maf,hgfile, 'cosmic3_DBS',composite=True)
        maf_id,id_df = get_spectra_from_maf(maf,hgfile,'cosmic3_ID',composite=True)
        # concatenate spectra
        spectra = pd.concat([sbs_df,dbs_df,id_df]).fillna(0)
        maf = pd.concat([maf_sbs,maf_dbs,maf_id])
    else:
        raise NotImplementedError()
    return maf, spectra
