import itertools
import pandas as pd
from twobitreader import TwoBitFile
from typing import Union
from sys import stdout
import numpy as np
from .utils import compl, get_true_snps_from_maf, get_dnps_from_maf
from .context import context96, context1536, context78, context83, context_composite, context_polymerase, context_polymerase_id

def get_spectra_from_maf(
    maf: pd.DataFrame,
    hgfile: Union[str,None] = None,
    reference: str = 'cosmic2',
    real_snps: bool = False,
    ):
    """
    Attaches context categories to maf and gets counts of contexts for each sample
    ---------------------------
    Args:
        * maf: Pandas DataFrame of maf
        * hgfile: path to 2bit genome build file for computing reference context
        * ref: reference signatures to decompose to

    Returns:
        * Pandas DataFrame of maf with context category attached
        * Pandas DataFrame of counts with samples as columns and context as rows
    """
    maf = maf.copy()

    if 'Start_Position' in list(maf):
        maf = maf.rename(columns={'Start_Position':'Start_position'})

    maf['sample'] = maf['Tumor_Sample_Barcode']

    if reference in ['cosmic2', 'cosmic3', 'cosmic3_exome', 'pcawg_SBS']:
        # Context type
        if reference in ['cosmic2', 'cosmic3', 'cosmic3_exome']: context_num, context_form, context_use = 'context96.num', 'context96.word', context96
        else: context_num, context_form, context_use = 'context1536.num', 'context1536.arrow', context1536
        
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

            chr_contig = all('chr{}'.format(i) in hg for i in list(range(1, 23)) + ['X', 'Y'])

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
                if chr_contig and not chromosome.startswith('chr'):
                    chromosome = 'chr' + chromosome
                if not chr_contig and chromosome.startswith('chr'):
                    chromosome = chromosome[3:]

                # 96 context, get reference [pos-1, pos, pos+1]
                if reference != 'pcawg_SBS':
                    _contexts.append(hg[chromosome][pos-2:pos+1].lower())
                # 1536 context, get refernece [pos-2, pos-1, pos, pos+1, pos+2]
                else:
                    _contexts.append(hg[chromosome][pos-3:pos+2].lower())

            maf['ref_context'] = _contexts
            stdout.write("\n")
            context = maf['ref_context'].str.upper()

        n_context = context.str.len()
        mid = n_context // 2

        if reference != 'pcawg_SBS':
            contig = pd.Series([r + a + c[m - 1] + c[m + 1] if r in 'AC' \
                                else compl(r + a + c[m + 1] + c[m - 1]) \
                                for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)
        else:
            contig = pd.Series([c[m-2:m] + "[" + r + ">" + a + "]" + c[m+1:m+3] if r in 'TC' \
                                else compl(c[::-1][m-2:m] + "[" + r + ">" + a + "]" + c[::-1][m+1:m+3]) \
                                for r, a, c, m in zip(ref, alt, context, mid)], index=maf.index)
        # Common bug: MAF has mutations with malformed sequence contexts. Rather than throwing an error, just print warning
        exclude_l = [c for c in contig if c not in context_use]
        print(f"WARNING: Dropping {len(exclude_l)} / {maf.shape[0]} contexts that are not included in the reference. Ensure proper formating of MAF")
        contig = contig[contig.isin(context_use)]
        maf = maf.loc[contig.index]
        try:
            maf[context_num] = contig.apply(context_use.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        
        maf[context_form] = contig
        spectra = maf.groupby([context_form, 'sample']).size().unstack().fillna(0).astype(int)
        for c in context_use:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra = spectra.loc[context_use.keys()]   
            
    elif reference == 'cosmic3_DBS':
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
           maf['context78.num'] = contig.apply(context78.__getitem__)
        except KeyError as e:
            raise KeyError('Unusual context: ' + str(e))

        maf['context78.word'] = contig
        spectra = maf.groupby(['context78.word', 'sample']).size().unstack().fillna(0).astype(int)
        for c in context78:
            if c not in spectra.index:
                spectra.loc[c] = 0
        spectra = spectra.loc[context78.keys()]     

    elif reference == 'cosmic3_ID':

        maf = maf.loc[(maf['Reference_Allele'] == '-') ^ (maf['Tumor_Seq_Allele2'] == '-')]

        ref = maf['Reference_Allele'].str.upper()
        alt = maf['Tumor_Seq_Allele2'].str.upper()

        assert hgfile is not None, 'Please provide genome build file.'

        try:
            hg = TwoBitFile(hgfile)
        except:
            raise Exception("{} not a valid 2bit file.".format(hgfile))

        chr_contig = all('chr{}'.format(i) in hg for i in list(range(1, 23)) + ['X', 'Y'])

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
            if chr_contig and not chromosome.startswith('chr'):
                chromosome = 'chr' + chromosome
            if not chr_contig and chromosome.startswith('chr'):
                chromosome = chromosome[3:]

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
        spectra = spectra.loc[context83.keys()]

        stdout.write("\n")
    
    elif reference in ["pcawg_COMPOSITE","pcawg_COMPOSITE96"]:
        """
        Concatenate 1536 or 96 SBS, DBS, and ID spectra
        """
        maf_dbs,dbs_df = get_spectra_from_maf(maf,hgfile, 'cosmic3_DBS')
        maf_id,id_df = get_spectra_from_maf(maf,hgfile,'cosmic3_ID')
        if reference == "pcawg_COMPOSITE":
            maf_sbs,sbs_df = get_spectra_from_maf(maf,hgfile,'pcawg_SBS',real_snps)
            maf = pd.concat([maf_sbs,maf_dbs,maf_id])
            maf['context.pcawg'] = maf['context1536.arrow'].fillna('') + maf['context78.word'].fillna('') +  maf['context83.word'].fillna('')            
        else:
            maf_sbs,sbs_df = get_spectra_from_maf(maf,hgfile,'cosmic3_exome',real_snps)
            maf = pd.concat([maf_sbs,maf_dbs,maf_id])
            maf['context.pcawg'] = maf['context96.word'].fillna('') + maf['context78.word'].fillna('') +  maf['context83.word'].fillna('')            
            
        # concatenate spectra
        spectra = pd.concat([sbs_df,dbs_df,id_df]).fillna(0)
        spectra.index.name = "context.pcawg"
    elif reference in ["pcawg_SBS_ID","pcawg_SBS96_ID"]:
        """
        Concatenate 1536 or 96 SBS + ID spectra
        """
        maf_id,id_df = get_spectra_from_maf(maf,hgfile,'cosmic3_ID')
        if reference == "pcawg_SBS96_ID":
            maf_sbs, sbs_df = get_spectra_from_maf(maf,hgfile,'cosmic3_exome',real_snps)
            maf = pd.concat([maf_sbs, maf_id])
            maf['context.pcawg'] = maf['context96.word'].fillna('') + maf['context83.word'].fillna('')
        else:
            maf_sbs, sbs_df = get_spectra_from_maf(maf,hgfile,'pcawg_SBS',real_snps)
            maf = pd.concat([maf_sbs, maf_id])
            maf['context.pcawg'] = maf['context1536.arrow'].fillna('') + maf['context83.word'].fillna('')
        spectra = pd.concat([sbs_df, id_df]).fillna(0)
        spectra.index.name = "context.pcawg"

    elif reference in ['polymerase_msi','polymerase_msi96']:
        """
        Concatenate 1536 or 96 SBS + POLE/POLD-MSI ID spectra
        """

        maf_id = maf[maf['Variant_Type'].isin(['DEL','INS'])].copy()
        def get_indel_len(x):
            if len(x) >= 4:
                return("4")
            else:
                return(str(len(x)))

        maf_id['context.polymerase'] = maf_id.apply(lambda x: '' if x['Variant_Type'] not in ['DEL','INS'] else
                                                    ('DEL' + get_indel_len(x['Reference_Allele']) if x['Variant_Type']=='DEL' else
                                                     'INS' + get_indel_len(x['Tumor_Seq_Allele2'])),1)
        id_df = maf_id.groupby(['context.polymerase','sample']).size().unstack().fillna(0).astype(int)

        if reference == "polymerase_msi":
            sbs_context = 'pcawg_SBS'
            context_form = 'context1536.arrow'
        else:
            sbs_context = 'cosmic3_exome'
            context_form = 'context96.word'

        maf_sbs, sbs_df = get_spectra_from_maf(maf, hgfile, sbs_context, real_snps)
        maf = pd.concat([maf_sbs, maf_id])
        maf['context.polymerase'] = maf[context_form].fillna('') + maf['context.polymerase'].fillna('')
        maf = maf.drop(columns=context_form)
        
        spectra = pd.concat([sbs_df, id_df]).fillna(0)
        spectra.index.name = "context.polymerase"
    else:
        raise NotImplementedError()
    return maf, spectra
