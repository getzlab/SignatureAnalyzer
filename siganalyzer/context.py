



def calc_muts_spectra(input_maf, hgfile='../../dat/hg19.2bit', out_file='sampleMutSpectra.txt'):
    """
    Calculate mutational context.
    Finds each type of mutation and creates a categorical variable. Saves this output as
    the mutational spectra file.
    """
    hg = TwoBitFile(hgfile)
    sample_muts_context = defaultdict(lambda: [0]*96)

    df = pd.read_csv(input_maf, sep='\t').loc[:,['Tumor_Sample_Barcode','Variant_Type','Chromosome','Start_position','Reference_Allele','Tumor_Seq_Allele2','ttype']]
    df = df[df['Variant_Type']=='SNP']

    for idx,row in df.iterrows():
        ref_base = row['Reference_Allele'].lower()
        new_base = row['Tumor_Seq_Allele2'].lower()

        if ref_base == '-' or new_base == '-':
            continue
        if len(ref_base) > 1 or len(new_base) > 1:
            continue

        pos = int(row['Start_position'])
        chromosome = str(row['Chromosome'])

        if chromosome == '23':
            chromosome = 'X'
        elif chromosome == '24':
            chromosome = 'Y'
        elif chromosome == 'MT':
            chromosome = 'M'

        abc = hg['chr'+chromosome][pos-2:pos+1].lower()

        if abc[1] != ref_base and ref_base != '--':
            print(abc, ref_base)
            print('non-matching reference.')
            continue

        pat = (row['Tumor_Sample_Barcode'], row['ttype'])

        try:
            sample_muts_context[pat][encode(abc, new_base)] += 1
        except:
            ## because of Ns
            print("Because of Ns")
            continue

    hdr = list(MUTATION_INDICES.items())
    hdr.sort(key=lambda x:x[1])
    index_col = ['patient', 'ttype'] + [i[0][0]+'-'+i[0][1] for i in hdr]

    df = pd.DataFrame.from_dict(sample_muts_context).T
    df = df.reset_index()
    df.columns = index_col
    df.to_csv(out_file, sep='\t',index=None)
