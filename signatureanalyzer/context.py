import itertools

# Cartesian product for 96 SNV
acontext = itertools.product('A', 'CGT', 'ACGT', 'ACGT')
ccontext = itertools.product('C', 'AGT', 'ACGT', 'ACGT')

# Cartesian product for 1536 SNV
a_1536 = itertools.product('ACGT','ACGT','[','T','>','ACG',']','ACGT','ACGT')
c_1536 = itertools.product('ACGT','ACGT','[','C','>','AGT',']','ACGT','ACGT')

# Define dictionary for all contexts
context96 = dict(zip(map(''.join, itertools.chain(acontext, ccontext)), range(1, 97)))
context1536 = dict(zip(map(''.join, itertools.chain(a_1536, c_1536)), range(1, 1537)))
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
context_composite = {**context1536, **({k:v+1536 for k,v in context78.items()}), **({k:v+1614 for k,v in context83.items()})}
