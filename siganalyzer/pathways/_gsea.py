from gprofiler import gprofiler
import pandas as pd
from typing import Union
import numpy as np
from tqdm import tqdm

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def gprof(
    query: Union[list, np.ndarray, pd.DataFrame],
    custom_bg: Union[list, None] = None,
    ordered_query: bool = True,
    src_filter: Union[None, str] = None,
    subtype_key: Union[None, str] = None,
    **kwargs
    ):
    """
    G-profiler wrapper for gene set enrichment analysis.
    See here: https://biit.cs.ut.ee/gprofiler/page/apis
    --------------------------------
    Gprofiler query: list of genes/proteins; can be mixed types of gene IDs, SNP IDs,
        chromosomal intervals, or term IDs

    Args:
        * query: iterable of possible Gprofiler queries or pd.DataFrame signatures output from
            ARD-NMF:

                signatures = pd.read_hdf("nmf_output.h5", "signatures")

        * custom_bg: custom background set for hypergeometric tests
            if None   --> use default background selection for gProfiler
            if (list) --> use provided genes as background
        * ordered_query (bool): g:Profiler gene lists may be interpreted as ordered
            lists where elements (genes, proteins, probesets) are in the order
            of decreasing importance. The ordered query option is useful when
            the genes can be placed in some biologically meaningful order.
            For instance, according to differential expression in a given
            RNA-seq experiment. g:Profiler then performs incremental enrichment analysis
            with increasingly larger numbers of genes starting from the top of
            the list. This procedure identifies specific functional terms that
            associate to most dramatic changes in given experimental setup, as
            well as broader terms that characterise the gene set as a whole;
            True <default>
        * src_filter (list): list of the following database sources:
            Source ID:
                GO:MF	molecular function
                GO:CC	cellular component
                GO:BP	biological process
                KEGG	Kyoto Encyclopedia of Genes and Genomes
                REAC	Reactome
                WP	    WikiPathways
                TF	    Transfac
                MIRNA	miRTarBase
                HPA	    Human Protein Atlas
                CORUM	CORUM protein complexes
                HP	    Human Phenotype Ontology
        * subtype_key (str): column to use to split input markers by; for default
            output from ARD-NMF, this will be 'max_id'; NOTE: only relevant for
            running method with pd.DataFrame

    Kwargs:
        * organism (str): 'hsapiens' <default>
        * signifcant (bool): True <default>
        * exclude_iea (bool): include electronic GO annotations: annotations assigned to
            genes using in silico curation methods adn will ahve an IEA evidence
            code to back this up; ecluding IEA annotations may help reduce bias towards
            abundant and ubiquitous housekeeping genes like ribosomal genes;
            False <default>
        * max_p_value (float): 1.0 <default>
        * max_set_size (int): 0 <default>
        * correction_method (str): 'analytical': gcustom g:SCS for reducing signifance -
            traditional multiple tests corrections are for tests independent of one another
            which is not correct for GO analysis; 'fdr': benjamini-hochberg, 'bc': bonferroni correction; 'analytical' <default>

    Returns:
        * pd.DataFrame of enrichment results
            ** Over all groupings if passed DataFrame
    """
    AVAIL_DATABASES = ("GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC", "WP", "TF", "MIRNA", "HPA", "CORUM", "HP")

    if src_filter is not None:
        for db in src_filter:
            assert db in AVAIL_DATABASES, "{} is not an available database.".format(db)

    # Gprofiler syntax
    if custom_bg is None:
        custom_bg = list()

    # ----------------------------
    # Bare-bones method
    # ----------------------------
    if isinstance(query, np.ndarray) or isinstance(query, list):
        enrichment =  gprofiler(
            query,
            custom_bg=custom_bg,
            ordered_query=ordered_query,
            src_filter=src_filter,
            **kwargs
        )

        if enrichment is not None:
            return enrichment.loc[:,['p.value','term.size','overlap.size','term.name','domain','intersection']].sort_values('p.value')
        else:
            print("   * no overlaps found.")

    # ----------------------------
    # Wrapper for ARD-NMF output
    # ----------------------------
    elif isinstance(query, pd.DataFrame):
        assert subtype_key in query, "{} NOT a column in input query. Set *subtype_key* to a valid column name.".format(subtype_key)

        enrichments = list()
        for cat in tqdm(query[subtype_key].astype("category").cat.categories, desc="Computing overlaps"):
            enrichment = gprofiler(
                query[query[subtype_key]==cat].index,
                custom_bg=custom_bg,
                ordered_query=ordered_query,
                src_filter=src_filter,
                **kwargs
            )

            if enrichment is not None:
                enrichment['subtype'] = cat
                enrichment = enrichment.loc[:,['subtype','p.value','term.size','overlap.size','term.name','domain','intersection']].sort_values('p.value')
                enrichments.append(enrichment)
            else:
                print("   * no overlaps found in {}.".format(cat))

        return pd.concat(enrichments)
