from wolf import Task

class SignatureAnalyzer(Task):
    inputs = {
      "maf" : None,
      "hg_build" : None,
      "type" : "maf",
      "reps" : 10,
      "reference" : "",
      "objective" : ""
    }
    script = """
FLAGS=""
if [! -z "$reference"]
then
    FLAGS=" --referenc ${reference}"
fi
if [! -z "$objective"]
    FLAGS+=" --objective ${objective}"
fi
signatureanalyzer ${maf} --hg_build ${hg_build} -n ${reps} -t ${type} $FLAGS
    """
    output_patterns = {
        "cosine_plot" : "cosine_similarity_plot.pdf",
        "kdist" : "k_dist.pdf",
        "nmfout" : "nmf_output.h5",
        "contributions" : "signature_contributions.pdf",
        "stacked_bar" : "signature_stacked_barplot.pdf",
        "weighted_maf" : "signature_weighted_maf.tsv"
        }
    docker = "gcr.io/broad-getzlab-workflows/signatureanalyzer_kprasad:v9"
    resources = { "mem" : "4G" }
