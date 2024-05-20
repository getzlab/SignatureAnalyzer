from wolf import Task

class SignatureAnalyzer(Task):
    inputs = {
      "maf" : None,
      "reps" : 10,
      "type" : "maf",
      "reference" : "",
      "hg_build" : "",
      "objective" : ""
    }
    script = """
FLAGS=""
if [! -z "$reference"]
then
    FLAGS=" --referenc ${reference}"
fi
if [! -z "$hg_build"]
then
    FLAGS+=" --hg_build ${hg_build}"
fi
if [! -z "$objective"]
    FLAGS+=" --objective ${objective}"
fi
signatureanalyzer ${maf} -n ${reps} -t ${type} $FLAGS
    """
    output_patterns = {
        "cosine_plot" : "cosine_similarity_plot.pdf",
        "kdist" : "k_dist.pdf",
        "nmfout" : "nmf_output.h5",
        "contributions" : "signature_contributions.pdf",
        "stacked_bar" : "signature_stacked_barplot.pdf",
        "weighted_maf" : "signature_weighted_maf.tsv"
        }
    docker = "gcr.io/broad-cga-sanand-gtex/signatureanalyzer:latest"
    resources = { "mem" : "4G" }

