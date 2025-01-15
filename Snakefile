import pandas as pd

configfile: "config.yaml"

MANIFEST = config.get("MANIFEST", "manifest.tab")
REF_SITE = config.get("REF_SITE", "/net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/ntsm/db_source/human_sites_n10.fa")
CITES_CENTER = config.get("CITES_CENTER","/net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/ntsm/db_source/human_sites_center.txt")
PCA_MATRIX = config.get("PCA_MATRIX","/net/eichler/vol28/projects/autism_genome_assembly/nobackups/qc/ntsm/db_source/human_sites_rotationalMatrix.tsv")

df = pd.read_csv(MANIFEST, sep="\t", header=0).set_index("ID", drop=True)
ref_df = df[df["TYPE"].isin(["REF", "BOTH"])]
target_df = df[df["TYPE"].isin(["TARGET", "BOTH"])]


def find_fofn(wildcards):
    return df.loc[wildcards.id, "FOFN"]


rule all:
    input:
        expand(
            "ntsm_summary/{id}.tsv", 
            id = target_df.index,
        )

rule make_count:
    input:
        fofn = find_fofn,
    output:
        count_file = "ntsm_inv_counts/{id}.count",
    params:
        ref_site = REF_SITE,
    threads: 4,
    benchmark: "benchmark/ntsm_inv_counts/{id}.txt",
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=24,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell:
        """
        ntsmCount -t {threads} -s {params.ref_site} $(cat {input.fofn}) > {output.count_file}
        """

rule summary_target:
    input:
        ref_all = expand(
            "ntsm_inv_counts/{ref_id}.count",
            ref_id = ref_df.index,
            ),
        target_count = "ntsm_inv_counts/{id}.count"
    output:
        summary = "ntsm_summary/{id}.tsv",
    threads: 4,
    benchmark: "benchmark/ntsm_summary/{id}.txt",
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=96,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell:
        """
        ntsmEval -a -t {threads} {input.target_count} {input.ref_all} | grep "{wildcards.id}\|homConcord" | sed -e 's/ntsm_inv_counts\//g' -e 's/\.count//g' > {output.summary}
        """    

rule all_pairwise:
    input:
        all_count = expand(
            "ntsm_inv_counts/{id}.count",
            id = target_df.index,
            ),
    output:
        summary = "ntsm_summary/all_pairwise.tsv",
    threads: 4,
    benchmark: "benchmark/ntsm_summary/all_pairwise.txt",
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=96,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell:
        """
        ntsmEval -a -t {threads} {input.all_count} | sed -e 's/ntsm_inv_counts\///g' -e 's/\.count//g' > {output.summary}
        """    


