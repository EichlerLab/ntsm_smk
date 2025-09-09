import pandas as pd
import os
import glob
import re

configfile: "config.yaml"

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
MANIFEST = config.get("MANIFEST", "manifest.tab")
REF_SITE = config.get("REF_SITE", f"{SNAKEMAKE_DIR}/db_source/human_sites_n10.fa")
EXTERNAL_COUNTS_DIR = config.get("EXTERNAL_COUNTS_DIR", "")
COUNT_FILE_EXP = config.get("COUNT_FILE_EXP", "count")

df = pd.read_csv(MANIFEST, sep="\t", header=0).set_index("ID", drop=True)

def find_fofn(wildcards):
    return df.loc[wildcards.id, "FOFN"]

def get_external_counts(wildcards):
    return sorted(list(glob.iglob("count_files/external/*.count")))

wildcard_constraints:
    id = r"[^/.]+"

localrules:
    all,
    link_external_count,
    check_external_links,

rule all:
    input:
        "summary/all_pairwise.pdf"

rule get_all_count_files:
    input:
        expand("count_files/{id}.count",
            id = df.index,
        ),

rule link_external_count:
    output:
        link_external_done = temp(".external_link_done")
    threads: 1,
    run:
        external_counts = sorted(glob.glob(os.path.join(EXTERNAL_COUNTS_DIR, f"*.{COUNT_FILE_EXP}")))
        external_counts_dir = "count_files/external"
        os.makedirs(external_counts_dir, exist_ok=True)
        for external_count in external_counts:
            link_name = os.path.basename(external_count).replace(f".{COUNT_FILE_EXP}","-EXT.count")
            link_dest_path = f"count_files/external/{link_name}"
            if not os.path.lexists(link_dest_path):
                os.symlink(os.path.abspath(external_count), link_dest_path)
        open(output.link_external_done, "w").close()
        
rule make_count:
    input:
        fofn = find_fofn,
    output:
        count_file = "count_files/{id}.count",
    params:
        ref_site = REF_SITE,
    threads: 8,
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=12,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell: """
        ntsmCount -t {threads} -s {params.ref_site} $(cat {input.fofn}) > {output.count_file}
        """

rule check_external_links:
    input:
        rules.link_external_count.output.link_external_done
    output:
        check_link_done = temp(".check_external_links_done")
    threads: 1,
    shell: """ 
        touch {output.check_link_done}
    """

rule all_pairwise:
    input:
        all_counts = expand("count_files/{id}.count",
                id = df.index,
            ),
        check_external_links = rules.check_external_links.output.check_link_done,
        external_all_counts = get_external_counts,
    output:
        summary = "summary/all_pairwise.tsv",
    threads: 4,
    resources:
        mem=lambda wildcards, attempt: 4 * attempt,
        hrs=24,
    singularity:
        "docker://eichlerlab/ntsm:1.2.1",
    shell: """
        ntsmEval -a -t {threads} {input.all_counts} {input.external_all_counts} | sed -e 's/external\///g' -e 's/count_files\///g' -e 's/\.count//g' > {output.summary}
        """    

rule plot_ntsm_summary:
    input:
        summary = rules.all_pairwise.output.summary,
    output:
        summary_plot = "summary/all_pairwise.pdf",
    threads: 1,
    resources:
        mem=16,
        hrs=4,
    run:
        import numpy as np
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, Normalize

        df = pd.read_csv(input.summary, sep="\t", header=0, usecols=["sample1","sample2","score"])
        
        samples = sorted(set(df["sample1"]).union(df["sample2"]))
        
        distance_matrix = pd.DataFrame(np.inf, index=samples, columns=samples)


        for idx, row in df.iterrows():
            distance_matrix.loc[row['sample1'], row['sample2']] = row['score']
            distance_matrix.loc[row['sample2'], row['sample1']] = row['score']

        np.fill_diagonal(distance_matrix.values, 0)

        data_max = df["score"].max()
        vmax = max(4.0, data_max)

        ## color code
        cool_blue = "#3b4cc0"
        middle_color = "#f7f7f7"
        warm_orange = "#f07c52"
        warm_red= "#b40426"
        
        colors = [
            (0.0, cool_blue),
            (0.5, middle_color),
            (1.0, warm_orange),
            (max(data_max, 4.0), warm_red)
        ]


        color_stops = [(min(v / vmax, 1.0), color) for v, color in sorted(colors, key=lambda x: x[0])]

        custom_cmap = LinearSegmentedColormap.from_list("my_map", color_stops)
        norm = Normalize(vmin=0.0, vmax=vmax)

        plt.figure(figsize=(30, 24))

        sns.heatmap(distance_matrix, annot=True, fmt=".2f", cmap=custom_cmap, annot_kws={"size":15}, cbar=True, norm=norm)

        plt.title(f"PCA Distance Heatmap", fontsize=30, pad=20)
        plt.xlabel("")
        plt.ylabel("")
        plt.xticks(rotation=90, fontsize=15)
        plt.yticks(rotation=0, fontsize=15)
        
        plt.savefig(output.summary_plot, format="pdf", bbox_inches="tight")
        plt.close()
