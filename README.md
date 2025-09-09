# ntsm_smk
Snakemake pipeline for generating pairwise PCA distance matrices and visual summaries between sequencing datasets using NTSM(https://github.com/JustinChu/ntsm).


## Requirements
- [Snakemake 7+](https://snakemake.readthedocs.io/)
- Singularity (or Apptainer) for running the ntsm Docker image
- Python libraries (for plotting step):
  - `pandas`
  - `numpy`
  - `seaborn`
  - `matplotlib`

---

## Input Data
### Manifest (`manifest.tab`)
The manifest must include:
- `ID` — Unique identifier for each dataset
- `FOFN` — File-of-filenames (list of FASTQ/FASTA files)

Example:
```
ID	FOFN
SampleA	fofn/SampleA.fofn
SampleB	fofn/SampleB.fofn
SampleC	fofn/SampleC.fofn
```


### Config (`config.yaml`)
Important keys:
- `MANIFEST`: Path to the manifest file
- `REF_SITE`: Reference sites fasta (default: `db_source/human_sites_n10.fa`)
- `CITES_CENTER`: Reference center positions
- `PCA_MATRIX`: Rotational matrix file
- `EXTERNAL_COUNTS_DIR`: Directory containing external count files (optional)
- `COUNT_FILE_EXP`: File extension for the external count files(default: `count`)

---

## Running the Pipeline
1. Edit `config.yaml` and `manifest.tab` to reflect your datasets.
2. Run Snakemake:
```bash
   ln -s /net/eichler/vol28/software/pipelines/ntsm_smk/runcluster .
   ./runcluster 30
```
