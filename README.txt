Anti-CRISPR Gene Prediction

Developed a pipeline that predicts putative anti-CRISPR genes by utilizing tools that have been developed for 
predicting CRISPR-Cas systems and prophage regions. 


Installation
For the pipeline to run, it requires four tools to be installed:
(i) Virsorter2: https://github.com/jiarong/VirSorter2
(ii) CRISPRone: https://github.com/mgtools/CRISPRone
(iii) AcRanker: https://github.com/amina01/AcRanker 
(iv) BLAST: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

These tools can be installed seperately and their path can be added in the pipeline or you can just run the pipline as the paths for these tools
are already present in the pipeline. 

Few points to remember
(i) Virsorter2 requires its own conda environment. It is better to run your files on virsorter2 and get the phage regions before using the pipeline.

(ii) AcRanker requires packages with specific versions for the tool to work.They are:
    scikit-learn==0.23.1
    scipy==1.5.2
    xgboost==0.90


Usage
python Final_pipeline.py –f [Genome.fna]  -b [Outdir] -tsv  [Virsorter2 result.out/final-viral-boundary.tsv]

-f: Path and name of the input genome file
-b: Output directory
-tsv: Directory to Virsorter2 result- final-viral-boundary.tsv


Test Dataset
Genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/508/765/GCA_000508765.1_ASM50876v1/GCA_000508765.1_ASM50876v1_genomic.fna.gz
Command: python /home/sithata/Prediction_pipeline/Acr_prediction.py –f  GCA_000508765.1_ASM50876v1_genomic.fna -b [Outdir] -tsv  GCA_000508765.1_ASM50876v1/final-viral-boundary.tsv


