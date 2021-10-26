# SubDiv
A pipeline for subgenome dividing

# Purpose
SubDiv pipeline is used for subgenome division within polyploid genomes

# Features
The using of this pipeline is very sample and easy, and the running time in the test is less than one day.

The mutiple CPUs running is supported and suggested, which enable the analysis completed within several hours.

# Install

The SubDiv need not install processing and users who want to use this pipeline just need to download the three scripts (two Python scripts and one R script) to your workplace or somewhere else.

# Usage

Before you use this pipeline for subgenome dividing, make sure 'makeblastdb' and 'blastn' can be find in your envionmental path. If not, please add them to your path by

export PATH=/Path/to/your/makeblastdb-blastn:$PATH

The R package 'factoextra' had been installed in your R project and Rscript command should be able to use in the window. If not, please add Rscript path to your environment by

export PATH=/Path/to/your/Rscript:$PATH

The flowing of running as the like as the examples of "work.example.sh"

###################

chromosome_scale_genome='XXXXXXXXXX'

path='path/to/bin/'

python3 $path/identify_homo_chrs.py $chromosome_scale_genome 250 250 30

python3 $path/divide_distinct_kmers_from_homochrs.py $chromosome_scale_genome pairs_res_file 30 13 2

Rscript $path/clustering_chrs.R distictive_kmer_and_counts cluster_center.pdf dendrogram.pdf

###################

The first two lines indicate the input chromosome-scale genome and bin path

The third line indicate the step1 (homoeologous chromosome pair identification) of running that includes four parameters which is (1) input genome (2) step length of break short-length segments (3) window length of break short-length segments (4) CPU numbers

The fourth line indicate the step2 (subgenome-specific repeat K-mer identification and counts) of running that includes four parameters which in order is (1) input genome (2) CPU numbers (3) K-mer length (11,13 or 15 is suggested) (4) the lowest different times of subgenome-specific repeat K-mer counts within each homoeologous chromosome pair

The fifth line indicate the step3 of (Normalization of counts and clustering of chromosomes within each subgenome), which uses R language to perform the clustering analysis. This step requires the result of the identified  subgenome-specific repeat K-mer couts matix as input (the first parameter), and the other two parameter give the figures of chromosomes cluster result within subgenomes.
