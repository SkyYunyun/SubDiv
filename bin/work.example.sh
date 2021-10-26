chromosome_scale_genome='XXXXXXXXXX'

path='path/to/bin/'

python3 $path/identify_homo_chrs.py $chromosome_scale_genome 250 250 30

python3 $path/divide_distinct_kmers_from_homochrs.py $chromosome_scale_genome pairs_res_file 30 13 2

Rscript $path/clustering_chrs.R distictive_kmer_and_counts cluster_center.pdf dendrogram.pdf
