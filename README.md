# Finding-TF-in-genome
Find TF binding sites for given gene_id_file from jaspar database based on genome and gff file。

python3 find_TF_sites.py -database [jasper database file] -genome [genome_file] -gff [genome_gff_file] -id [find gene_id file] -thread [numbers of thread] -threshold [similarity threshold] -out [out_file_name_prefix] 
(database、genome、gff、id is required，thread、threshold、out----Default is 10、0.9、out_file）

out_file:
gene_id1 TF_id1 TF_id2 TF_id3 ...
TF_id1 : start end direction（+、-、+R、-R:（+）5‘-3’、（+）3‘-5’、（-）3‘-5’、（-）5‘-3’） seq similarity
.
.
.
gene_id2 TF_id1 TF_id2 TF_id3 ...
.
.
.
