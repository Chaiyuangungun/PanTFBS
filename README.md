# Dependencies

software

        GMAP
  
# Usage

1、Mapping reference cds to all genomes in pan-genome

        # For each sample, run command below
        gmap_build -D . -d DB genomeN.fasta
        gmap -D . -d DB -f 2 -n 1 -t 20 ref.cds > genomeN.gff3
2、Search TF binding sites 

        python3 TFSsearch -g genomeN.fasta -a genomeN.gff3 -i geneid.file -k 0.8 -t 10 -o outfile.prefix -d jaspar_database.csv
                -g genome file
                -a GMAP gff file
                -i geneid file for search
                -k threshold (default=0.8)
                -t thread (default=10)
                -d jaspar database (default=jaspar_database.csv)
                -0 outfile prefix 
3、 number different of TF in Samples
       
        python3 DTFN.py -s Sample_type_file -i geneid.file -d jaspar_database.csv -o outfile.prefix 
        

