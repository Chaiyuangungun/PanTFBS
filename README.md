# Dependencies

software

        GMAP
  
# Usage

1、Mapping reference cds to all genomes in pan-genome

        # For each sample, run command below
        gmap_build -D . -d DB genomeN.fasta
        gmap -D . -d DB -f 2 -n 1 -t 20 ref.cds > genomeN.gff3
2、Searching TF binding sites 

        python3 TFBS.py -g genomeN.fasta -a genomeN.gff3 -i geneid.file -k 0.8 -t 10 -o outfile.prefix -d jaspar_database.csv
                -g genome file
                -a GMAP gff file
                -i geneid file for searching
                -k threshold (default=0.8)
                -t thread (default=10)
                -d jaspar database (default=jaspar_database.csv)
                -0 outfile prefix 
3、 numerical differences of TF in sample
       
        python3 NDTF.py -s Sample_type_file -i geneid.file -d jaspar_database.csv -o outfile.prefix 
                -s Classification file for samples
                -i geneid file for searching
                -d jaspar database (default=jaspar_database.csv)
                -o outfile prefix 
input file:

        geneid.file :
        example:
          geneid1
          geneid2
          ...
          geneidn
        Sample_type_file:
        example:
         A        sample1
         A        sample2
         A        sample3        
         B        sample4
         B        sample5
         B        sample6

output file:

        out.TFids.genomesites : Transcription factor binding site location on the genome
        out.TFname.promotersites : Transcription factor binding sites are based on the location of gene start sites
        out.Pvalue : Transcription factors with numerical differences

        

