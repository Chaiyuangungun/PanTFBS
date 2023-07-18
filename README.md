# Dependencies

software

        GMAP
  
# Usage

1、Mapping reference cds to all genomes in pan-genome

        # For each sample, run command below
        gmap_build -D . -d DB genomeN.fasta
        gmap -D . -d DB -f 2 -n 1 -t 20 ref.cds > genomeN.gff3
