#!/bin/bash

## @author     Pilar Rentero Garrido
## @email      prentero@incliva.es
## @github     https://github.com/pirega
## @twitter    https://twitter.com/pirega83

output_dir="/media/scratch/3_anotacion/conteo/"
GENOME="/media/scratch1/09MBIF/1_datos/2_reference_genome/genome.gtf"
for FILE in $(ls ./*sorted.bam)
do
    #echo $FILE
    ID=$(basename $FILE)
    SEQ=$(echo $ID | cut -d'_' -f1)
    OUTPUT="$output_dir$(echo $SEQ)_counts.tsv"
    echo "Procesando $ID"
    echo "Escribiendo $OUTPUT"
    htseq-count -i gene_id -n 24 -f bam --stranded=yes -r name -s no $FILE $GENOME > $OUTPUT
    echo "Finalizado conteo"
done
