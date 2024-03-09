#!/bin/bash

## @author     Pilar Rentero Garrido
## @email      prentero@incliva.es
## @github     https://github.com/pirega
## @twitter    https://twitter.com/pirega83

output_dir="/media/scratch1/09MBIF/1_datos/3_anotacion/"
for FILE in $(ls /media/scratch1/09MBIF/1_datos/4_procesados/2_trimmed/*.fq)
do
    #echo $FILE
    ID=$(basename $FILE)
    SEQ=$(echo $ID | cut -d'_' -f1)
    SAM="$output_dir$(echo $SEQ)_hisat2.sam"
    BAM="$output_dir$(echo $SEQ)_hisat2.bam"
    SBAM="$output_dir$(echo $SEQ)_sorted.bam"
    echo "Procesando $ID"
    hisat2 -k1 -p 20 -U $FILE -x /media/scratch1/09MBIF/1_datos/2_reference_genome/genome -S $SAM
    echo "Finalizada anotacion"
    samtools view -Sbh $SAM > $BAM
    samtools sort $BAM -o $SBAM
    samtools index $SBAM
    samtools idxstats $SBAM    
done
