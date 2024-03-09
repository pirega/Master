#!/bin/bash
## @author     Pilar Rentero Garrido
## @email      prentero@incliva.es
## @github     https://github.com/pirega
## @twitter    https://twitter.com/pirega83

for f1 in $(ls /media/scratch1/09MBIF/1_datos/1_dato_crudo/PRJNA542148/*_1.fastq.gz)
do
    echo $f1
    ID=$(basename $f1)
    SEQ=$(echo $ID | cut -d'_' -f1)
    f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
    trim_galore --illumina -q 30 -j 8 --paired --fastqc -o /media/scratch1/09MBIF/1_datos/4_procesados/2_trimmed $f1 $f2 
done
