#!/bin/bash
## @author     Pilar Rentero Garrido
## @email      prentero@incliva.es
## @github     https://github.com/pirega
## @twitter    https://twitter.com/pirega83

for f1 in $(ls /media/scratch1/09MBIF/1_datos/1_dato_crudo/PRJNA542148/*.fastq.gz)
do
    echo $f1
    ID=$(basename $f1)
    SEQ=$(echo $ID | cut -d'_' -f1)
    trim_galore --illumina -q 30 -j 8 --fastqc -o /media/scratch1/09MBIF/1_datos/4_procesados/2_trimmed $f1 
done
