#!/bin/bash

while read i
do
    fasterq-dump --split-files --skip-technical --outdir data/fun_data --threads 8 --progress $i
	gzip fun_data/$i*.fastq
done < $1
