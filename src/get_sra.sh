#!/bin/bash

while read i
do
    fasterq-dump --split-files --skip-technical --outdir data --threads 8 --progress $i
	gzip data/$i*.fastq
done < $1
