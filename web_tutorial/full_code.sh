#!/bin/bash

# script takes infile_1 and infile_2 fastq files and trimmes. 
# It produces fastqc on trimmed and original data. 
# usage: full_code.sh data/liver_1.fq data/liver_1.fq


infile_1=$1
infile_1=`basename $1`
infile_2=`basename $2`

outprefix=`basename ${infile_1%%_1.fq}.html ${infile_1%%_1.fq} `
datadir_short=`dirname $1`
datadir=`realpath $datadir_short`

mkdir test_results
docker run -v $datadir:/data -it biocontainers/fastp:v0.20.1_cv1 fastp --in1 ${infile_1} \
      --in2 ${infile_2} \
      --out1 ${infile_1%%.fq}.trimmed.fq \
      --out2 ${infile_2%%.fq}.trimmed.fq \
      -h test_results/${outprefix}

docker run -v .:/data -it biocontainers/fastqc:v0.11.9_cv8 fastqc ${infile_1%%.fq}.trimmed.fq -o ./test_results
docker run -v .:/data -it biocontainers/fastqc:v0.11.9_cv8 fastqc ${infile_1} -o ./test_results 
docker run -v .:/data -it biocontainers/fastqc:v0.11.9_cv8 fastqc ${infile_2} -o ./test_results
docker run -v .:/data -it biocontainers/fastqc:v0.11.9_cv8 fastqc ${infile_2%%.fq}.trimmed.fq -o ./test_results
