
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'
  tag "$sampleid"

  input:
  tuple val(sampleid),path(infiles)
  
  output:
  path '*html' 
  publishDir params.outdir

  
  script:
  """
  fastqc ${infiles[0]} ${infiles[1]}
  """
}

process TRIM {
    
  container 'biocontainers/fastp:v0.20.1_cv1'
  tag "$sampleid"

  input:
  tuple val(sampleid),path(infiles)
  
  output:
  path '*trimmed.fq', emit: trimmed
  tuple val(sampleid), path('*trimmed.fq'), emit: tup, optional: true

  publishDir params.outdir

  
  script:
  """
  fastp\
  --in1 ${infiles[0]}\
  --in2 ${infiles[1]}\
  --out1 ${sampleid}_1.trimmed.fq\
  --out2 ${sampleid}_2.trimmed.fq\
  -h ${sampleid}.html 

  """
}

process MAP {
    
  // use appropriate container:  

  // add tag - on sampleid

  input:
  // input tuple same as in TRIM
  // another input is transcriptome file, which is a constant value 
  // (needs to be passed as a absolute path from the workflow)
  path transcriptome

  output:
  path '*' 
  // add publishDir to be provided by outdir, subfolder mapped
 
  
  script:
  """
  # index the transcriptome file

  # call bwa mem transcriptome, read1, read2 -o sampleid.sam

  """
}


workflow {

  // parse input:
  infile_channel = Channel.fromFilePairs( params.infile )
                    .view()
  // run FASTQC:
  FASTQC(infile_channel) 
  trimmed_channel = TRIM(infile_channel)
  FASTQC_SINGLE(trimmed_channel.trimmed)
  // call MAP process on trimmed reads:  

}
