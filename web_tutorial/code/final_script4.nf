
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'
  tag "running on $sampleid"
  
  input:
  tuple val(sampleid),path(infiles)
  
  output:
  path '*html' 
  publishDir "results"

  
  script:
  """
  fastqc ${infiles[0]} ${infiles[1]}
  """
}

process TRIM {
    
  container 'biocontainers/fastp:v0.20.1_cv1'
  tag "running on $sampleid"
  
  input:
  tuple val(sampleid),path(infiles)
  
  output:
  path '*trimmed.fq' , emit: trimmed
  publishDir "results"

  
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


process FASTQC_SINGLE {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'
  tag "$sampleid"

  input:
  path infile
  
  output:
  path '*html' 
  publishDir "results"

  
  script:
  """
  fastqc ${infile}
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

}
