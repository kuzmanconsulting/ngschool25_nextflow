
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'
  tag "$sampleid"

  input:
  tuple val(sampleid),path(infiles)
  
  output:
  path '*html' 
  
  script:
  """
  fastqc ${infiles[0]} ${infiles[1]}
  """
}


workflow {

  // parse input:
  infile_channel = Channel.fromFilePairs( params.infile )
                    .view()
  // run FASTQC:
  FASTQC(infile_channel) 
}
