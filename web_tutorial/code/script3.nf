
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'

  input:
  path infilename

  output:
  path "*html", emit: html
  path "*zip", emit: zip

  script:
  """
  fastqc ${infilename}  
  """
}

workflow {

  // parse input:
  infile_channel = Channel.fromPath( params.infile )
                          .view()
}