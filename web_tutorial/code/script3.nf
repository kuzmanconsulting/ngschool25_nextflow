
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'

  input:
  path infilename

  output:
  path '*html' 
  publishDir "results"

  
  script:
  """
  fastqc ${infilename}  
  """
}

workflow {

  // parse input:
  infile_channel = Channel.fromFilePairs( params.infile )
                    .view()

}