
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'

  input:
  path infilename

  output:
  path '*' 
      
  script:
  """
  fastqc ${infilename}  
  """
}

workflow {

  // parse input:
  infile_channel = Channel.fromPath( "../../data/liver_*.fq" ).view()

  // run FASTQC:
  FASTQC(infile_channel)

}