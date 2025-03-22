
nextflow.enable.dsl=2

process FASTQC {
    
  container 'biocontainers/fastqc:v0.11.9_cv8'
  publishDir "results"

  input:
  path infilename

  output:
  path '*html'
  path '*zip', emit: zip
      
  script:
  """
  fastqc ${infilename}  
  """
}

workflow {

  // parse input:
  infile_channel = Channel.fromPath( params.infile )
                    .view()
  // run FASTQC:
  out_channel = FASTQC(infile_channel)
  println ("out_html: ") 
  out_channel.html.view()
 // println ("out_zip: ") 
 // out_channel.zip.view()


}
