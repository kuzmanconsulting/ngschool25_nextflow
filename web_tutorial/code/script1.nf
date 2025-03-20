
nextflow.enable.dsl=2

workflow {

  // parse input:
  infile_channel_1 = 
    Channel.fromPath( "../../data/liver_*.fq" )
    .view()


}