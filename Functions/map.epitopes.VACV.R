map.epitopes.VACV = function(GENE, epitopes){

  VACV_gb <- genbankr::readGenBank(file = here("Data", "Refseqs", "NC_006998.gb"))
  
  VACV_transcripts <- VACV_gb %>% genbankr::transcripts()
  
  VACV_genome <- VACV_gb %>% genbankr::getSeq()
  
  prot_gene <- VACV_transcripts$translation[which(VACV_transcripts$gene == GENE)] %>% as.character()

  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot_gene, NUM_MAX_MISMATCH = 0)
  
  out$Description <- VACV_transcripts$product[which(VACV_transcripts$gene == GENE)]

out

}