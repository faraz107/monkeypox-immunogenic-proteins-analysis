map.epitopes.MPVCA = function(PROTEIN, epitopes){

  prot_msa <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_msa.fasta")), 
                                      format = "fasta")
  
  prot <- prot_msa[str_which(names(prot_msa), "NC_003310")] %>% str_remove_all("-")
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  # prot_gene <- MPVCA_transcripts$translation[which(MPVCA_transcripts$gene == GENE)] %>% as.character()
  # 
  # out <- map.REFepitopes(Tepitopes = epitopes, REF = prot_gene, NUM_MAX_MISMATCH = 0)
  # 
  # out$Description <- MPVCA_transcripts$product[which(MPVCA_transcripts$gene == GENE)]
  
  out
  
}