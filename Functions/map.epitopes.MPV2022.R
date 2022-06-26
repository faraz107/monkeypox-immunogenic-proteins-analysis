map.epitopes.MPV2022 = function(PROTEIN, epitopes){
  
  prot <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_msa.fasta")), 
                                      format = "fasta")
  
  Biostrings::writeXStringSet(x = prot[1:(length(prot)-4)], 
                              filepath = here("Data", "2022", "msa_prots.fasta"))
  
  setwd(here("Data", "2022"))
  system(paste0("mafft --thread -1 msa_prots.fasta > ", PROTEIN, "_2022_msa.fasta"))
  setwd(here())
  
  prot_msa <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_2022_msa.fasta")), 
                                          format = "fasta") 
  
  # prot_consensus <- Biostrings::consensusString(prot_msa)
  prot_consensus <- Biostrings::consensusString(prot_msa) %>% str_remove_all(pattern = "-")
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot_consensus, NUM_MAX_MISMATCH = 0)
  
  out
  
}