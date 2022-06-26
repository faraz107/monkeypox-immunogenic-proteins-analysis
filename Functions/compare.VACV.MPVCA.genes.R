compare.VACV.MPVCA.genes = function(PROTEIN){
  
  prot_all <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_msa.fasta")), 
                                          format = "fasta") 
  
  prot <- c(prot_all[str_which(names(prot_all), "NC_003310")] %>% str_remove_all("-"),
            prot_all[str_which(names(prot_all), "NC_006998")] %>% str_remove_all("-"))
  
  paired <- Biostrings::pairwiseAlignment(pattern = prot[2], subject = prot[1])
  
  # this reports the number of differences between each pair of sequences, ignoring any indels
  
  num_subs <- paired %>% Biostrings::nmismatch()
  
  dels <- paired %>% Biostrings::deletion()
  
  ins <- paired %>% Biostrings::insertion()
  
  num_indels <- ins[[1]] %>% length + dels[[1]] %>% length
  
  t <- Biostrings::width(prot[1]) 
  
  return(list(num_subs = num_subs,
              num_indels = num_indels,
              width = t))
  
}