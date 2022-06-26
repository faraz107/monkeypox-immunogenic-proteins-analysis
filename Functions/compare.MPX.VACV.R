compare.MPX.VACV = function(PROTEIN){
  
  prot_all <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_msa.fasta")), 
                                          format = "fasta") 
  
  Biostrings::writeXStringSet(x = prot_all[1:(length(prot_all)-3)], 
                              filepath = here("Data", "2022", "msa_prots.fasta"))
  
  setwd(here("Data", "2022"))
  system(paste0("mafft --thread -1 msa_prots.fasta > ", PROTEIN, "_2022_msa.fasta"))
  setwd(here())
  
  prot_msa <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(PROTEIN, "_2022_msa.fasta")), 
                                          format = "fasta") 
  
  prot_consensus <- Biostrings::consensusString(prot_msa)
  
  prot <- c(prot_consensus %>% str_remove_all("-"),
            prot_all[str_which(names(prot_all), "NC_006998")] %>% str_remove_all("-"))
  
  paired <- Biostrings::pairwiseAlignment(pattern = prot[2], subject = prot[1])
  
  # this reports the number of differences between each pair of sequences, ignoring any indels
  
  num_subs <- paired %>% Biostrings::nmismatch()
  
  dels <- paired %>% Biostrings::deletion()
  
  ins <- paired %>% Biostrings::insertion()
  
  num_indels <- ins[[1]] %>% length + dels[[1]] %>% length
  
  t <- Biostrings::width(prot[2]) 
  
  pos_subs <- which(pattern(paired) %>% as.character() %>% seqinr::s2c() != subject(paired) %>% as.character() %>% seqinr::s2c())
  
  aa_refs <- (pattern(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  aa_subs <- (subject(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  Substitutions <- vector(mode = "character")
  
  for(i in seq_along(aa_refs)){
    Substitutions <- c(Substitutions,
                       paste0(aa_refs[i], pos_subs[i], aa_subs[i])) 
  }
  
  if(length(Substitutions) == 0){
    pos_subs <- NA
    aa_refs <- NA
    aa_subs <- NA
    Substitutions <- NA
  }
  
  df <- data.frame(Protein = PROTEIN,
                   Case = "MPV2022-VACV",
                   Substitutions = Substitutions,
                   Positions = paste(pos_subs, collapse = ","),
                   RefsAA = paste(aa_refs, collapse = ","),
                   SubsAA = paste(aa_subs, collapse = ","),
                   Length = t, 
                   Indels = num_indels,
                   Dissimilarity = (num_subs + num_indels)/t*100)
  
  prot <- c(prot_consensus %>% str_remove_all("-"),
            prot_all[str_which(names(prot_all), "AY753185")] %>% str_remove_all("-"))
  
  paired <- Biostrings::pairwiseAlignment(pattern = prot[2], subject = prot[1])
  
  # this reports the number of differences between each pair of sequences, ignoring any indels
  
  num_subs <- paired %>% Biostrings::nmismatch()
  
  dels <- paired %>% Biostrings::deletion()
  
  ins <- paired %>% Biostrings::insertion()
  
  num_indels <- ins[[1]] %>% length + dels[[1]] %>% length
  
  t <- Biostrings::width(prot[2]) 
  
  pos_subs <- which(pattern(paired) %>% as.character() %>% seqinr::s2c() != subject(paired) %>% as.character() %>% seqinr::s2c())
  
  aa_refs <- (pattern(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  aa_subs <- (subject(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  Substitutions <- vector(mode = "character")
  
  for(i in seq_along(aa_refs)){
    Substitutions <- c(Substitutions,
                       paste0(aa_refs[i], pos_subs[i], aa_subs[i])) 
  }
  
  if(length(Substitutions) == 0){
    pos_subs <- NA
    aa_refs <- NA
    aa_subs <- NA
    Substitutions <- NA
  }
  
  df <- df %>% bind_rows(data.frame(Protein = PROTEIN,
                                    Case = "MPV2022-MPVCA",
                                    Substitutions = Substitutions,
                                    Positions = paste(pos_subs, collapse = ","),
                                    RefsAA = paste(aa_refs, collapse = ","),
                                    SubsAA = paste(aa_subs, collapse = ","),
                                    Length = t,
                                    Indels = num_indels,
                                    Dissimilarity = (num_subs + num_indels)/t*100))
  
  prot <- c(prot_consensus %>% str_remove_all("-") %>% str_remove_all("-"),
            prot_all[str_which(names(prot_all), "NC_003310")] %>% str_remove_all("-"))
  
  paired <- Biostrings::pairwiseAlignment(pattern = prot[2], subject = prot[1])
  
  # this reports the number of differences between each pair of sequences, ignoring any indels
  
  num_subs <- paired %>% Biostrings::nmismatch()
  
  dels <- paired %>% Biostrings::deletion()
  
  ins <- paired %>% Biostrings::insertion()
  
  num_indels <- ins[[1]] %>% length + dels[[1]] %>% length
  
  t <- Biostrings::width(prot[2]) 
  
  pos_subs <- which(pattern(paired) %>% as.character() %>% seqinr::s2c() != subject(paired) %>% as.character() %>% seqinr::s2c())
  
  aa_refs <- (pattern(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  aa_subs <- (subject(paired) %>% as.character() %>% seqinr::s2c())[pos_subs]
  
  Substitutions <- vector(mode = "character")
  
  for(i in seq_along(aa_refs)){
    Substitutions <- c(Substitutions,
                       paste0(aa_refs[i], pos_subs[i], aa_subs[i])) 
  }
  
  if(length(Substitutions) == 0){
    pos_subs <- NA
    aa_refs <- NA
    aa_subs <- NA
    Substitutions <- NA
  }
  
  df <- df %>% bind_rows(data.frame(Protein = PROTEIN,
                                    Case = "MPV2022-MPVWA",
                                    Substitutions = Substitutions,
                                    Positions = paste(pos_subs, collapse = ","),
                                    RefsAA = paste(aa_refs, collapse = ","),
                                    SubsAA = paste(aa_subs, collapse = ","),
                                    Length = t,
                                    Indels = num_indels,
                                    Dissimilarity = (num_subs + num_indels)/t*100))
  
  df
  
}