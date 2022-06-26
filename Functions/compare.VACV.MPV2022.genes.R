compare.VACV.MPV2022.genes = function(GENE){
  
  # DNA_STRING <-str_sub(string = GENE, start = -1)
  # 
  # if (DNA_STRING %in% c("R", "L")) {
  #   
  #   VACV_gb <- genbankr::readGenBank(file = here("Data", "Refseqs", "NC_006998.gb"))
  #   
  #   VACV_genes <- VACV_gb %>% genbankr::genes()
  #   
  #   genes_IR <- (VACV_genes %>% ranges())
  #   
  #   VACV_genome <- VACV_gb %>% genbankr::getSeq()
  #   
  #   prot_gene <- VACV_genome$`WR (Western Reserve)`[genes_IR[which(VACV_genes$gene == GENE)]]
  #   
  #   Biostrings::writeXStringSet(x = prot_gene %>% Biostrings::DNAStringSet(), 
  #                               filepath = here("Data", "2022", "prot_gene.fasta"))
  #   
  #   setwd(here("Data", "2022"))
  #   system("mafft --thread -1 --add prot_gene.fasta out.fasta > out_prot_gene_msa.fasta")
  #   setwd(here())
  #   
  #   out_prot_gene_msa <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "out_prot_gene_msa.fasta"), 
  #                                                     format = "fasta") %>% toupper() %>% Biostrings::DNAStringSet()
  #   
  #   prot_gene_aligned <- out_prot_gene_msa %>% tail(1) %>% Biostrings::DNAMultipleAlignment()
  #   
  #   prot_gene_aligned <- prot_gene_aligned %>% Biostrings::maskGaps() 
  #   
  #   gaps <- prot_gene_aligned %>% Biostrings::colmask()
  #   
  #   out_prot_gene_msa_aligned <- out_prot_gene_msa %>% Biostrings::DNAMultipleAlignment()
  #   
  #   Biostrings::colmask(out_prot_gene_msa_aligned) <- gaps
  #   
  #   msa_gene <- out_prot_gene_msa_aligned %>% as("DNAStringSet")
  #   
  #   msa_prots <- Biostrings::AAStringSet()
  #   
  #   for (i in seq_along(msa_gene)) {
  #     temp <- msa_gene[i] %>% str_remove_all(pattern = "-") %>% Biostrings::DNAString()
  #     
  #     if(DNA_STRING == "L"){
  #       msa_prots[[i]] <- Biostrings::reverseComplement(x = temp) %>% Biostrings::translate(if.fuzzy.codon = "solve")
  #     }
  #     
  #     if(DNA_STRING == "R"){
  #       msa_prots[[i]] <- temp %>% Biostrings::translate(if.fuzzy.codon = "solve")
  #     }
  #     
  #   }
  #   
  #   names(msa_prots) <- names(msa_gene)
  #   
  #   Biostrings::writeXStringSet(x = msa_prots, 
  #                               filepath = here("Data", "2022", "msa_prots.fasta"))
  #   
  #   setwd(here("Data", "2022"))
  #   system(paste0("mafft --thread -1 msa_prots.fasta > ", GENE, "_msa.fasta"))
  #   setwd(here())
  #   
  #   print(paste0("Done with ", GENE))
  #   
  # }
  # else{
  #   print("Protein not found!")
  # }
  
  prot_all <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(GENE, "_msa.fasta")), 
                                          format = "fasta") 
  
  prot <- c(prot_all[str_which(names(prot_all), "NC_006998")] %>% str_remove_all("-"),
            prot_all[1:(length(prot_all) - 3)] %>% str_remove_all("-"))
  
  num_subs = vector(mode = "numeric")
  num_indels = vector(mode = "numeric")

  for (i in 2:length(prot)) {
    
    paired <- Biostrings::pairwiseAlignment(pattern = prot[i], subject = prot[1])
    
    num_subs <- c(num_subs, paired %>% Biostrings::nmismatch())
    
    dels <- paired %>% Biostrings::deletion()
    
    ins <- paired %>% Biostrings::insertion()
    
    num_indels <- c(num_indels, ins[[1]] %>% length + dels[[1]] %>% length)
  }
  
  t <- Biostrings::width(prot[1]) 
  
  return(
    data.frame(Protein = GENE,
               num_subs = num_subs,
               num_indels = num_indels,
               width = t, stringsAsFactors = FALSE)
    )
  
  
  
  
  
  # prot <- Biostrings::readAAStringSet(filepath = here("Data", "2022", paste0(GENE, "_msa.fasta")), 
  #                                     format = "fasta") 
  # 
  # # this reports the number of differences between each pair of sequences, ignoring any gap position
  # dists_SNP_prot <- ape::dist.dna(x = ape::as.AAbin(as.matrix(prot[-1])), 
  #                                 model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
  # 
  # 
  # a <- dists_SNP_prot[str_which(colnames(dists_SNP_prot), "COP"), 
  #                       str_which(colnames(dists_SNP_prot), "Vaccinia")]
  # 
  # b <- dists_SNP_prot[str_which(colnames(dists_SNP_prot), "NC_003310.1 Monkeypox"), 
  #                       str_which(colnames(dists_SNP_prot), "Vaccinia")]
  # 
  # c <- dists_SNP_prot[str_which(colnames(dists_SNP_prot), "NC_003310.1 Monkeypox"), 
  #                       str_which(colnames(dists_SNP_prot), "COP")]
  # 
  # t <- Biostrings::width(prot) %>% unique()
  # 
  # d <- dists_SNP_prot[1:(nrow(dists_SNP_prot)-3), str_which(colnames(dists_SNP_prot), "Vaccinia")] %>% summary()
  # 
  # e <- dists_SNP_prot[1:(nrow(dists_SNP_prot)-3), str_which(colnames(dists_SNP_prot), "COP")] %>% summary()
  # 
  # f <- dists_SNP_prot[1:(nrow(dists_SNP_prot)-3), str_which(colnames(dists_SNP_prot), "NC_003310.1 Monkeypox")] %>% summary()
  # 
  # temp_SNPs  <- matrix(data = list(NA, b/t*100, a/t*100, d/t*100, b, NA, c/t*100, f/t*100, a, c, NA, e/t*100), ncol = 4, byrow = TRUE)
  # 
  # colnames(temp_SNPs) <- c("VACV", "MPV-CA", "MPV-WA", "MPV-2022")
  # 
  # rownames(temp_SNPs) <- c("VACV", "MPV-CA", "MPV-WA")
  # 
  # dists_INDEL_prot <- ape::dist.dna(x = ape::as.AAbin(as.matrix(prot)), 
  #                                   model = "indel", as.matrix = TRUE)
  # 
  # a <- dists_INDEL_prot[str_which(colnames(dists_INDEL_prot), "COP"), 
  #                       str_which(colnames(dists_INDEL_prot), "Vaccinia")]
  # 
  # b <- dists_INDEL_prot[str_which(colnames(dists_INDEL_prot), "NC_003310.1 Monkeypox"), 
  #                       str_which(colnames(dists_INDEL_prot), "Vaccinia")]
  # 
  # c <- dists_INDEL_prot[str_which(colnames(dists_INDEL_prot), "NC_003310.1 Monkeypox"), 
  #                       str_which(colnames(dists_INDEL_prot), "COP")]
  # 
  # d <- dists_INDEL_prot[1:(nrow(dists_INDEL_prot)-3), str_which(colnames(dists_INDEL_prot), "Vaccinia")] %>% summary()
  # 
  # e <- dists_INDEL_prot[1:(nrow(dists_INDEL_prot)-3), str_which(colnames(dists_INDEL_prot), "COP")] %>% summary()
  # 
  # f <- dists_INDEL_prot[1:(nrow(dists_INDEL_prot)-3), str_which(colnames(dists_INDEL_prot), "NC_003310.1 Monkeypox")] %>% summary()
  # 
  # temp_INDELs  <- matrix(data = list(NA, b/t*100, a/t*100, d/t*100, b, NA, c/t*100, f/t*100, a, c, NA, e/t*100), ncol = 4, byrow = TRUE)
  # 
  # colnames(temp_INDELs) <- c("VACV", "MPV-CA", "MPV-WA", "MPV-2022")
  # 
  # rownames(temp_INDELs) <- c("VACV", "MPV-CA", "MPV-WA")
  # 
  # return(list(all.SNPs.comp = temp_SNPs,
  #             all.INDELs.comp = temp_INDELs,
  #             msa.width = t))
  
}