compare.MPX.genomes = function(){

MPX_WA_ref <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "AY753185refseq_monkeypox.fasta"), 
                                           format = "fasta")

MPX_CA_ref <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "NC_003310refseq_monkeypox.fasta"), 
                                           format = "fasta")

VACV_ref <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "NC_006998_refseq_vaccinia.fasta"), 
                                         format = "fasta")


setwd(here("Data", "2022"))
system("mafft --thread -1 --add AY753185refseq_monkeypox.fasta MPX_2022_msa.fasta > out1.fasta")
system("mafft --thread -1 --add NC_003310refseq_monkeypox.fasta out1.fasta > out2.fasta")
system("mafft --thread -1 --add NC_006998_refseq_vaccinia.fasta out2.fasta > out.fasta")
setwd(here())

msa_gisaid_refseqs <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "out.fasta"), 
                                                   format = "fasta") %>% toupper() %>% Biostrings::DNAStringSet()

Biostrings::writeXStringSet(x = msa_gisaid_refseqs[1:(length(msa_gisaid_refseqs)-4)], 
                            filepath = here("Data", "2022", "msa_genomes.fasta"))

# setwd(here("Data", "2022"))
# system(paste0("mafft --thread -1 msa_genomes.fasta > ", "genomes_2022_msa.fasta"))
# setwd(here())

genomes_msa <- msa_gisaid_refseqs[1:(length(msa_gisaid_refseqs)-4)]

genomes_consensus <- Biostrings::consensusString(genomes_msa, ambiguityMap = "N") %>% Biostrings::DNAStringSet()

names(genomes_consensus) <- "MPV2022_consensus"

Biostrings::writeXStringSet(x = Biostrings::DNAStringSet(x = c(genomes_consensus, VACV_ref, MPX_CA_ref, MPX_WA_ref)), 
                            filepath = here("Data", "2022", "MPXs_VACV.fasta"))

setwd(here("Data", "2022"))
system("mafft --thread -1 MPXs_VACV.fasta > out_refs.fasta")
setwd(here())

MPXs_VACV_refs <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "out_refs.fasta"), 
                                               format = "fasta") %>% toupper() %>% Biostrings::DNAStringSet()

# this reports the number of indels between each pair of sequences

dist_indels <- ape::dist.dna(x = ape::as.DNAbin(as.matrix(MPXs_VACV_refs)), 
                             model = "indel", as.matrix = TRUE)

dist_indels_1 <- ape::dist.dna(x = ape::as.DNAbin(as.matrix(MPXs_VACV_refs)), pairwise.deletion = TRUE,
                             model = "indel", as.matrix = TRUE)

# this reports the number of indel blocks (consecutive gaps) between each pair of sequences

dist_indel_blks <- ape::dist.dna(x = ape::as.DNAbin(as.matrix(MPXs_VACV_refs)), 
                        model = "indelblock", pairwise.deletion = FALSE, as.matrix = TRUE)

# this reports the number of differences between each pair of sequences, ignoring any gap position

dist_muts <- ape::dist.dna(x = ape::as.DNAbin(as.matrix(MPXs_VACV_refs)), 
                        model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)

dist_muts_1 <- ape::dist.dna(x = ape::as.DNAbin(as.matrix(MPXs_VACV_refs)), 
                           model = "N", pairwise.deletion = FALSE, as.matrix = TRUE)

t <- MPXs_VACV_refs %>% Biostrings::width() %>% unique()


# Among the 2022 sequnces & wrt the refseqs

d <- dist_muts[1:(nrow(dist_muts)-4), str_which(colnames(dist_muts), "Vaccinia")] %>% summary()
e <- dist_muts[1:(nrow(dist_muts)-4), str_which(colnames(dist_muts), "COP")] %>% summary()
f <- dist_muts[1:(nrow(dist_muts)-4), str_which(colnames(dist_muts), "NC_003310.1 Monkeypox")] %>% summary()
t <- Biostrings::width(msa_gisaid_refseqs) %>% unique()

temp  <- matrix(data = list(NA, b/t*100, a/t*100, d/t*100, b, NA, c/t*100, f/t*100, a, c, NA, e/t*100), ncol = 4, byrow = TRUE)

colnames(temp) <- c("VACV", "MPV-CA", "MPV-WA", "MPV-2022")
rownames(temp) <- c("VACV", "MPV-CA", "MPV-WA")

GENOMIC_dists <- list()

GENOMIC_dists[[1]] <- temp


d <- dist_indels[1:(nrow(dist_indels)-4), str_which(colnames(dist_indels), "Vaccinia")] %>% summary()
e <- dist_indels[1:(nrow(dist_indels)-4), str_which(colnames(dist_indels), "COP")] %>% summary()
f <- dist_indels[1:(nrow(dist_indels)-4), str_which(colnames(dist_indels), "NC_003310.1 Monkeypox")] %>% summary()
t <- Biostrings::width(msa_gisaid_refseqs) %>% unique()

temp  <- matrix(data = list(NA, b/t*100, a/t*100, d/t*100, b, NA, c/t*100, f/t*100, a, c, NA, e/t*100), ncol = 4, byrow = TRUE)

colnames(temp) <- c("VACV", "MPV-CA", "MPV-WA", "MPV-2022")
rownames(temp) <- c("VACV", "MPV-CA", "MPV-WA")

GENOMIC_dists[[2]] <- temp

names(GENOMIC_dists) <- c("SNPs", "Indels")

}