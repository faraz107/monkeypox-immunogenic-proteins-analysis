
library(tidyverse)
library(here)
library(Biostrings)
library(magrittr)

source("~/Rprojects/monkeypox/Functions/map.REFepitopes.R")

MPV_PROTEINS <- c("L6R", "E5R", "O1L", "E1R", "I8R", "D4L", "B19R", "B7R", "B9R", "E12L", "D10L", "A9R", "Q1L", "A11L", "E13L", "A4L", "M4R", "B8R", "A35R", "M1R", "H3L", "E8L", "B6R", "A30L", "A29L", "A18L")

MPV_MSA_PROTEINS <- c("L6", "E5", "O1", "E1", "I8", "D4", "B19", "B7", "B9", "E12", "D10", "A9", "Q1", "A11", "E13", "A4", "M4", "B8", "A35", "M1", "H3", "E8", "B6", "A30", "A29", "A18")

VACV_PROTEINS <- c("J6R", "D5R", "M1L", "D1R", "I8R", "C10L", "C12L", "B6R", "B8R", "D12L", "C7L", "A8R", "O1L", "A10L", "D13L", "A3L", "L4R", "B7R", "A33R", "L1R", "H3L", "D8L", "B5R", "A28L", "A27L", "A17L")

PROTEINS <- data.frame(MPV_PROTEINS, VACV_PROTEINS, MPV_MSA_PROTEINS, stringsAsFactors = FALSE)

# msa_path <- here("Data", "UGENE_Data", "Proteins", paste0(MPV_MSA_PROTEINS[i], "_nuc_MSA_transl.aln"))
# 
# MSA <- seqinr::read.alignment(file = msa_path, format = "clustal", forceToLower = FALSE)
# 
# prot <- MSA$seq[MSA$nam == "NC_003310.1:83399-87259_Monkeypox_virus_Zaire-96-I-16,_complete_genome(translated)"]

tcell_I <- read_csv("Data/IEDB/tcell_I.csv", skip = 1)

tcell_II <- read_csv("Data/IEDB/tcell_II.csv", skip = 1)

mhc_I <- read_csv("Data/IEDB/mhc_I.csv", 
                  skip = 1)

mhc_II <- read_csv("Data/IEDB/mhc_II.csv", 
                  skip = 1)

##### T cell assay data #####

##### T cell assay class I

epitopes <- tcell_I$Description...12 %>% stringr::str_squish() %>% unique()

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_PROTEINS)) {

  prot <- here("Data", "Refseqs", paste0(MPV_PROTEINS[i], "_refprot.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
    
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_PROTEINS"))

y$Virus <- "MPV"

df_tcell_I <- y

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(VACV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(VACV_PROTEINS[i], "_prot_vaccinia.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- VACV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "VACV_PROTEINS"))

y$Virus <- "VACV"

df_tcell_I <- bind_rows(df_tcell_I, y)

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_MSA_PROTEINS)) {
  
  msa_path <- here("Data", "UGENE_Data", "Proteins", paste0(MPV_MSA_PROTEINS[i], "_nuc_MSA_transl.aln"))
  
  MSA <- seqinr::read.alignment(file = msa_path, format = "clustal", forceToLower = FALSE)
  
  prot <- MSA$seq[MSA$nam == "ON563414.2_Monkeypox_virus_isolate_MPXV_USA_2022_MA001,_complete_genome(translated)"]
  
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_MSA_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_MSA_PROTEINS"))

y$Virus <- "MPV2022"

df_tcell_I <- bind_rows(df_tcell_I, y)

df_tcell_I$Class <- "I"

##### T cell assay class II

epitopes <- tcell_II$Description...12 %>% stringr::str_squish() %>% unique()

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(MPV_PROTEINS[i], "_refprot.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_PROTEINS"))

y$Virus <- "MPV"

df_tcell_II <- y

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(VACV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(VACV_PROTEINS[i], "_prot_vaccinia.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- VACV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "VACV_PROTEINS"))

y$Virus <- "VACV"

df_tcell_II <- bind_rows(df_tcell_II, y)

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_MSA_PROTEINS)) {
  
  msa_path <- here("Data", "UGENE_Data", "Proteins", paste0(MPV_MSA_PROTEINS[i], "_nuc_MSA_transl.aln"))
  
  MSA <- seqinr::read.alignment(file = msa_path, format = "clustal", forceToLower = FALSE)
  
  prot <- MSA$seq[MSA$nam == "ON563414.2_Monkeypox_virus_isolate_MPXV_USA_2022_MA001,_complete_genome(translated)"]
  
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_MSA_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_MSA_PROTEINS"))

y$Virus <- "MPV2022"

df_tcell_II <- bind_rows(df_tcell_II, y)

df_tcell_II$Class <- "II"

z <- bind_rows(df_tcell_I %>% group_by(Epitope) %>% summarise(n=n()) %>% filter(n==1), 
               df_tcell_II %>% group_by(Epitope) %>% summarise(n=n()) %>% filter(n==1))


df_tcell <- bind_rows(df_tcell_I, df_tcell_II)

temp <- tcell_I %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) 
temp <- bind_rows(temp, tcell_II %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")))
colnames(temp) <- c("Epitope", "HLA")
df_tcell <- df_tcell %>% left_join(y = temp)

###### MHC assays data #####

##### MHC assay class I

epitopes <- mhc_I$Description...12 %>% stringr::str_squish() %>% unique()

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(MPV_PROTEINS[i], "_refprot.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_PROTEINS"))

y$Virus <- "MPV"

df_mhc_I <- y

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(VACV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(VACV_PROTEINS[i], "_prot_vaccinia.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- VACV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "VACV_PROTEINS"))

y$Virus <- "VACV"

df_mhc_I <- bind_rows(df_mhc_I, y)

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_MSA_PROTEINS)) {
  
  msa_path <- here("Data", "UGENE_Data", "Proteins", paste0(MPV_MSA_PROTEINS[i], "_nuc_MSA_transl.aln"))
  
  MSA <- seqinr::read.alignment(file = msa_path, format = "clustal", forceToLower = FALSE)
  
  prot <- MSA$seq[MSA$nam == "ON563414.2_Monkeypox_virus_isolate_MPXV_USA_2022_MA001,_complete_genome(translated)"]
  
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_MSA_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_MSA_PROTEINS"))

y$Virus <- "MPV2022"

df_mhc_I <- bind_rows(df_mhc_I, y)

df_mhc_I$Class <- "I"


##### MHC assay class II

epitopes <- mhc_II$Description...12 %>% stringr::str_squish() %>% unique()

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(MPV_PROTEINS[i], "_refprot.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_PROTEINS"))

y$Virus <- "MPV"

df_mhc_II <- y

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(VACV_PROTEINS)) {
  
  prot <- here("Data", "Refseqs", paste0(VACV_PROTEINS[i], "_prot_vaccinia.fasta")) %>% 
    seqinr::read.fasta() %>% unlist %>% as.character() %>% seqinr::c2s() %>% toupper()
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- VACV_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "VACV_PROTEINS"))

y$Virus <- "VACV"

df_mhc_II <- bind_rows(df_mhc_II, y)

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (i in seq_along(MPV_MSA_PROTEINS)) {
  
  msa_path <- here("Data", "UGENE_Data", "Proteins", paste0(MPV_MSA_PROTEINS[i], "_nuc_MSA_transl.aln"))
  
  MSA <- seqinr::read.alignment(file = msa_path, format = "clustal", forceToLower = FALSE)
  
  prot <- MSA$seq[MSA$nam == "ON563414.2_Monkeypox_virus_isolate_MPXV_USA_2022_MA001,_complete_genome(translated)"]
  
  
  out <- map.REFepitopes(Tepitopes = epitopes, REF = prot, NUM_MAX_MISMATCH = 0)
  
  out$Protein <- MPV_MSA_PROTEINS[i]
  
  y <- bind_rows(y, out)
  
}

y <- y %>% filter(!is.na(length))

y <- y %>% left_join(PROTEINS, keep = TRUE, by = c("Protein" = "MPV_MSA_PROTEINS"))

y$Virus <- "MPV2022"

df_mhc_II <- bind_rows(df_mhc_II, y)

df_mhc_II$Class <- "II"

z <- bind_rows(df_mhc_I %>% group_by(Epitope) %>% summarise(n=n()) %>% filter(n==1), 
               df_mhc_II %>% group_by(Epitope) %>% summarise(n=n()) %>% filter(n==1))


df_mhc <- bind_rows(df_mhc_I, df_mhc_II)

temp <- mhc_I %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) 
temp <- bind_rows(temp, mhc_II %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")))
colnames(temp) <- c("Epitope", "HLA")
df_mhc <- df_mhc %>% left_join(y = temp)


##### Misc analyses ####


HLA_mhc_I <- mhc_I %>% select(Description...12, `Allele Name`) %>% distinct() 

colnames(HLA_mhc_I) <- c("Epitope", "HLA")

HLA_mhc_I <- HLA_mhc_I %>%
  group_by(Epitope) %>% 
  summarise(HLA = paste(HLA, collapse = ",")) %>% 
  ungroup()

HLA_mhc_I <- HLA_mhc_I %>% mutate(HLA = str_remove_all(HLA, "HLA-"))

HLA_mhc_II <- mhc_II %>% select(Description...12, `Allele Name`) %>% distinct()

colnames(HLA_mhc_II) <- c("Epitope", "HLA")

HLA_mhc_II <- HLA_mhc_II %>%
  group_by(Epitope) %>% 
  summarise(HLA = paste(HLA, collapse = ",")) %>% 
  ungroup()

HLA_mhc_II <- HLA_mhc_II %>% mutate(HLA = str_remove_all(HLA, "HLA-"))


tcell_I %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) %>% dim()

tcell_II %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) %>% dim()

mhc_I %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) %>% dim()

mhc_II %>% distinct(across(c(Description...12, `Allele Name`))) %>% filter(str_detect(`Allele Name`, ":")) %>% dim()


df_tcell %>% 
  filter(Epitope %in% z$Epitope) %>% 
  distinct(across(Epitope), .keep_all = TRUE) %>% 
  group_by(MPV_PROTEINS) %>% 
  summarise(n=n())

df_tcell %>% 
  distinct(across(Epitope), .keep_all = TRUE) %>% 
  group_by(MPV_PROTEINS) %>% 
  summarise(n=n())


df_tcell %>% 
  distinct(across(c(Epitope, Virus)), .keep_all = TRUE) %>% 
  group_by(Virus) %>% 
  summarise(n=n())

df_tcell %>% 
  filter(Epitope %in% z$Epitope) %>% 
  distinct(across(Epitope), .keep_all = TRUE) %>% 
  group_by(MPV_MSA_PROTEINS) %>% 
  summarise(n=n()) 

df_mhc %>% 
  distinct(across(c(Epitope, Virus)), .keep_all = TRUE) %>% 
  group_by(Virus) %>% 
  summarise(n=n())

df_mhc %>% 
  filter(Epitope %in% z$Epitope) %>% 
  distinct(across(Epitope), .keep_all = TRUE) %>% 
  group_by(MPV_MSA_PROTEINS) %>% 
  summarise(n=n()) 


c <- df_tcell %>% 
  distinct(across(c(Epitope, Virus)), .keep_all = TRUE) %>% filter(Virus == "MPV2022", is.na(HLA)) %>% 
  pull(Epitope) %>% unique()

df_mhc %>% 
  filter(Epitope %in% c) %>% distinct(across(c(Epitope, HLA))) %>% 
  group_by(Epitope) %>% summarise(n=n()) %>% view()


