
# Loading required R packages

library(here)
library(tidyverse)
library(magrittr)
library(lubridate)
library(RColorBrewer)
library(ggrepel)

# Sourcing functions

source("~/Rprojects/monkeypox/Functions/compare.MPX.VACV.R")
source("~/Rprojects/monkeypox/Functions/prepare.protein.MSAs.R")
source("~/Rprojects/monkeypox/Functions/compare.VACV.MPV2022.genes.R")
source("~/Rprojects/monkeypox/Functions/compare.VACV.MPVCA.genes.R")
source("~/Rprojects/monkeypox/Functions/map.REFepitopes.R")
source("~/Rprojects/monkeypox/Functions/map.epitopes.VACV.R")
source("~/Rprojects/monkeypox/Functions/map.epitopes.MPVCA.R")
source("~/Rprojects/monkeypox/Functions/map.epitopes.MPV2022.R")

# Loading sequence data

VACV_gb <- genbankr::readGenBank(file = here("Data", "Refseqs", "NC_006998.gb"))
VACV_genes <- VACV_gb %>% genbankr::genes()
genes_IR <- (VACV_genes %>% ranges())
VACV_genome <- VACV_gb %>% genbankr::getSeq()

VACV_GENES <- na.omit(VACV_genes$gene)

meta_df <- read_tsv(file = here("Data", "2022", "gisaid_pox_2022_06_23.tsv"))

meta_df <- meta_df %>% 
  mutate(`Collection date` = lubridate::parse_date_time(x = `Collection date`, 
                                                        orders = c("ymd", "ym"), 
                                                        truncated = 3) %>% as_date()) %>%
  mutate(`Submission date` = lubridate::parse_date_time(x = `Submission date`, 
                                                        orders = c("ymd", "ym"), 
                                                        truncated = 3) %>% as_date())

msa_gisaid <- Biostrings::readDNAStringSet(filepath = here("Data", "2022", "MSA_full_24_06_2022.fasta"), 
                                           format = "fasta")

seqs_names <- msa_gisaid %>% names() 

seqs_ids <- meta_df %>% filter(`Collection date` >= ymd("2022-01-01")) %>% pull(`Accession ID`)

msa_gisaid[msa_gisaid %>% names() %in% seqs_ids]

msa_genomes <- Biostrings::DNAStringSet()

for (i in seqs_ids) {
  msa_genomes <- c(msa_genomes,
                   msa_gisaid[str_which(string = seqs_names, pattern = i)])
}


# Writing the MPXV-2022 full genome MSA and GISAID accession IDs

Biostrings::writeXStringSet(x = msa_genomes, 
                            filepath = here("Data", "2022", "MPX_2022_msa.fasta"))

GISAID_IDs <- names(msa_genomes) %>% str_split(pattern = "\\|") %>% unlist

GISAID_IDs <- GISAID_IDs[GISAID_IDs %>% str_which(pattern = "EPI_ISL_")]

write_csv(x = GISAID_IDs %>% as.data.frame(), col_names = FALSE, file = here("Data", "GISAID_IDs.csv"))

# Adding and aligning the full genomes of VACV and MPXV-CB to MPXV-2022 MSA

setwd(here("Data", "2022"))
system("mafft --thread -1 --add NC_006998_refseq_vaccinia.fasta MPX_2022_msa.fasta > out1.fasta")
system("mafft --thread -1 --add NC_003310refseq_monkeypox.fasta out1.fasta > out.fasta")
setwd(here())


# Preparing MSAs for all VACV proteins and orthologs in MPXV-2022 and MPXV-CB

for (GENE in VACV_GENES) {
  prepare.protein.MSAs(GENE = GENE)
}


# Comparing the translated protein sequences of MPXV-CB and MPXV-2022 with VACV 

df_protein_comparisons <- data.frame(Protein = vector(mode = "character"),
                                     num_subs = vector(mode = "numeric"),
                                     num_indels = vector(mode = "numeric"),
                                     width = vector(mode = "numeric"), 
                                     stringsAsFactors = FALSE)

df_protein_comparisons_MPV2022 <- data.frame(Protein = vector(mode = "character"),
                                             num_subs = vector(mode = "numeric"),
                                             num_indels = vector(mode = "numeric"),
                                             width = vector(mode = "numeric"), 
                                             stringsAsFactors = FALSE)


for (i in seq_along(VACV_GENES)) {
  
  print(VACV_GENES[i])
  
  df_protein_comparisons <- bind_rows(df_protein_comparisons, 
                                      data.frame(Protein = VACV_GENES[i],
                                                 MPVCA = compare.VACV.MPVCA.genes(PROTEIN = VACV_GENES[i]),
                                                 # MPVWA = compare.VACV.MPVWA.genes(GENE = VACV_GENES[i]),
                                                 stringsAsFactors = FALSE)) 
  
  df_protein_comparisons_MPV2022 <- bind_rows(df_protein_comparisons_MPV2022,
                                              compare.VACV.MPV2022.genes(GENE = VACV_GENES[i]))
  
  
}

saveRDS(object = list(df_protein_comparisons, df_protein_comparisons_MPV2022), 
        file = here("Data", "comparison-dataframes.RDS"))

temp1 <- df_protein_comparisons %>% 
  select(Protein, starts_with("MPVCA")) %>% 
  mutate(case = "MPVCA")

colnames(temp1) <- c("Protein", "num_subs", "num_indels", "width", "case")

test0 <- temp1 %>% 
  mutate(Dissimilarity = (num_subs+num_indels)/width*100) %>% 
  group_by(Protein, case) %>% 
  summarise(meanDis = mean(Dissimilarity),
            width = unique(width),
            medianDis = median(Dissimilarity)) %>% ungroup()

test1 <- df_protein_comparisons_MPV2022 %>% rowwise() %>% 
  mutate(Dissimilarity = (num_subs+num_indels)/width*100) %>% 
  group_by(Protein) %>% 
  summarise(meanDis = mean(Dissimilarity),
            width = unique(width),
            medianDis = median(Dissimilarity),
            p05Dis = quantile(Dissimilarity, probs = 0.05),
            p95Dis = quantile(Dissimilarity, probs = 0.95),
            p10Dis = quantile(Dissimilarity, probs = 0.1),
            p90Dis = quantile(Dissimilarity, probs = 0.9)) %>% ungroup() %>% 
  mutate(case = "MPV2022")  

test <- bind_rows(test0, test1) %>% 
  filter(Protein %in% VACV_IMM_PROTEINS) %>%
  distinct() %>% rowwise() %>% 
  mutate(prot_label = paste0(Protein, " (", width, ")")) %>% 
  ungroup()

PROTEIN_FACTOR <- test %>% filter(case == "MPV2022") %>% arrange(desc(medianDis)) %>% pull(Protein) 
names(PROTEIN_FACTOR) <- test %>% filter(case == "MPV2022") %>% arrange(desc(medianDis)) %>% pull(prot_label) 

# Plotting to check mean similarity of proteins

## boxplot

g2 <- test %>% filter(Protein %in% PROTEIN_FACTOR_t1) %>% 
  ggplot(mapping = aes(x = case %>% factor(levels = c("MPVCA", "MPV2022"), 
                                           labels = c("MPXV-CA", "MPXV-2022")), 
                       y = 100-meanDis, color = case))

g2 <- g2 + geom_boxplot(outlier.alpha = 0.75, show.legend = FALSE)

g2 <- g2 + theme_minimal() + theme(legend.position = "none", title = element_text(size = 11))

g2 <- g2 + 
  ggtitle("VACV immunogenic proteins") +
  ylab("Similarity (%)") + 
  xlab("")


## point plot

g1 <- test %>% filter(Protein %in% PROTEIN_FACTOR) %>% 
  ggplot(mapping = aes(x = Protein %>% factor(levels = PROTEIN_FACTOR, 
                                              labels = PROTEIN_FACTOR %>% names()), 
                       y = 100-medianDis, color = case))

g1 <- g1 + geom_errorbar(aes(ymax = 100-p90Dis, ymin = 100-p10Dis))
g1 <- g1 + geom_point(alpha = 0.55)
# g1 <- g1 + geom_text(aes(label = width, x = -0.1), check_overlap = TRUE)
g1 <- g1 + 
  scale_color_brewer(palette = "Set1", name = "Comparison with:") +
  scale_y_continuous(limits = c(10,100), breaks = seq.int(10,100, 10), expand = expansion(add = c(0.2,0.5)))
g1 <- g1 + theme_minimal() + theme(legend.position = "none", 
                                   axis.text.x = element_text(angle = 90, size = 8))
g1 <- g1 + 
  # ggtitle("VACV immunogenic proteins comparison with MPXV") + 
  ylab("Similarity (%)") + 
  xlab("Protein")

##### T cell epitope data analysis ####

# Mapping of T cell epitopes onto: (i) MPV2022, (ii) MPVCA refseq, (iii) MPVWA refseq

tcell_I <- read_csv("Data/IEDB/tcell_table_class_I/tcell_table_export_1653555280.csv", 
                    skip = 1)

tcell_II <- read_csv("Data/IEDB/tcell_table_class_II/tcell_table_export_1653555365.csv", 
                     skip = 1)

epitopes_I <- tcell_I$Description...12 %>% stringr::str_squish() %>% unique()

epitopes_II <- tcell_II$Description...12 %>% stringr::str_squish() %>% unique()

df_epitopes <- data.frame(Epitopes = c(epitopes_I, epitopes_II) %>% unique(), stringsAsFactors = FALSE)

df_epitopes <- df_epitopes %>% rowwise() %>% mutate(length = str_length(Epitopes))

y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (g in VACV_GENES) {
  
  out <- map.epitopes.VACV(GENE = g, epitopes = c(epitopes_I, epitopes_II))
  
  out$Protein <- g
  
  y <- bind_rows(y, out)
}

df_epi_map_VACV <- y

df_epi_map_VACV <- df_epi_map_VACV %>% filter(!is.na(start))

VACV_IMM_PROTEINS <- df_epi_map_VACV %>% pull(Protein) %>% unique()


y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (g in VACV_GENES) {
  
  out <- map.epitopes.MPV2022(PROTEIN = g, epitopes = c(epitopes_I, epitopes_II))
  
  out$Protein <- g
  
  y <- bind_rows(y, out)
}

df_epi_map_MPV2022 <- y

df_epi_map_MPV2022 <- df_epi_map_MPV2022 %>% filter(!is.na(start))


y <- data.frame(Epitope = vector(mode = "character"), 
                Protein = vector(mode = "character"), 
                start = vector(mode = "numeric"),
                end = vector(mode = "numeric"),
                length = vector(mode = "numeric"),
                stringsAsFactors = FALSE)

for (g in VACV_GENES) {
  
  out <- map.epitopes.MPVCA(PROTEIN = g, epitopes = c(epitopes_I, epitopes_II))
  
  out$Protein <- g
  
  y <- bind_rows(y, out)
}

df_epi_map_MPVCA <- y

df_epi_map_MPVCA <- df_epi_map_MPVCA %>% filter(!is.na(start))

saveRDS(object = list(df_epi_map_VACV, df_epi_map_MPVCA, df_epi_map_MPV2022), 
        file = here("Data", "epitopes-mappings.RDS"))


df_epi <- bind_rows(df_epi_map_VACV %>% select(Epitope, Protein, length) %>% distinct() %>% mutate(Virus = "VACV"),
                    df_epi_map_MPVCA %>% select(Epitope, Protein, length) %>% distinct() %>% mutate(Virus = "MPVCA"),
                    df_epi_map_MPV2022 %>% select(Epitope, Protein, length) %>% distinct() %>% mutate(Virus = "MPV2022"))

t <- df_epi %>% filter(Virus == "VACV") %>% filter(length <= 21) %>% pull(Epitope) %>% unique() 
t <- union(t, c("IDNESGWKTLVSRAIDLSSKK", "CIDGKWNPILPTCVR", "NSWNVIPSCQQKCDI"))