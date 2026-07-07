library(tidyverse)
library(httr2)
library(jsonlite)
library(rentrez)

setwd("~/WorkForaging/Teaching/iResearch/Summer2026/zaiden")

options(timeout = 300000)

#### Download COGs

#description of COGS
ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.def.tab" # 
local_file <- "data/COGDefs.tsv" #
download.file(url = ftp_url, 
              destfile = local_file, 
              method = "auto", 
              quiet = FALSE)
cog_desc<- read_tsv(local_file,quote = "", col_names = F)
colnames(cog_desc) <- c("COGID","ProtCat","Description","Symbol","Pathway","PubMedID","PDBID")

#master list of all proteins
ftp_url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/cog-24.cog.csv" # cog-24.def.tab
local_file <- "data/COGLocs.csv" #
download.file(url = ftp_url, 
              destfile = local_file, 
              method = "auto", 
              quiet = FALSE)
cog_data <- read_csv(local_file,quote = "", col_names = F)

#quick search for dismutases
cog_desc %>% filter(grepl("dismutase",Description))
icog <- "COG3253" #the chlorite dismutase is right here


#### Download sequences
prot_ids <- cog_data %>% filter(X7==icog) %>% pull(X3) %>% unique() 

fasta_out <- paste0("data/COGSeqs/",icog, ".fasta")
file.create(fasta_out)
batch_size <- 10
for(seq_start in seq(0, length(prot_ids)-1, by = batch_size)) {
  
  fasta_batch <- entrez_fetch(db = "protein", 
                              id = prot_ids[seq_start:min(seq_start+batch_size, length(prot_ids))],
                              rettype = "fasta", 
                              retmode = "text")
  
  
  write(fasta_batch, file = fasta_out, append = TRUE)
  
  Sys.sleep(2) #politeness check
}


