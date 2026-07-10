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

aa_out <- paste0("data/COGSeqs/",icog, ".fasta")
nuc_out <- paste0("data/COGSeqs/",icog, ".nuc.fasta")
file.create(nuc_out)
file.create(aa_out)
tmet <- c()
batch_size <- 10
for(seq_start in seq(0, length(prot_ids)-1, by = batch_size)) {
  print(paste0("Downloading seqs: ", seq_start+1, "-", seq_start+10))
  
  #select a batch
  iprots <- prot_ids[(seq_start+1):min(seq_start+batch_size, length(prot_ids))]
  
  #translate protein ID to nucleotide
  gene_ids <- entrez_fetch(db = "protein", 
                           id=iprots, rettype = 'ipg') 
  gid <- read.table(text=gene_ids, sep = '\t', header = T,fill=T) %>% 
    filter(Source=="RefSeq") %>% 
    group_by(Protein) %>% 
    mutate(order=row_number()) %>% 
    filter(order==1) %>% select(-order) #sekect ONE representative per protein seq
  tmet <- rbind(tmet, gid)
  
  #download seqs
  nuc_batch <- c()
  aa_batch <- c()
  for (i in 1:nrow(gid)){
  rep <- entrez_fetch(db = "nuccore", 
                              id = gid[i,'Nucleotide.Accession'],
                              seq_start = gid[i,'Start'],
                              seq_stop = gid[i,'Stop'],
                              strand=gid[i,'Strand'],  
                              rettype = "fasta", 
                              retmode = "text")
  nuc_batch <- c(nuc_batch, rep)
  
  repp <- entrez_fetch(db = "protein", 
                       id = gid[i,'Protein'],  
                       rettype = "fasta", 
                       retmode = "text")
  aa_batch <- c(aa_batch, repp)
  }
  
  #write files
  write(aa_batch, file = aa_out, append = TRUE)
  write(nuc_batch, file = nuc_out, append = TRUE)
  
  Sys.sleep(2) #politeness check
}
write_csv(tmet, "data/COGSeqs/COG3253.meta.csv")
