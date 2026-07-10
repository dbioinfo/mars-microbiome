library(tidyverse)
library(Biostrings)
setwd("~/WorkForaging/Teaching/iResearch/Summer2026/zaiden")


#read in specific COG
icog <- "COG3253" #the chlorite dismutase is right here

aa_out <- readAAStringSet(paste0("data/COGSeqs/",icog, ".fasta"))
nuc_out <- readDNAStringSet(paste0("data/COGSeqs/",icog, ".nuc.fasta"))

#remove duplicates
#aa_out <- aa_out[-c(435,145,79,80),]
#nuc_out <- nuc_out[-c(435,145,79,80),]

meta <- read_csv("data/COGSeqs/COG3253.meta.csv")[-c(435),] %>% 
  mutate(aa_header = paste(Protein, Protein.Name, Organism, sep=" "),
         nuc_header = paste0(Nucleotide.Accession,':',Start,'-',Stop),
         final_header = gsub(" ","_",paste(Organism, Protein.Name, Protein, sep="..")) )

nuc_names <- sub(" .*", "",names(nuc_out))
aa_names <- sub(" .*", "",names(aa_out))

new_nucs <- meta[match(nuc_names, meta$nuc_header),]$final_header
new_aa <- meta[match(aa_names, meta$Protein),]$final_header

#sanity check, do the fastas have the same order?
all(new_nucs==new_aa)

#if so, proceed and rewrite
names(aa_out) <- new_aa
names(nuc_out) <- new_nucs


#actually, need to translate the reverse compliment of the neg stranded dna
negstrands <- meta %>% filter(Strand=='-') %>% pull(final_header)
newnuc_out <- nuc_out
newnuc_out[negstrands,]  <- reverseComplement(nuc_out[negstrands,])

#remove more dups
newnuc_out <- newnuc_out[-c(145,79,80),]
newaa_out <- aa_out[-c(145,79,80),]
meta <- meta[-c(145,79,80),]

#remove all excessively long or short versions
meta <- meta %>% filter(Stop-Start>=650, Stop-Start<=800)

writeXStringSet(newaa_out[meta$final_header,], paste0("data/COGSeqs/",icog, ".trim.fasta"))
writeXStringSet(newnuc_out[meta$final_header,], paste0("data/COGSeqs/",icog, ".nuc.comp.fasta"))
write.csv(meta, "data/COGSeqs/COG3253.meta.final.csv")
tmp <- meta %>% mutate(group=case_when(grepl("chlor", final_header)~"ChloriteDismutase", grepl("heme", final_header)~"HemeQ", .default="HemeQ"))
write_csv(tmp, 'data/COGSeqs/COG3253.groups.txt')
