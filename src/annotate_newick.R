library(tidyverse)
library(ape)
setwd("~/WorkForaging/Teaching/iResearch/Summer2026/zaiden")

#here we need to annotate the newick with the foreground labels
#we also need to rewrite the names of the taxa to be <30 chars
#then use fastcodeml and hope it works quickly

atr <- readLines("data/COGSeqs/COG3253.nucaln.paml")
gtr <- read.tree("data/COGSeqs/COG3253.nwk")
gtext <- readLines("data/COGSeqs/COG3253.nwk")
meta <- read_csv("data/COGSeqs/COG3253.groups.txt")
foreground_tips <- meta$final_header[meta$group == "ChloriteDismutase"]
meta <- dplyr::bind_rows(meta, data.frame(final_header=setdiff(gtr$tip.label, tips), group="HemeQ"))
meta <- meta %>% mutate(SeqID=paste0("ID",row_number()))
tips <- meta$final_header %>% unique()


for (tip in tips){
    rtip <- meta %>% filter(final_header==tip) %>% pull(SeqID)
    exact_match_pattern <- paste0("\\b", tip, "\\b")
    gtext <- gsub(pattern = exact_match_pattern, 
                  replacement = rtip, 
                  x = gtext)
    atr <- gsub(pattern = exact_match_pattern, 
                replacement = rtip, 
                x = atr)
    
    if (tip %in% foreground_tips) {
      exact_match_pattern <- paste0("\\b", rtip, "\\b")
      paml_tag <- paste0(rtip, "#1")
      
      gtext <- gsub(pattern = exact_match_pattern, 
                        replacement = paml_tag, 
                        x = gtext)
    }
}

#internal nodes labels
mrcanode <- getMRCA(gtr, foreground_tips)
clade_nodes <- phytools::getDescendants(gtr, mrcanode)
for (inode in clade_nodes){
  if(inode>Ntip(gtr)){
    node_idx <- inode - Ntip(gtr)
    gtr$node.label[node_idx] <- paste0(inode, "#1")
    }
}

test <- tidytree::read.tree(text=gtext)
test$tip.label[which(!grepl("ID", test$tip.label))] <- c("ID339","ID340","ID341",
                                                         "ID342","ID343")
ggtree(test, layout='circular')+geom_tiplab()


writeLines(gtext, "data/COGSeqs/COG3253.cld_labelled.nwk")
write_csv(meta, "data/COGSeqs/COG3253.meta.labelled.csv")
writeLines(atr, "data/COGSeqs/COG3253.nucaln.labelled.paml")



#internal nodes labels
gtr <- read.tree("data/COGSeqs/COG3253.cld_labelled.nwk")
gtr$node.label <- character(gtr$Nnode)
foreground_nodes <- paste0(meta$SeqID[meta$group == "ChloriteDismutase"],'#1')
mrcanode <- getMRCA(gtr, foreground_nodes)
clade_nodes <- phytools::getDescendants(gtr, mrcanode)
for (inode in clade_nodes){
  if(inode>Ntip(gtr)){
    node_idx <- inode - Ntip(gtr)
    gtr$node.label[node_idx] <- paste0(inode, "#1")
  }
}
write.tree(gtr, "data/COGSeqs/COG3253.cld_labelled.int.nwk")
