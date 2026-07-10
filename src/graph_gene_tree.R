library(tidyverse)
library(ggtree)
library(ape)
library(tidytree)
setwd("~/WorkForaging/Teaching/iResearch/Summer2026/zaiden")

gtr <- read.tree("data/codeml/COG3253/COG3253.fastcodeml.branches.nwk")
meta <- read_csv("data/COGSeqs/COG3253.meta.labelled.csv")
tmp <- meta %>% select(SeqID, final_header, group, Organism) %>% 
  mutate(SeqID=case_when(group=="ChloriteDismutase"~paste0(SeqID,"#1"),
                         .default=SeqID))

p<-ggtree(gtr, layout="radial") %<+% tmp  +
  geom_tiplab(aes(label=Organism),align = TRUE, size=1, linesize = 0.2)+
  geom_tippoint(aes(color=group)) + 
  scale_color_manual(values=c("#7ACFB1","#D85071"), name="Gene Function")
ggtree::rotate(p, 340) #fix groupings as needed
ggsave("figs/HemeQTreeCodeML.png")
