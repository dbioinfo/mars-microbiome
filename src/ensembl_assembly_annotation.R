library(tidyverse)
library(httr2)
library(jsonlite)
library(rentrez)

setwd("~/WorkForaging/Teaching/iResearch/Summer2026/zaiden")

options(timeout = 300000)

#use the rest API url 
url <- "https://rest.ensembl.org/info/genomes/division/EnsemblBacteria"

#building the query  (return all bacterial genomes)
req <- request(url) %>% 
  req_headers("Content-Type" = "application/json") %>% 
  req_retry(max_tries = 3)

#sends the query 
res <- req_perform(req)

#unpack the query
dat <- res %>% resp_body_string() %>% fromJSON(flatten=T)

#search for a microbe
search_pats <- c() #'Azotobacter vinelandii', 'Chroococcidiopsis'
search_pat <- 'Azotobacter vinelandii' #Chroococcidiopsis
dat %>% filter(grepl(search_pat, scientific_name)) %>% select(scientific_name, assembly_name, base_count, species_taxonomy_id)

search_pats <- c(search_pats, search_pat)


#search all patterns 
final_search <- paste(search_pats, collapse='|')

search_results <- dat %>% 
  filter(grepl(final_search, scientific_name)) %>% 
  select(scientific_name, assembly_name, assembly_accession, base_count, species_taxonomy_id)


#write_csv(search_results, file = "data/ensembl_assemblies.csv")
blacklist<-c(573,72407,384,1423)
search_results <- search_results %>% filter(!(species_taxonomy_id %in% blacklist))


#quick graph of genome sizes
ggplot(search_results)+
geom_bar(aes(x=as.character(species_taxonomy_id), y=base_count, fill=assembly_accession),stat='identity',position = 'dodge')+
  theme_bw()+
  scale_y_continuous(expand=c(0,0))+
  guides(fill='none')+
  theme(axis.text.x=element_text(angle=90))+
  xlab("Species Tax ID")+
  ylab("Bases in Genome")+
  ggtitle("Search Results by Base Count")

#################
####NOTES
#################


#add in the environmental stressors that each organism resists / is adapted for
#reduce search window and produce pruned dataset
#create graph for poster of genome sizes in search results
