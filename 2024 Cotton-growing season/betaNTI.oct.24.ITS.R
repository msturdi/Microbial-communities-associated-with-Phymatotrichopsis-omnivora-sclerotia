library(readxl)
library(phyloseq)
library(magrittr)
library(picante)
library(microeco)
library(tidyverse)
library(glue)
library(file2meco)
library(ggtext)
library(dplyr)

set.seed(19980609)

otu_count <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.fungal.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  column_to_rownames(var = "sample")

otu_count_t <- as.data.frame(t(otu_count))

otumat <- as.matrix(otu_count_t) #use for phyloseq

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.oct.24.fungal.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="k__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="p__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="c__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="o__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="f__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="g__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="s__", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus"), 
           sep=";") %>%
  mutate(kingdom = str_replace_all(kingdom,
                             "(.*)_unclassified", "Unclassified \\1"),
         kingdom = str_replace_all(kingdom,
                             "unclassified_(.*)", "Unclassified \\1"),
         kingdom = str_replace_all(kingdom, 
                                 "_", " ")) %>%
  mutate(phylum = str_replace_all(phylum,
                                   "(.*)_unclassified", "Unclassified \\1"),
         phylum = str_replace_all(phylum,
                                   "unclassified_(.*)", "Unclassified \\1"),
         phylum = str_replace_all(phylum, 
                                   "_", " ")) %>%
  mutate(class = str_replace_all(class,
                                  "(.*)_unclassified", "Unclassified \\1"),
         class = str_replace_all(class,
                                  "unclassified_(.*)", "Unclassified \\1"),
         class = str_replace_all(class, 
                                  "_", " ")) %>%
  mutate(order = str_replace_all(order,
                                 "(.*)_unclassified", "Unclassified \\1"),
         order = str_replace_all(order,
                                 "unclassified_(.*)", "Unclassified \\1"),
         order = str_replace_all(order, 
                                 "_", " ")) %>%
  mutate(family = str_replace_all(family,
                                 "(.*)_unclassified", "Unclassified \\1"),
         family = str_replace_all(family,
                                 "unclassified_(.*)", "Unclassified \\1"),
         family = str_replace_all(family, 
                                 "_", " ")) %>%
  mutate(genus = str_replace_all(genus,
                                 "(.*)_unclassified", "Unclassified \\1"),
         genus = str_replace_all(genus,
                                 "unclassified_(.*)", "Unclassified \\1"),
         genus = str_replace_all(genus, 
                                 "_", " ")) %>%
  column_to_rownames(var = "otu")

taxmat <- as.matrix(taxonomy) #use for phyloseq

metadata <- read_excel("oct.24.ITS.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}")) %>%
  column_to_rownames(var = "sample")

sampledata = sample_data(data.frame(metadata))

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX, sampledata)
physeq

tree <- read_tree("final.opti_mcc.0.03.rep.phylip.oct.24.tre") #representative tree make using get.oturep in mothur

tree_taxa <- tree$tip.label #get current tree tip labels
otu_taxa <- taxa_names(OTU) #get OTU names from the table
tree$tip.label <- otu_taxa #rename tree tips to match OTU table

physeq = phyloseq(OTU, TAX, sampledata, tree)
physeq

meco_dataset <- phyloseq2meco(physeq)

meco_dataset$filter_taxa(rel_abund = 0.00001, freq = 0.1)
meco_dataset

dt <- trans_norm$new(meco_dataset)
dt2 <- dt$norm(method = "CSS")

t1 <- trans_nullmodel$new(dt2, filter_thres = 0.000005)
t1$data_tree

t1$cal_ses_betamntd(runs = 999, abundance.weighted = TRUE, null.model = "sample.pool")
t1$res_ses_betamntd

#View and store beta NTI results
betaNTI_results <- t1$res_ses_betamntd
print(betaNTI_results)

#Add beta NTI matrix to the beta_diversity list in the data set
dt2$beta_diversity[["betaNTI"]] <- betaNTI_results
print(names(dt2$beta_diversity))

t2 <- trans_beta$new(dataset = dt2, group = "combo", measure = "betaNTI")

#Calculate group distances
t2$cal_group_distance()

t2$res_group_distance$Process #NULL means done

#for combo
t2$sample_table$combo <- factor(t2$sample_table$combo, 
                                   levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))

g1 <- t2$plot_group_distance(boxplot_add = "mean")

plot_data <- g1$data

write.csv(plot_data, file = "group.distance.oct.24.ITS.csv", row.names = FALSE)
##

beta_NTI_ave <- read.csv("group.distance.oct.24.ITS.csv") %>%
  rename_all(tolower) %>%
  group_by(combo) %>%
  summarize(mean = mean(value),
            se = sd(value) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")
