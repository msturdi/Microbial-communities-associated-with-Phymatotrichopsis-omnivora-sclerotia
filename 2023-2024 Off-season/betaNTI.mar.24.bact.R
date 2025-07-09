library(readxl)
library(phyloseq)
library(magrittr)
library(picante)
library(microeco)
library(tidyverse)
library(glue)
library(file2meco)
library(dplyr)
library(metagenomeSeq)
library(iCAMP)
library(MoMAColors)
library(ggtext)

set.seed(19980609)

##
#filtering out single OTUs in shared file
##
#shared <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.bact.shared",
#col_types = cols(Group = col_character(),
#.default = col_double()))

#meta_cols <- c("label", "Group", "numOtus") #define columns to keep from shared file
#otu_sums <- colSums(shared[, !(names(shared) %in% meta_cols)]) #sum OTU abundances across all samples
#otus_to_keep <- names(otu_sums[otu_sums > 1]) #keep OTUs with total counts above 1

#shared_filtered <- shared %>%
#select(all_of(meta_cols), all_of(otus_to_keep))

otu_count <- read_tsv("final.opti_mcc.0.03.subsample.mar.24.bact.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  column_to_rownames(var = "sample")

otu_count_t <- as.data.frame(t(otu_count))

otumat <- as.matrix(otu_count_t) #use for phyloseq

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.mar.24.bact.taxonomy") %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern="\\(\\d*\\)", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern='["]', replacement="")) %>%
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

metadata <- read_excel("mar.24.bact.metadata.xlsx") %>%
  mutate(combo = glue("{location} {treatment}")) %>%
  column_to_rownames(var = "sample")

sampledata = sample_data(data.frame(metadata))

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX, sampledata)
physeq

tree <- read_tree("final.opti_mcc.0.03.rep.phylip.mar.24.bact.tre")

tree_taxa <- tree$tip.label #get current tree tip labels
otu_taxa <- taxa_names(OTU) #get OTU names from the table
tree$tip.label <- otu_taxa #rename tree tips to match OTU table

physeq = phyloseq(OTU, TAX, sampledata, tree)
physeq

meco_dataset <- phyloseq2meco(physeq)

##
#Code from Vanessa to filter out very low abundance OTUs
##
meco_dataset$filter_taxa(rel_abund = 0.00001, freq = 0.1)
meco_dataset

#allow you to calculate beta NTI and other phylogenetic metrics. Use a threshold to filter low-abundance taxa if needed.
#beta:method = "CSS" (Cumulative Sum Scaling)
dt <- trans_norm$new(meco_dataset)
dt2 <- dt$norm(method = "CSS")

t1 <- trans_nullmodel$new(dt2, filter_thres = 0.000005)
t1$data_tree

t1$cal_ses_betamntd(runs = 999, abundance.weighted = TRUE, null.model = "sample.pool")
t1$res_ses_betamntd

#view and store beta NTI results
betaNTI_results <- t1$res_ses_betamntd
print(betaNTI_results)

#add beta NTI matrix to the beta_diversity list in the data set
dt2$beta_diversity[["betaNTI"]] <- betaNTI_results
print(names(dt2$beta_diversity))

t2 <- trans_beta$new(dataset = dt2, group = "combo", measure = "betaNTI")

#calculate group distances
t2$cal_group_distance()

t2$res_group_distance$Process #NULL means done

#for combo
t2$sample_table$combo <- factor(t2$sample_table$combo, 
                                levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))

#for treatment
#t2$sample_table$treatment <- factor(t2$sample_table$treatment, 
                                    #levels = c("live", "dead", "none"))

##
#make CSV file so that you don't have to run the 999 step again
##
g1 <- t2$plot_group_distance(boxplot_add = "mean")

plot_data <- g1$data

write.csv(plot_data, file = "group.distance.mar.24.csv", row.names = FALSE)

beta_NTI_ave <- read.csv("group.distance.mar.24.csv") %>%
  rename_all(tolower) %>%
  group_by(combo) %>%
  summarize(mean = mean(value),
            se = sd(value) / sqrt(n()),
            max = mean+se,
            min = mean-se,
            N = n(),
            .groups = "drop")

g2 = g1 +
  scale_color_manual(name=NULL,
                     breaks=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                     labels=c("BF live<br>sclerotia", "BF heat-killed<br>sclerotia", "BF bulk", 
                              "SF live<br>sclerotia", "SF heat-killed<br>sclerotia", "SF bulk"),
                     values = c(moma.colors("Sidhu", 6, direction=1, type="discrete"))) +
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("BF live<br>sclerotia", "BF heat-killed<br>sclerotia", "BF bulk", 
                            "SF live<br>sclerotia", "SF heat-killed<br>sclerotia", "SF bulk")) +
  labs(x=NULL,
       y="betaNTI",
       title="Stochastic or Deterministic Selection<br>Off-season 2023-24 **Bacterial**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size = 14, color = "black"),
        axis.text.y = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 14),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size = 16))

ggsave("betaNTI.mar24.16S.png", width=12, height=8)

##
#proportion analysis
##
t2$res_group_distance$Process <- ifelse(t2$res_group_distance$Value > 2.0 | t2$res_group_distance$Value < -2.0, 
                                        "Deterministic", "Stochastic")

proportions <- t2$res_group_distance %>%
  group_by(combo, Process) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count*100 / sum(Count)) %>%
  arrange(combo, Process)

proportions$combo <- factor(proportions$combo, 
                            levels = c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"))

proportions %>%
  ggplot(aes(x = combo, y = Proportion, fill = Process)) +
  geom_col() +
  scale_fill_manual(name="Selection Type",
                    breaks=c("Deterministic", "Stochastic"),
                    values=c('cadetblue3','coral2'),
                    labels=c("Deterministic", "Stochastic")) +
  scale_x_discrete(limits=c("BF live", "BF dead", "BF none", "Stiles live", "Stiles dead", "Stiles none"),
                   labels=c("BF live<br>sclerotia", "BF heat-killed<br>sclerotia", "BF bulk", 
                            "SF live<br>sclerotia", "SF heat-killed<br>sclerotia", "SF bulk")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL,
       y="Stochastic vs. Deterministic Selection",
       title="Proportion of Stochastic vs. Deterministic Selection<br>Off-season 2023-24 **Fungal**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size = 14, color = "black"),
        axis.text.y = element_markdown(size = 10),
        axis.title.y = element_markdown(size = 14),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10),
        legend.title = element_text(size=11),
        legend.key.size = unit(10, "pt"),
        plot.title = element_markdown(hjust=0.5, size = 16))

ggsave("porportion.selection.mar24.ITS.png", width=12, height=8)