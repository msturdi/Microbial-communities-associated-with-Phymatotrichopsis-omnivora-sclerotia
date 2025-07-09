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

display.all.moma(6, colorblind_only=T)

##
#filtering out single otus in shared file
##
#shared <- read_tsv("final.opti_mcc.0.03.subsample.oct.24.bacterial.shared",
#col_types = cols(Group = col_character(),
#.default = col_double()))

#meta_cols <- c("label", "Group", "numOtus") #define columns to keep from shared file
#otu_sums <- colSums(shared[, !(names(shared) %in% meta_cols)]) #sum otu abundances across all samples
#otus_to_keep <- names(otu_sums[otu_sums > 1]) #keep otus with total counts above 1

#shared_filtered <- shared %>%
#select(all_of(meta_cols), all_of(otus_to_keep))

#otu_count <- shared_filtered %>%
#select(Group, starts_with("Otu")) %>%
#rename(sample = Group) %>%
#column_to_rownames(var = "sample")

otu_count <- read_tsv("final.opti_mcc.0.03.subsample.shared") %>%
  select(Group, starts_with("Otu")) %>%
  rename(sample = Group) %>%
  column_to_rownames(var = "sample")

otu_count_t <- as.data.frame(t(otu_count))

otumat <- as.matrix(otu_count_t) #use for phyloseq

taxonomy <- read_tsv(file="final.opti_mcc.0.03.cons.taxonomy") %>%
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

metadata <- read_excel("16S.2023.metadata.xlsx") %>%
  column_to_rownames(var = "sample")

sampledata = sample_data(data.frame(metadata))

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

physeq = phyloseq(OTU, TAX, sampledata)
physeq

tree <- read_tree("final.opti_mcc.0.03.rep.phylip.oct.23.bact.tre") #representative tree make using get.oturep in mothur

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

t1$cal_ses_betamntd(runs = 999, abundance.weighted = TRUE, null.model = "sample.pool") #runs should = 999, trying 99 to see how long it takes
t1$res_ses_betamntd

#View and store beta NTI results
betaNTI_results <- t1$res_ses_betamntd
print(betaNTI_results)

#Add beta NTI matrix to the beta_diversity list in the data set
dt2$beta_diversity[["betaNTI"]] <- betaNTI_results
print(names(dt2$beta_diversity))

t2 <- trans_beta$new(dataset = dt2, group = "trt", measure = "betaNTI")

#Calculate group distances
t2$cal_group_distance()

t2$res_group_distance$Process #NULL means done

#for combo
t2$sample_table$trt <- factor(t2$sample_table$trt, 
                                levels = c("live", "dead", "none"))

#for treatment
t2$sample_table$treatment <- factor(t2$sample_table$treatment, 
                                    levels = c("live", "dead", "none"))

##
#make CSV file so that you don't have to run the 999 step again
##
g1 <- t2$plot_group_distance(boxplot_add = "mean") #in g1 figure, above -1 is stochastic and below is deterministic

plot_data <- g1$data

write.csv(plot_data, file = "group.distance.oct.23.bact.csv", row.names = FALSE)
##

g2 = g1 +
  scale_color_manual(name=NULL,
                     breaks=c("live", "dead", "none"),
                     labels=c("Live sclerotia", "Heat-killed sclerotia", "Bulk soil"),
                     values = c(moma.colors("Sidhu", 3, direction=1, type="discrete"))) +
  scale_x_discrete(limits=c("live", "dead", "none"),
                   labels=c("Live sclerotia", "Heat-killed sclerotia", "Bulk soil")) +
  labs(x=NULL,
       y="betaNTI",
       title="Stochastic or Deterministic Selection<br>In-season 2023 Stiles Farm **Bacterial**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size = 14, color = "black"),
        axis.text.y = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 14),
        legend.position = "none",
        plot.title = element_markdown(hjust=0.5, size = 16))

ggsave("betaNTI.oct23.16S.png", width=12, height=8)

##
#from Vanessa's code, but doesn't seem to take into account the -2 threshold
##
#t2$res_group_distance$Process <- ifelse(
#abs(t2$res_group_distance$Value) > 2.0, "Deterministic", "Stochastic")

t2$res_group_distance$Process <- ifelse(t2$res_group_distance$Value > 2.0 | t2$res_group_distance$Value < -2.0, 
                                        "Deterministic", "Stochastic")

proportions <- t2$res_group_distance %>%
  group_by(trt, Process) %>%
  summarize(Count = n()) %>%
  mutate(Proportion = Count*100 / sum(Count)) %>%
  arrange(trt, Process)

proportions$trt <- factor(proportions$trt, 
                          levels = c("live", "dead", "none"))

proportions %>%
  ggplot(aes(x = trt, y = Proportion, fill = Process)) +
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
       title="Proportion of Stochastic vs. Deterministic Selection<br>In-season 2024 **Fungal**") +
  theme_classic() +
  theme(axis.text.x = element_markdown(size = 14, color = "black"),
        axis.text.y = element_markdown(size = 10),
        axis.title.y = element_markdown(size = 14),
        legend.position = "bottom",
        legend.text = element_markdown(size = 10),
        legend.title = element_text(size=11),
        legend.key.size = unit(10, "pt"),
        plot.title = element_markdown(hjust=0.5, size = 16))

ggsave("porportion.selection.oct24.ITS.png", width=12, height=8)