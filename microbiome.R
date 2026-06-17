#New analysis of microbiome data from first experiment

#Set working directory
setwd("U:/Adriana_Tables_For_Analysis/Sandra/Analysis/Sandra2/Run2/Phylo_M_filtered")

#Load necessary libraries
library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(vegan)        #Stats
library(rncl)         #graphics
library(ggpubr)       #graphics
library(rstatix)      #Need for alpha diversity bar plots



#1. Import data####
otu_mat<- read_excel("ASVs.xlsx")
tax_mat<- read_excel("TAX_K.xlsx")
samples_df <- read_excel("MET.xlsx") 

#2. Convert in phyloseq object########################
#Transform into matrixes otu and tax tables (sample table can be left as data frame)

samples_df <- as.data.frame(samples_df)
otu_mat <- as.data.frame(otu_mat)
tax_mat <- as.matrix(tax_mat)

samples_df1 <- samples_df[,-1]
rownames(samples_df1) <- samples_df[,1]

otu_mat1 <- otu_mat[,-1]
rownames(otu_mat1) <- otu_mat[,1]

tax_mat1 <- tax_mat[,-1]
rownames(tax_mat1) <- tax_mat[,1]

otu_mat <- data.matrix(otu_mat1)


#transform to phyloseq objects

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat1)
samples = sample_data(samples_df1)


#import tree

tree <- read_tree("tree.nwk", errorIfNULL=FALSE)


Aphids <- phyloseq(OTU, TAX, samples, tree)
Aphids


#3. Filter mitochondria and chloroplast reads#########


filtered16S <- Aphids %>%
  subset_taxa(!(Family %in% c("f_Mitochondria"))) %>%
  subset_taxa( !(Order %in% c("o_Chloroplast")))

filtered16S


#4. Alpha rarefaction for alpha diversity####

#Rarefy to minumin library size
ps_rarefied <- rarefy_even_depth(filtered16S, sample.size = min(sample_sums(filtered16S)), replace = FALSE, rngseed = 123)

print(ps_rarefied)

#check the depth
depth <- min(sample_sums(ps_rarefied)) #Minimum sample size
depth

depth2 <- max(sample_sums(ps_rarefied)) #Maximum sample size
depth2



# Transpose the OTU matrix so samples are rows and ASVs are columns
otu_matrix3 <- t((as(otu_table(ps_rarefied), "matrix"))) 


#Plot
rarecurve(otu_matrix3, step = 20, sample = min(colSums(otu_matrix3)), col = "#998ec3", cex = 0.6)


#Check the number of reads for all samples
sample_sums(ps_rarefied)

#Save plot

# Set up the PDF output
pdf("271124_rarefaction.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
rarecurve(otu_matrix3, step = 20, sample = min(colSums(otu_matrix3)), col = "#998ec3", cex = 0.6)

# Close the PDF device
dev.off()


#Plot before alpha rarefaction and store

# Transpose the OTU matrix so samples are rows and ASVs are columns
otu_matrix3 <- t((as(otu_table(Aphids), "matrix"))) 


# Set up the PDF output
pdf("271124_no_rarefaction.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
rarecurve(otu_matrix3, step = 20, sample = min(colSums(otu_matrix3)), col = "#998ec3", cex = 0.6)

# Close the PDF device
dev.off()

#5. Alpha diversity analysis#########

#Start with basic richness analysis using the phyloseq object and vegan


richness <- estimate_richness(ps_rarefied, split = TRUE)

richness

head(richness)

# Phyloseq function: "plot richness" generates dot plots showing the results from above

plot_richness(ps_rarefied)

#Export the richness file to adapt for analysis
#Organise the levels in Rep for the graph
#export the richness file and organise to make the ggplot instead of making individual ones
write.csv(richness, "271124_richness.csv") #I had already done this, but leaving it now as I will need to remember how to do this

#Upload the modified file
richness <- read.csv("181124_richness.csv")


#I will only keep Shannon and Observed for the big analysis
rich <- subset(richness, Index %in% c("Observed"))

#Factor
rich$Rep <- factor(rich$Rep, levels = c("Before", "Herb2", "NoH2", "Herb4", "NoH4", "Bulk"))


#use the package rstatix and ggpubr to add p.adj values to the plots
rich

#Test statistics
library(dplyr)
library(rstatix)

stat.test <- rich %>%
   t_test(Value ~ Rep) 
stat.test

#Only show significant
stat.test <- stat.test %>%
  filter(p.adj < 0.05)


stat.test <- stat.test %>% add_xy_position(x = "Rep") #Adjust position

stat.test$y.position <- 2000 #The bar was too high, so adjust manually



observed <- ggplot(rich, aes(x = Rep, y = Value)) +
  geom_boxplot(aes(fill = Insect)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_fill_manual(values = c("#61576A","#fc8d59", "#91bfdb")) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray", linewidth = 0.8) + # Add vertical dashed line
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray", linewidth = 0.8) + 
  ylab("Observed") + xlab("")+  
  theme(plot.title = element_blank()) +
  theme_pubr() +
  ggtitle("Observed")

#Remove the legend and create a new object 
obs <- observed  + theme(legend.position = "none") # Remove legend from this plot

#If annotations inside the graph are wanted

observed + annotate("text", x = c(2,3), y = -5, label = "Group A", size = 4, color = "black")  # Add text below the first box


#now shannon


#I will only keep Shannon and Observed for the big analysis
rich <- subset(richness, Index %in% c("Shannon"))

#Factor
rich$Rep <- factor(rich$Rep, levels = c("Before", "Herb2", "NoH2", "Herb4", "NoH4", "Bulk"))


#use the package rstatix and ggpubr to add p.adj values to the plots
rich

#Test statistics
library(dplyr)
library(rstatix)

stat.test <- rich %>%
  group_by(Index) %>%
  t_test(Value ~ Rep) 
stat.test

#Only show significant
stat.test <- stat.test %>%
  filter(p.adj < 0.05)


stat.test <- stat.test %>% add_xy_position(x = "Rep")

rich$Insect

shannon <- ggplot(rich, aes(x = Rep, y = Value)) +
  geom_boxplot(aes(fill = Insect)) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_fill_manual(values = c("#61576A","#fc8d59", "#91bfdb")) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "gray", linewidth = 0.8) + # Add vertical dashed line
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "gray", linewidth = 0.8) + 
  ylab("Shannon") + xlab("")+
  theme_pubr() +
  ggtitle("Shannon")

sha <- shannon + theme(legend.position = "right")


# Combine using patchwork
library(patchwork)
combined_plot <- obs + sha  # Side by side


# Show combined plot
combined_plot
#save the plots

# Set up the PDF output
pdf("091224_alphadiv.pdf", width = 12, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF

combined_plot

# Close the PDF device
dev.off()


#ANOVA OF ALPHA DIVERSITY

#Observed
data <- cbind(sample_data(ps_rarefied), richness) #This richness comes from the phyloseq object
data$Depth <- sample_sums(ps_rarefied)

kinetic.richness.anova <- aov(Observed ~ SampleType*Insect*Time, data)
summary(kinetic.richness.anova) ## Depth is very significant

library(broom)
anova_df <- tidy(kinetic.richness.anova)

#export the table
library(flextable)

# Create the flextable
ft <- flextable(anova_df) %>%
  set_caption("Anove of Observed Alpha Diversity") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  

# Save to Word
save_as_docx(ft, path = "anovaobserved_table.docx")


#Shannon

kinetic.richness.anova <- aov(Shannon ~ SampleType*Insect*Time, data)
summary(kinetic.richness.anova) ## Depth is very significant

library(broom)
anova_df <- tidy(kinetic.richness.anova)

#export the table
library(flextable)

# Create the flextable
ft <- flextable(anova_df) %>%
  set_caption("Anova of Shannon Alpha Diversity") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  

# Save to Word
save_as_docx(ft, path = "anovashannon_table.docx")






#6. Taxonomy bar plot################

library(mia)
library(miaViz)


#filter by presence 

filter_by_presence<-function(physeq,sampling_group,threshold){
  if(class(physeq)=="phyloseq"){
    if(!taxa_are_rows(physeq)){
      phyloseq::otu_table(physeq)<-t(phyloseq::otu_table(physeq))
    }
    sampling_group<-sample_data(physeq)[[sampling_group]]    
  }
  groups<-as.character(unique(sampling_group))
  cols_per_group<-sapply(groups,function(x) which(sampling_group==x),simplify=FALSE)
  ncols_per_group<-lengths(cols_per_group)
  
  # x is the row of OTU, y is the farm id
  
  if(class(physeq)=="phyloseq"){
    matrix<-t(apply(otu_table(physeq),1,function(x) sapply(cols_per_group,function(y) 100-(sum(x[y]==0)/length(y)*100 ))))
  } else{
    matrix<-t(apply(physeq,1,function(x) sapply(cols_per_group,function(y) 100-(sum(x[y]==0)/length(y)*100 ))))
  }
  n_groups_over_threshold       <- apply(matrix,1,function(x) length(which(x>threshold)))
  keep<-names(which(n_groups_over_threshold!=0))
  if(class(physeq)=="phyloseq"){
    physeq<-prune_taxa(keep,physeq)
  } else {
    physeq<-physeq[which(n_groups_over_threshold!=0),]
  }
  return(physeq)
}


#filter

filtered16S <- filter_by_presence(Aphids, "Rep", 50) #Note that this is not rarefied data
filtered16S

#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

#Agglomerate by Rep
tse_sub <- agglomerateByVariable(tse, by = "cols", group = "Rep")
tse_sub


# Computing relative abundance
tse_sub <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")

# Getting top taxa on a Phylum level
tse_sub <- agglomerateByRank(tse_sub, rank ="Phylum")
top_taxa <- getTop(tse_sub, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse_sub)$Phylum, function(x){
  if (x %in% top_taxa) {x} else {"Other"}
})
rowData(tse_sub)$Phylum_sub <- as.character(phylum_renamed)
# Agglomerate the data based on specified taxa
tse_sub2 <- agglomerateByVariable(tse_sub, by = "rows", f = "Phylum_sub")

# Visualizing the composition barplot, with samples order by "Bacteroidetes"

plotAbundance(
  tse_sub2, assay.type = "relabundance",
  order.row.by = "abund", add_x_text = TRUE)

#Save the abundance plot

# Set up the PDF output
pdf("021224_abundance_phylum.pdf", width = 5, height = 5)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF


plotAbundance(
  tse_sub2, assay.type = "relabundance",
  order.row.by = "abund", add_x_text = TRUE)


# Close the PDF device
dev.off()


#Make bar plot to the Class level

#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

#Agglomerate by Rep
tse_sub <- agglomerateByVariable(tse, by = "cols", group = "Rep")
tse_sub


# Computing relative abundance
tse_sub <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")


#Agglomerate by class
tse_sub <- agglomerateByRank(tse_sub, rank ="Class")
top_taxa <- getTop(tse_sub, top = 10, assay.type = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
class_renamed <- lapply(rowData(tse_sub)$Class, function(x){
  if (x %in% top_taxa) {x} else {"Other"}
})
rowData(tse_sub)$Class_sub <- as.character(class_renamed)
# Agglomerate the data based on specified taxa
tse_sub2 <- agglomerateByVariable(tse_sub, by = "rows", f = "Class_sub")



# Visualizing the composition barplot, with samples order by "Bacteroidetes"

plotAbundance(
  tse_sub2, assay.type = "relabundance",
  order.row.by = "abund", add_x_text = TRUE)


#Save the abundance plot

# Set up the PDF output
pdf("021224_abundance_class.pdf", width = 5, height = 5)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF


plotAbundance(
  tse_sub2, assay.type = "relabundance",
  order.row.by = "abund", add_x_text = TRUE)


# Close the PDF device
dev.off()


#Get the relabundance table

library(SummarizedExperiment)

# Extract the relative abundance matrix
relabundance_table <- assay(tse_sub2, "relabundance")

# Combine metadata with the relative abundance table
relabundance_table_with_metadata <- cbind(
  as.data.frame(colData(tse_sub2)),  # Metadata
  t(as.data.frame(relabundance_table))  # Transposed relative abundances
)

#save
write.csv(relabundance_table_with_metadata, "relative_abundance_table.csv", row.names = TRUE)

#7. Beta diversity#####################
remove.packages(c("mia", "miaViz"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(ask = FALSE)
BiocManager::install(c("mia", "miaViz"), ask = FALSE, update = TRUE)
#filter

filtered16S <- filter_by_presence(Aphids, "Rep", 50) #Note that this is not rarefied data
filtered16S

#Convert to mia object
tse <- mia::convertFromPhyloseq(filtered16S)



# Load package to plot reducedDim
library(scater)



# Beta diversity metrics like Bray-Curtis are often
# applied to relabundances
tse <- transformAssay(
  tse, assay.type = "counts", method = "relabundance")

# Other metrics like Aitchison to clr-transformed data
tse <- transformAssay(
  tse, assay.type = "relabundance", method = "clr", pseudocount = TRUE)


# Run PCoA on relabundance assay with Bray-Curtis distances
tse <- runMDS(
  tse,
  FUN = getDissimilarity,
  method = "bray",
  assay.type = "relabundance",
  name = "MDS_bray")

# Create ggplot object
p <- plotReducedDim(tse, "MDS_bray", colour_by = "Insect",  shape_by = "Time")

# Calculate explained variance
e <- attr(reducedDim(tse, "MDS_bray"), "eig")
rel_eig <- e / sum(e[e > 0])

# Add explained variance for each axis
p <- p +  labs(
  x = paste("PCoA 1 (", round(100 * rel_eig[[1]], 1), "%", ")", sep = ""),
  y = paste("PCoA 2 (", round(100 * rel_eig[[2]], 1), "%", ")", sep = "")
) +  scale_color_manual(values = c("#61576A","#fc8d59", "#91bfdb")) +
  scale_fill_manual(values = c("#61576A","#fc8d59", "#91bfdb")) 

p


#Store PCoA


#Save
# Set up the PDF output
pdf("021224_PCoA_all.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
p

# Close the PDF device
dev.off()

# Run NMDS on relabundance assay with Bray-Curtis distances
tse <- runNMDS(
  tse,
  FUN = getDissimilarity,
  method = "bray",
  assay.type = "relabundance",
  name = "NMDS_bray")


# Load package for multi-panel plotting
library(patchwork)


# Generate multi-panel plot
plots <- lapply(
  c("NMDS_bray"),
  function(dim) {
    p <- plotReducedDim(
      object = tse,
      dim = dim,
      colour_by = "Insect",
      shape_by = "Time"
    )
    # Add custom color scale
    p + scale_color_manual(values = c("#61576A","#fc8d59", "#91bfdb")) +
      scale_fill_manual(values = c("#61576A","#fc8d59", "#91bfdb"))
  }
)

# Combine plots if needed
combined_plot <- wrap_plots(plots)
print(combined_plot)


#Save the NMDS

#Save
# Set up the PDF output
pdf("021224_NMDS_all.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
print(combined_plot)

# Close the PDF device
dev.off()


#Run dbRDA

#For this, I need to remove the bulk soil as I want to see the effect of herbivory and
#sample type

#Remove bulk soil
# Subset by sample
tse_sub <- tse[ , tse$Time %in% c("Second", "Third")]



# Perform RDA
tse_sub <- addRDA(
  tse_sub,
  assay.type = "relabundance",
  formula = assay ~ Insect + Condition(Time),
  distance = "bray",
  na.action = na.exclude)

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse_sub, "RDA"), "significance")

#PERMANOVA table
rda_info$permanova |>
  knitr::kable()


# Load packages for plotting function
library(miaViz)

# Generate RDA plot coloured by clinical status
test <- plotRDA(tse_sub, "RDA", colour.by = "Insect")

test
test + scale_color_manual(values = c("Herbivory" = "#fc8d59", "NoHerb" = "#91bfdb"))


#Save
# Set up the PDF output
pdf("021224_dbRDA_Insect.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
test + scale_color_manual(values = c("Herbivory" = "#fc8d59", "NoHerb" = "#91bfdb"))

# Close the PDF device
dev.off()



#Check for homogeneity

rda_info$homogeneity |>
  knitr::kable()


#Make and export the tables

library(broom)


#export the table
library(flextable)


# Create the flextable
ft <- flextable(rda_info$permanova) %>%
  set_caption("PERMANOVA_dbRDA") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  

# Save to Word
save_as_docx(ft, path = "PERMANOVA_dbRDA_table.docx")


#Now, assuming Time as a covariate


# Perform RDA
tse_sub <- addRDA(
  tse_sub,
  assay.type = "relabundance",
  formula = assay ~ Insect + Condition(Time),
  distance = "bray",
  na.action = na.exclude)

# Store results of PERMANOVA test
rda_info <- attr(reducedDim(tse_sub, "RDA"), "significance")

#PERMANOVA table
rda_info$permanova |>
  knitr::kable()


#Make and export the tables

library(broom)


#export the table
library(flextable)


# Create the flextable
ft <- flextable(rda_info$permanova) %>%
  set_caption("PERMANOVA_dbRDA") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  

# Save to Word
save_as_docx(ft, path = "PERMANOVA_partialdbRDA_table.docx")



######Perform the PERMANOVA in phyloseq, as it's easier to test the interactions
#between SampleType and Insect

nosoil <- subset_samples(filtered16S, SampleType %in% c("Rhizosphere"))

#Calculate the relative abundance - total sum normalization (TSS)

relative_emilio <- transform_sample_counts(filtered16S, function(x) x/sum(x)*100)


pcoa_bray <- ordinate(relative_emilio, "PCoA", "bray")

plot_ordination (relative_emilio, pcoa_bray, type="samples", color="Insect", shape = "Time", title="PCoA") + geom_point(size=3) + scale_colour_manual(values = c("#61576A","#fc8d59", "#91bfdb")) + scale_shape_manual(values=c(15,3,16,4))


#PERMANOVA

#calculate bray curtis distance matrix
relative_emilio_bray <- phyloseq::distance(relative_emilio, method="bray")

#make a emilio frame from the relative abundance emilio
relative_emiliodf <- data.frame(sample_data(relative_emilio))

#adonis test (Test adonis on Treatment and Niche and interaction of both (*). If I did not want to test the interaction, choose +). 
adonis2(relative_emilio_bray ~ Insect*Time, data=relative_emiliodf, by= "terms")


#test homogeneity of variances (group dispersion) - because sometimes adonis result can be not because of variables but because of variance of data itself
betadispersion <- betadisper(relative_emilio_bray, relative_emiliodf$Insect)

# run a permutation test to get a statistic and a significance score
# if significant, (p < 0.01), we can reject the null hypothesis that our groups have the same dispersions
# thus, we have to accept that our emiliosets do not have the same variances. 
# So we cannot be so confident with our adonis results, it could be due to differences in group dispersions
# So we should be cautious in interpreting our Bray-Curtis ordination. 
# We can not simply say that the communities consist of different OTUs. 
# Among the groups there is a difference in the abundances between the OTUs. 
# Or better the evenness is different between these (two?) groups.

permutest(betadispersion)




#8. Venn diagram#####################

#data
filtered16S

bsoil <- subset_samples(filtered16S, Rep =="BS")
first <- subset_samples(filtered16S, Rep =="NH")
second_herb <- subset_samples(filtered16S, Rep =="RH")
second_noherb <- subset_samples(filtered16S, Rep =="RN")
third_herb <- subset_samples(filtered16S, Rep =="TH")
third_noherb <- subset_samples(filtered16S, Rep =="TN")


library(VennDiagram)

# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(first))
asv_table_treatment2 <- as.data.frame(otu_table(second_herb))
asv_table_treatment3 <- as.data.frame(otu_table(second_noherb))

asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))
asv_binary_treatment3 <- as.data.frame(ifelse(asv_table_treatment3 > 0, 1, 0))

# Get ASVs for each treatment
asv_list <- list(
  First = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  Herbivory = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0],
  NoHerb = rownames(asv_binary_treatment3)[rowSums(asv_binary_treatment3) > 0]
)

# Draw Venn Diagram for treatments
root <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c('#ffea95',"#fc8d59", "#91bfdb"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot the diagram
grid.draw(root)

dev.off()


dev.new()


#Only the second time


# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(second_herb))
asv_table_treatment2 <- as.data.frame(otu_table(second_noherb))

asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))


# Get ASVs for each treatment
asv_list <- list(
  Herbivory = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  NoHerb = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0]
)

# Draw Venn Diagram for treatments
second <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c("#fc8d59", "#91bfdb"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot the diagram
grid.draw(second)

dev.off()


dev.new()





#Do the third time


# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(third_herb))
asv_table_treatment2 <- as.data.frame(otu_table(third_noherb))

asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))


# Get ASVs for each treatment
asv_list <- list(
  Herbivory = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  NoHerb = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0]
)

# Draw Venn Diagram for treatments
third <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c("#fc8d59", "#91bfdb"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
) 

# Plot the diagram
grid.draw(third)

dev.off()


dev.new()


#Put together and save
grid.arrange(second, third, ncol =2)

#Save
# Set up the PDF output
pdf("021224_venn_secondthird.pdf", width = 8, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
grid.arrange(second, third, ncol =2)

# Close the PDF device
dev.off()


#Venn diagram only first and herbivory second


# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(first))
asv_table_treatment2 <- as.data.frame(otu_table(second_herb))


asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))

# Get ASVs for each treatment
asv_list <- list(
  First = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  Herbivory = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0]
)

# Draw Venn Diagram for treatments
comp <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c('#a7b039',"#fc8d59"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot the diagram
grid.draw(comp)

dev.off()


dev.new()




#Make the comparison of first and second no herbivory


# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(first))
asv_table_treatment2 <- as.data.frame(otu_table(second_noherb))


asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))

# Get ASVs for each treatment
asv_list <- list(
  First = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  NoHerb = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0]
)

# Draw Venn Diagram for treatments
comp2 <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c('#a7b039',"#91bfdb"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot the diagram
grid.draw(comp2)

dev.off()


dev.new()


#save the two latest graphs

grid.arrange(comp, comp2, ncol =2)



#Diagram for weeks 2 and 4 only herbivory treatments



#Make the comparison of first and second no herbivory


# Create a binary presence/absence for each treatment
asv_table_treatment1 <- as.data.frame(otu_table(second_herb))
asv_table_treatment2 <- as.data.frame(otu_table(third_herb))


asv_binary_treatment1 <- as.data.frame(ifelse(asv_table_treatment1 > 0, 1, 0))
asv_binary_treatment2 <- as.data.frame(ifelse(asv_table_treatment2 > 0, 1, 0))

# Get ASVs for each treatment
asv_list <- list(
  Second_Herb = rownames(asv_binary_treatment1)[rowSums(asv_binary_treatment1) > 0],
  Third_Herb = rownames(asv_binary_treatment2)[rowSums(asv_binary_treatment2) > 0]
)

# Draw Venn Diagram for treatments
comp2 <- venn.diagram(
  x = asv_list,
  filename = NULL,
  fill = c('#a7b039',"#91bfdb"),  # or any colors you prefer
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Plot the diagram

dev.off() #If need to stop 


dev.new() #If need a new space


grid.draw(comp2)


#Get the dataframe of shared and unique


df <- phyloseq_to_df(
  physeq,
  addtax = T,
  addtot = F,
  addmaxrank = F,
  sorting = "abundance"
)


#9. Differential abundance################

#subset for second time only
# Subset by sample
#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

tse_sub <- tse[ , tse$Time %in% c("Second")]

# Agglomerate by genus and subset by prevalence
tse <- subsetByPrevalent(tse_sub,rank = "Genus", prevalence = 10/100)
tse

# Transform count assay to relative abundances
tse <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")


# Load package
library(ANCOMBC)


#Need to rearrange the levels
tse$Insect <- factor(tse$Insect, levels = c("NoHerb", "Herbivory"))
levels(tse$Insect)

test2 = ancombc(data = tse, tax_level = "Genus", 
                formula = "Insect", 
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                group = "Insect", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
                max_iter = 100, conserve = TRUE, alpha = 0.10, global = TRUE,
                n_cl = 1, verbose = TRUE)


# Create the LFC DataFrame with taxa IDs and log-fold changes
df_lfc <- data.frame(test2$res$lfc[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Create the SE DataFrame with taxa IDs and standard errors
df_se <- data.frame(test2$res$se[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())


# Update column names to include "SE" for standard error columns
colnames(df_se)[-1] <- paste0(colnames(df_se)[-1], "SE")

# Merge LFC and SE data, focusing on the "Insect" column for analysis
df_fig_InsectHerbivory <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, InsectHerbivory = InsectHerbivory) %>%
  dplyr::filter(InsectHerbivory != 0) %>%
  dplyr::arrange(desc(InsectHerbivory)) %>%
  dplyr::mutate(direct = ifelse(InsectHerbivory > 0, "Positive LFC", "Negative LFC"))


#Add the errobars as they are not in the dataset
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  dplyr::left_join(df_se, by = "taxon_id")


# Factorize taxon_id and direction for plotting purposes
df_fig_InsectHerbivory$taxon_id <- factor(df_fig_InsectHerbivory$taxon_id, levels = df_fig_InsectHerbivory$taxon_id)
df_fig_InsectHerbivory$direct <- factor(df_fig_InsectHerbivory$direct, levels = c("Positive LFC", "Negative LFC"))



#Plot

second = ggplot(data = df_fig_InsectHerbivory, 
               aes(x = taxon_id, y = InsectHerbivory, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = InsectHerbivory - InsectHerbivorySE, ymax = InsectHerbivory + InsectHerbivorySE), width = 0.2,
                position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes after two weeks of herbivory") + 
  scale_fill_manual(values =c("#fc8d59", "#91bfdb")) +
  scale_color_manual(values =c("#fc8d59", "#91bfdb")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1))

second+ coord_flip()



#Save
# Set up the PDF output
pdf("021224_diffabund_2week.pdf", width = 8, height = 6)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
second+ coord_flip()

# Close the PDF device
dev.off()




#Third week

#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

tse_sub <- tse[ , tse$Time %in% c("Third")]

# Agglomerate by genus and subset by prevalence
tse <- subsetByPrevalent(tse_sub,rank = "Genus", prevalence = 10/100)
tse

# Transform count assay to relative abundances
tse <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")


# Load package
library(ANCOMBC)


#Need to rearrange the levels
tse$Insect <- factor(tse$Insect, levels = c("NoHerb", "Herbivory"))
levels(tse$Insect)

test2 = ancombc(data = tse, tax_level = "Genus", 
                formula = "Insect", 
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                group = "Insect", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
                max_iter = 100, conserve = TRUE, alpha = 0.10, global = TRUE,
                n_cl = 1, verbose = TRUE)


# Create the LFC DataFrame with taxa IDs and log-fold changes
df_lfc <- data.frame(test2$res$lfc[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Create the SE DataFrame with taxa IDs and standard errors
df_se <- data.frame(test2$res$se[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())


# Update column names to include "SE" for standard error columns
colnames(df_se)[-1] <- paste0(colnames(df_se)[-1], "SE")

# Merge LFC and SE data, focusing on the "Insect" column for analysis
df_fig_InsectHerbivory <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, InsectHerbivory = InsectHerbivory) %>%
  dplyr::filter(InsectHerbivory != 0) %>%
  dplyr::arrange(desc(InsectHerbivory)) %>%
  dplyr::mutate(direct = ifelse(InsectHerbivory > 0, "Positive LFC", "Negative LFC"))


#Add the errobars as they are not in the dataset
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  dplyr::left_join(df_se, by = "taxon_id")


# Factorize taxon_id and direction for plotting purposes
df_fig_InsectHerbivory$taxon_id <- factor(df_fig_InsectHerbivory$taxon_id, levels = df_fig_InsectHerbivory$taxon_id)
df_fig_InsectHerbivory$direct <- factor(df_fig_InsectHerbivory$direct, levels = c("Positive LFC", "Negative LFC"))



third = ggplot(data = df_fig_InsectHerbivory, 
                aes(x = taxon_id, y = InsectHerbivory, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = InsectHerbivory - InsectHerbivorySE, ymax = InsectHerbivory + InsectHerbivorySE), width = 0.2,
                position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes after four weeks of herbivory") + 
  scale_fill_manual(values =c("#fc8d59", "#91bfdb")) +
  scale_color_manual(values =c("#fc8d59", "#91bfdb")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1))

third+coord_flip()



#Save
# Set up the PDF output
pdf("021224_diffabund_4week.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
third+coord_flip()

# Close the PDF device
dev.off()


#10. Differential abundance To the ASV level######################


#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

tse_sub <- tse[ , tse$Time %in% c("Second")]


# Transform count assay to relative abundances
tse <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")


# Load package
library(ANCOMBC)


#Need to rearrange the levels
tse$Insect <- factor(tse$Insect, levels = c("NoHerb", "Herbivory"))
levels(tse$Insect)

test2 = ancombc(data = tse, 
                formula = "Insect", 
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                group = "Insect", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)


# Create the LFC DataFrame with taxa IDs and log-fold changes
df_lfc <- data.frame(test2$res$lfc[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Create the SE DataFrame with taxa IDs and standard errors
df_se <- data.frame(test2$res$se[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())


# Update column names to include "SE" for standard error columns
colnames(df_se)[-1] <- paste0(colnames(df_se)[-1], "SE")

# Merge LFC and SE data, focusing on the "Insect" column for analysis
df_fig_InsectHerbivory <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, InsectHerbivory = InsectHerbivory) %>%
  dplyr::filter(InsectHerbivory != 0) %>%
  dplyr::arrange(desc(InsectHerbivory)) %>%
  dplyr::mutate(direct = ifelse(InsectHerbivory > 0, "Positive LFC", "Negative LFC"))


##Try to add the taxonomy to the table, so I can make a dotplot
# Extract the taxonomy table from the TSE object
taxonomy <- rowData(tse)

#convert to dataframe
test<- as.data.frame(taxonomy)

test$taxon_id <- rownames(taxonomy)  # Add row names as a new column


#Add the errobars and taxonomy as they are not in the dataset
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::left_join(test, by = "taxon_id")
  


# Factorize taxon_id and direction for plotting purposes
df_fig_InsectHerbivory$taxon_id <- factor(df_fig_InsectHerbivory$taxon_id, levels = df_fig_InsectHerbivory$taxon_id)
df_fig_InsectHerbivory$direct <- factor(df_fig_InsectHerbivory$direct, levels = c("Positive LFC", "Negative LFC"))



second_ASV = ggplot(data = df_fig_InsectHerbivory, 
               aes(x = taxon_id, y = InsectHerbivory, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = InsectHerbivory - InsectHerbivorySE, ymax = InsectHerbivory + InsectHerbivorySE), width = 0.2,
                position = position_dodge(0.05), color = "black") +
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes after two weeks of herbivory at ASV level") + 
  scale_fill_manual(values =c("#fc8d59", "#91bfdb")) +
  scale_color_manual(values =c("#fc8d59", "#91bfdb")) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1))

second_ASV +coord_flip()



#Save
# Set up the PDF output
pdf("021224_diffabund_2week_ASV.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
third+coord_flip()

# Close the PDF device
dev.off()



#Dot plot


library(ggplot2)

# Create the dot plot
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c('#00429d', '#3e67ae', '#618fbf', '#85b7ce', '#b1dfdb', '#ffcab9', '#fd9291', '#e75d6f', '#c52a52', '#93003a'))  # Customize colors for significance


#Save

# Set up the PDF output
pdf("021224_diffabund_2week_ASV_dot.pdf", width = 7, height = 6)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
# Create the dot plot
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c('#00429d', '#3e67ae', '#618fbf', '#85b7ce', '#b1dfdb', '#ffcab9', '#fd9291', '#e75d6f', '#c52a52', '#93003a'))  # Customize colors for significance



# Close the PDF device
dev.off()



##################Try to reorder

# Assign genus-level behavior based on ASV-level LFCs
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  group_by(Genus) %>%
  mutate(
    genus_behavior = case_when(
      all(InsectHerbivory > 0) ~ "Increase",
      all(InsectHerbivory < 0) ~ "Decrease",
      TRUE ~ "Both"
    )
  ) %>%
  ungroup()

# Arrange by genus behavior and average LFC magnitude
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  group_by(Genus) %>%
  mutate(avg_lfc = mean(InsectHerbivory)) %>%
  ungroup() %>%
  arrange(factor(genus_behavior, levels = c("Decrease", "Increase", "Both")), desc(abs(avg_lfc)))

# Update Genus factor levels based on the new order
df_fig_InsectHerbivory$Genus <- factor(df_fig_InsectHerbivory$Genus, levels = unique(df_fig_InsectHerbivory$Genus))


#Filter NA and unknown
# Filter out unwanted or unclassified genera
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  filter(
    !is.na(Genus),                           # Remove NA genera
    !Genus %in% c("g_", "g_uncultured", "g_uncultured bacterium")  # Remove unwanted genera
  )


#export CSV file of differential abundance
write.csv(test2$res$q_val, "181224_qval.csv")
write.csv(test2$res$W, "181224_w.csv")
write.csv(test2$res$p_val, "181224_pval.csv")
write.csv(df_fig_InsectHerbivory, "181224_all_lfc.csv")
write.csv(test2$res$lfc, "181224_lfc.csv")

#extract the w
w <- test2$res$W
names(w)[names(w) == "taxon"] <- "taxon_id" #change column name

q <- test2$res$q_val
names(q)[names(q) == "taxon"] <- "taxon_id" #change column name

p <- test2$res$p_val
names(p)[names(p) == "taxon"] <- "taxon_id" #change column name

lfc <- test2$res$lfc
names(lfc)[names(lfc) == "taxon"] <- "taxon_id" #change column name


#filter all tables by the selected ASVs for differential abundance

df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  dplyr::left_join(w, by = "taxon_id") %>%
  dplyr::left_join(p, by = "taxon_id") %>%
  dplyr::left_join(q, by = "taxon_id")



#As this has a lot of information, change names and keep only w, p, q phylum and genus

names(df_fig_InsectHerbivory)[names(df_fig_InsectHerbivory) == "InsectHerbivory.y"] <- "W" #change column name
names(df_fig_InsectHerbivory)[names(df_fig_InsectHerbivory) == "InsectHerbivory.x.x"] <- "p-val" #change column name
names(df_fig_InsectHerbivory)[names(df_fig_InsectHerbivory) == "InsectHerbivory.y.y"] <- "q-val" #change column name
names(df_fig_InsectHerbivory)[names(df_fig_InsectHerbivory) == "InsectHerbivory.x"] <- "LFC" #change column name

#Create a data frame with the variables I want to keep
subset_df <- c("taxon_id", "Phylum", "Genus","LFC", "W", "p-val", "q-val")

#Subset
new_df <- df_fig_InsectHerbivory[subset_df]




#Export the table


library(flextable)

# Create the flextable
ft <- flextable(new_df) %>%
  set_caption("ANCOMBC_2week") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  



# Save to Word
save_as_docx(ft, path = "ancombc_secondweek_filtered.docx")



# Create the dot plot with ordered genera
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c('#ffcab9', '#93003a', '#e75d6f', '#618fbf', '#00429d','#b1dfdb'))  # Customize colors for significance




#Save
# Set up the PDF output
pdf("181224_diffabund_2week_ASV.pdf", width = 7, height = 6)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c('#ffcab9', '#93003a', '#e75d6f', '#618fbf', '#00429d','#b1dfdb'))  # Customize colors for significance


# Close the PDF device
dev.off()



####Differential abundance week 4 at the ASV level


#Convert to mia object
tse <- convertFromPhyloseq(filtered16S)

tse_sub <- tse[ , tse$Time %in% c("Third")]


# Transform count assay to relative abundances
tse <- transformAssay(tse_sub, assay.type = "counts", method = "relabundance")


# Load package
library(ANCOMBC)


#Need to rearrange the levels
tse$Insect <- factor(tse$Insect, levels = c("NoHerb", "Herbivory"))
levels(tse$Insect)

test2 = ancombc(data = tse, 
                formula = "Insect", 
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
                group = "Insect", struc_zero = FALSE, neg_lb = TRUE, tol = 1e-5, 
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
                n_cl = 1, verbose = TRUE)



# Create the LFC DataFrame with taxa IDs and log-fold changes
df_lfc <- data.frame(test2$res$lfc[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())

# Create the SE DataFrame with taxa IDs and standard errors
df_se <- data.frame(test2$res$se[, -1] * test2$res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = test2$res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())


# Update column names to include "SE" for standard error columns
colnames(df_se)[-1] <- paste0(colnames(df_se)[-1], "SE")

# Merge LFC and SE data, focusing on the "Insect" column for analysis
df_fig_InsectHerbivory <- df_lfc %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, InsectHerbivory = InsectHerbivory) %>%
  dplyr::filter(InsectHerbivory != 0) %>%
  dplyr::arrange(desc(InsectHerbivory)) %>%
  dplyr::mutate(direct = ifelse(InsectHerbivory > 0, "Positive LFC", "Negative LFC"))


##Try to add the taxonomy to the table, so I can make a dotplot
# Extract the taxonomy table from the TSE object
taxonomy <- rowData(tse)

#convert to dataframe
test<- as.data.frame(taxonomy)

test$taxon_id <- rownames(taxonomy)  # Add row names as a new column


#Add the errobars and taxonomy as they are not in the dataset
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::left_join(test, by = "taxon_id")

#This is for the dotplot

# Factorize taxon_id and direction for plotting purposes
df_fig_InsectHerbivory$taxon_id <- factor(df_fig_InsectHerbivory$taxon_id, levels = df_fig_InsectHerbivory$taxon_id)
df_fig_InsectHerbivory$direct <- factor(df_fig_InsectHerbivory$direct, levels = c("Positive LFC", "Negative LFC"))


# Assign genus-level behavior based on ASV-level LFCs
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  group_by(Genus) %>%
  mutate(
    genus_behavior = case_when(
      all(InsectHerbivory > 0) ~ "Increase",
      all(InsectHerbivory < 0) ~ "Decrease",
      TRUE ~ "Both"
    )
  ) %>%
  ungroup()

# Arrange by genus behavior and average LFC magnitude
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  group_by(Genus) %>%
  mutate(avg_lfc = mean(InsectHerbivory)) %>%
  ungroup() %>%
  arrange(factor(genus_behavior, levels = c("Decrease", "Increase", "Both")), desc(abs(avg_lfc)))

# Update Genus factor levels based on the new order
df_fig_InsectHerbivory$Genus <- factor(df_fig_InsectHerbivory$Genus, levels = unique(df_fig_InsectHerbivory$Genus))




#Filter NA and unknown
# Filter out unwanted or unclassified genera
df_fig_InsectHerbivory <- df_fig_InsectHerbivory %>%
  filter(
    !is.na(Genus),                           # Remove NA genera
    !Genus %in% c("g_", "g_uncultured", "g_uncultured bacterium")  # Remove unwanted genera
  )


# Create the dot plot with ordered genera
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c( '#93003a',  '#e75d6f', '#618fbf', '#ffcab9', '#b1dfdb', '#00429d'))  # Customize colors for significance




#Save
# Set up the PDF output
pdf("031224_diffabund_4week_ASV.pdf", width = 8, height = 6)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
ggplot(df_fig_InsectHerbivory, aes(x = InsectHerbivory, y = Genus)) +
  geom_point(aes(color = Phylum), size = 3, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +  # Dashed line
  labs(
    x = "Log Fold Change",
    y = "Genus",
    title = "Differential Abundance Analysis by Genus and ASV"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = c( '#93003a',  '#e75d6f', '#618fbf', '#ffcab9', '#b1dfdb', '#00429d'))  # Customize colors for significance

# Close the PDF device
dev.off()



 write.csv(test2$res$diff_abn, "181224_fourweeks_diff.csv")


#11. Reads counts#################
#Convert to mia object
tse <- convertFromPhyloseq(Aphids)

library(scater)
# Get an overview of sample and taxa counts
summary(tse, assay.type= "counts")

#filter mitochondria and chloroplasts

filtered16S <- Aphids %>%
  subset_taxa(!(Family %in% c("f_Mitochondria"))) %>%
  subset_taxa( !(Order %in% c("o_Chloroplast")))

filtered16S

#Convert to mia object
tse_2 <- convertFromPhyloseq(filtered16S)

# Get an overview of sample and taxa counts
summary(tse_2, assay.type= "counts")