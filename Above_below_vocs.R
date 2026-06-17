
########################Final analysis aboveground volatiles#####################
#Date: 07-12-2024
#Created by: Sandra Cortes
setwd("C:/Users/cortess/OneDrive - Rothamsted Research/Rothamsted/First experiment/GC/Data for GC Aligner")

#Load libraries
library(GCalignR)
library(ggplot2)
library(gridExtra)
library(plot3D)
library(vegan)
library(dplyr)
library(ptw)
library(readxl)

#1.Load files##############

#Aphid herbivory second week

AG9 <- read.delim("AGA9.txt")
AG9

AG9A <- as.data.frame(AG9)
###

AG5 <- read.delim("AG5.txt")
AG5

AG5A <- as.data.frame(AG5)

###
AG12 <- read.delim("A12.txt")
AG12

AG12C <- as.data.frame(AG12)

###

AG14 <- read.delim("A14.txt")
AG14

AG14A <- as.data.frame(AG14)

###

AG2 <- read.delim("A2.txt")
AG2

AG2A <- as.data.frame(AG2)


###
AG6 <- read.delim("A6.txt")
AG6

AG6A <- as.data.frame(AG6)

###

AG15 <- read.delim("C15.txt")
AG15

AG15C <- as.data.frame(AG15)

###
AG10 <- read.delim("C10.txt")
AG10

AG10C <- as.data.frame(AG10)

###

AG4 <- read.delim("C4.txt")
AG4

AG4C <- as.data.frame(AG4)


#### Blank
Blank <- read_xlsx("Blank.xlsx")
Blank

Blank <- as.data.frame(Blank)


Blank2 <- read_xlsx("Blank.xlsx")
Blank2

Blank2 <- as.data.frame(Blank2)



#Load my data


LA7 <- read.delim("LA7.txt")
LA7

LA7 <- as.data.frame(LA7)
###

LA8 <- read.delim("LA8.txt")
LA8

LA8 <- as.data.frame(LA8)

###
LA12 <- read.delim("LA12.txt")
LA12

LA12 <- as.data.frame(LA12)


###

LA13 <- read.delim("LA13.txt")
LA13

LA13 <- as.data.frame(LA13)


###
LC6 <- read.delim("LC6.txt")
LC6

LC6 <- as.data.frame(LC6)


###
LC7 <- read.delim("LC7.txt")
LC7

LC7 <- as.data.frame(LC7)


###

LC8 <- read.delim("LC8.txt")
LC8

LC8 <- as.data.frame(LC8)

###
LC13 <- read.delim("LC13.txt")
LC13

LC13 <- as.data.frame(LC13)


###

LC14 <- read.delim("LC14.txt")
LC14

LC14 <- as.data.frame(LC14)


#2. Check that data passes GCalignR##########


aphid_all <- list(AG9A, AG5A, AG14A, AG2A, AG4C, AG10C, AG15C, AG12C, Blank, LA7, LA8, LA12, LA13, LC6, LC7, LC8, LC13, LC14)
aphid_all

names(aphid_all) <- c("AG9A", "AG5A", "AG14A", "AG2A", "AG4C", "AG10C", "AG15C", "AG12C", "Blank", "LA7", "LA8", "LA12", "LA13", "LC6", "LC7", "LC8", "LC13", "LC14")




check_input(aphid_all, plot =T)


#3. Check peak interspace to determine the best linear shift###############


peak_interspace(data = aphid_all, rt_col_name = "RT", quantile_range = c(0, 0.8), quantiles = 0.05)


#4. Align#######

peak_data_aligned <- align_chromatograms(data = aphid_all, # input data
                                         rt_col_name = "RT", # retention time variable name 
                                         rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                         rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                         reference = "AG10C", # choose automatically 
                                         max_linear_shift = 0.05, # max. shift for linear corrections
                                         max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                         min_diff_peak2peak = 0.03, # min. expected distance between peaks
                                         delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                         blanks = "Blank",
                                         write_output = NULL) # add variable names to write aligned data to text files


#Check the peaks that were aligned

gc_heatmap(peak_data_aligned)

plot(peak_data_aligned, which_plot = "all")

print(peak_data_aligned)



#5. Access the files from the alignment#################
#The aligned data matrices are now stored in data frames which can be accessed as follows:

peak_data_aligned$aligned$RT # to access the aligned retention times

peak_data_aligned$aligned$Area # to access the aligned areas

alignment_data <- peak_data_aligned$aligned$Area # to access the aligned area data

#Save the file

write.csv(alignment_data, "071224_aligneddata.csv", row.names = FALSE)




#6. Normalise the data ###############
norm_apcon <- norm_peaks(peak_data_aligned, conc_col_name = "Area", rt_col_name = "RT", out = "data.frame")
norm_apcon <- log(norm_apcon + 1)
norm_apcon

#Pareto scaling

pareto_scale <- function(df) {
  scaled_df <- scale(df, center = TRUE, scale = sqrt(apply(df, 2, sd, na.rm = TRUE)))
  return(as.data.frame(scaled_df))
}

# Assuming `norm_apcon` is your numeric data and `metadata` contains group information
numeric_data <- norm_apcon  # Replace with your numeric dataframe

# Pareto scale the data
scaled_data <- pareto_scale(numeric_data)

scaled_data

#7. Visualise using multivariate analysis###########################

#Upload metadata

metadata <- read.csv("metadata.csv")
##Keep order of rows consistent###

norm_apcon <- norm_apcon[match(row.names(metadata), row.names(norm_apcon)), ] 
norm_apcon


## NMDS using Bray-Curtis dissimilarities
norm_nmds <- vegan::metaMDS(comm = norm_apcon, distance = "bray")
norm_nmds


##I'm getting a warning, so I'll explore the stress plot
vegan::stressplot(norm_nmds)


## get x and y coordinates
norm_nmds <- as.data.frame(norm_nmds[["points"]])  
## add the colony as a factor to each sample
norm_nmds <- cbind(norm_nmds,Treatment = metadata[["Combined"]])
norm_nmds

## ordiplot with ggplot2
##This is not working, I need to use a different analysis
library(ggplot2)
ggplot(data = norm_nmds,aes(MDS1,MDS2)) +
  geom_point(aes(colour = Treatment)) + 
  theme_void() +
  scale_color_manual(values = c('#00429d', '#73a2c6', '#ffc0ad', '#ff7c03')) +
  theme_bw()
 

##I'll try now the multivariate analysis

adonis_result <- vegan::adonis2(
  norm_apcon ~ Treatment + Sampling + Treatment:Sampling,
  data = metadata,
  permutations = 999
)

adonis_result


# Load necessary libraries
library(vegan)
library(pairwiseAdonis)

# Example data: your distance matrix and group factor
distance_matrix <- vegdist(dist_matrix)  # replace with your distance matrix
group_factor <- factor(metadata$Combined)  # replace with your grouping variable

# Run PERMANOVA on the whole dataset
adonis_result <- adonis(distance_matrix ~ group_factor, data = metadata)

# Perform pairwise comparisons using pairwise.adonis
pairwise_result <- pairwise.adonis(distance_matrix, group_factor)

# Adjust p-values for multiple comparisons
pairwise_result$p.adjusted <- p.adjust(pairwise_result$p.value, method = "fdr")

# View the results with adjusted p-values
print(pairwise_result)

#export the table
library(flextable)

# Create the flextable
ft <- flextable(pairwise_result) %>%
  set_caption("Pairwise result with FDR") %>%  # Add a title
  colformat_double(digits = 3) %>%                  # Format numeric columns
  autofit()  

# Save to Word
save_as_docx(ft, path = "071224_pairwisepermanova_aboveVOCs_table.docx")




#####PCA with factomineR

library(FactoMineR)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")
library(factoextra)

res.pca <- PCA(norm_apcon,  graph = FALSE)
res.pca
get_eig(res.pca)

fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var

head(var$contrib)
fviz_pca_var(res.pca, col.var = "black")

fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)

fviz_contrib(res.pca, choice = "var", axes = 1, top = 20)

fviz_pca_biplot(res.pca, repel = TRUE)

fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

fviz_pca_ind(res.pca,
             habillage = metadata$Treatment, # color by groups
             addEllipses = TRUE # Concentration ellipses
)

library("corrplot")
corrplot(var$contrib, is.corr=FALSE) 

fviz_pca_ind(res.pca, pointsize = "cos2", 
             pointshape = 21, fill = "#E7B800",
             repel = TRUE # Avoid text overlapping (slow if many points)
)


pca_all <- fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = metadata$Combined, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups"
)

pca_all


pca_plot <- pca_all +
  aes(shape = metadata$Sampling, color = metadata$Treatment, fill = metadata$Treatment) +# Assign shapes based on another factor
  scale_shape_manual(values = c(17, 16, 17, 16, 16, 16)) + # Customize shape types
  scale_color_manual(values = c("#fc8d59", "#fc8d59", "#91bfdb", "#91bfdb")) + #Customise colours
  scale_fill_manual(values = c("#fc8d59", "#fc8d59", "#91bfdb", "#91bfdb")) +
  theme_minimal() + # Apply minimal theme
  theme(
    panel.grid = element_blank(), # Remove dashed lines
    axis.line = element_line(color = "black") # Add axis lines
  )


pca_plot


#8. Save plot############################

#Save
# Set up the PDF output
pdf("071224_PCA_abovegroundVOCs.pdf", width = 6, height = 4)  # Adjust size as needed

# Arrange the four Venn diagrams in a grid and save to PDF
pca_plot

# Close the PDF device
dev.off()



###############################################################################
#Heatmap compounds selected second week

#Load library for heatmap visualization

library(ComplexHeatmap)

#1. Load the file########

df_2 <- read.csv("100424_selectedpeaks.csv")

rownames(df_2) <-  df_2[,1]

#Remove the first 2 columns and convert to matrix

df_2 <- df_2[,-1] %>%
  as.matrix()

#Need to scale using zcore
scaled_matrix <- scale(df_2)


#Transpose
transposed_matrix <- t(scaled_matrix)



# Example annotation data frame
annotations <- data.frame(
  Metabolite = c("Unknown_739", "E-2Hexen-1-ol", "1,3 Dimethylbenzene", "5-Methyl-3-methylene-5-hexen-2-one", 
                 "Hexyl acetate", "Benzyl alcohol", "Limonene", "Z-Ocimene", "Acetophenone", "Isophorone", 
                 "α-Methylenephenylacetaldehyde", "Dodecene", "Myrtenol", "Z-3-Hexen-1-yl-3-methylbutyrate", 
                 "Indole", "α-cubebene", "3-Methylindole", "(Z)-Jasmone", "Citronellyl propionate", 
                 "Methyl jasmonate", "Unknown_1631", "Caryophyllene oxide", "Unknown_1867"),
  Pathway = c("Unknown", "LOX", "Shikimate", "LOX", "LOX", "LOX", "MEV/Non-MEV", "MEV/Non-MEV", 
              "Shikimate", "LOX", "Shikimate", "LOX", "MEV/Non-MEV", "LOX", "Shikimate", "MEV/Non-MEV", 
              "Shikimate", "LOX", "MEV/Non-MEV", "LOX", "Unknown", "MEV/Non-MEV", "Unknown"),
  Function = c("Unknown", "Signalling", "No", "Limited", "Signalling", "Antimicrobial", "Direct defence", 
               "Indirect defence", "Limited", "Limited", "Limited", "Limited", "Direct defence", "Signalling", 
               "Multifunctional", "Indirect defence", "Multifunctional", "Multifunctional", "Limited", "Signalling", 
               "Unknown", "Direct defence", "Unknown")
)



#Try a different heatmap color that is colorblind safe
library(circlize)

col_fun = colorRamp2(c(-1, 0, 1), c('#002a7f', '#ffffe0', '#e9002c'))

#Row annotation

class_colors <- c("Unknown" = '#99bc84', "LOX" = '#008ea8', "Shikimate"= '#ee6a66', "MEV/Non-MEV" ='#a54fc5')


 
#Check unique number of functions
unique(annotations$Function)
length(unique(annotations$Function))

function_colors <- c("Signalling" = '#5671d4', "No" = '#8a92a2', "Indirect defence" =  '#93003a', "Direct defence"= '#78dac9', "Multifunctional"= '#ffa77b', "Antimicrobial" ='#ee6a66', "Limited" = '#ca2f50', "Unknown"='#8db6a7')

#Heatmap object for row annotation
ha = rowAnnotation(
  Condition = annotations$Pathway,
  Function = annotations$Function,
  col = list(Condition = class_colors, 
             Function = function_colors)
  )

#Column name

metadata <- data.frame(
  Sample = c("Herb1", "Herb2", "Herb3", "Herb4", "NoHerb1", "NoHerb2", "NoHerb3", "NoHerb4"),
  Condition = c("Herbivory","Herbivory", "Herbivory","Herbivory", "NoHerb","NoHerb","NoHerb","NoHerb")
)

cond_colors <- c("Herbivory" = "#fc8d59", "NoHerb" ="#91bfdb")

ha_2 = HeatmapAnnotation(
  Condition = metadata$Condition,
  col = list(Condition = cond_colors
             ))




Heatmap(transposed_matrix, 
        name = "Metabolites",
        top_annotation = ha_2,
        left_annotation = ha,
        show_row_dend = FALSE,  # Hide the row dendrogram
        cluster_rows = FALSE,   # Disable clustering
        row_split = annotations$Pathway,  # Split rows by class
        column_split = metadata$Condition,
        show_row_names = TRUE,
        show_column_names = FALSE)

#test grayscale palette
# Define a grayscale palette
grayscale_palette <- colorRamp2(c(-2, 0, 2), c("black", "gray", "white"))

# Define a grayscale palette
red_palette <- colorRamp2(c(-2, 0, 2), c('white', '#ff7692','#870a65'))
green_palette <- colorRamp2(c(-2, 0, 2), c('white', '#6cbc81','#007000'))
otherpalette <- colorRamp2(c(-2, 0, 2), c('white', '#FFE1D4','#F53D2A'))



 
# Save heatmap as a file
pdf("heatmap_output_red_new.pdf", width = 8, height = 7)

#create image

Heatmap(transposed_matrix, 
        name = "Metabolites",
        top_annotation = ha_2,
        left_annotation = ha,
        col = otherpalette,
        show_row_dend = FALSE,  # Hide the row dendrogram
        cluster_rows = FALSE,   # Disable clustering
        row_split = annotations$Pathway,  # Split rows by class
        column_split = metadata$Condition,
        show_row_names = TRUE,
        show_column_names = FALSE)

#finish

dev.off()

##################ALKANES EXAMPLE#################################################

library(GCalignR)
library(ggplot2)
library(gridExtra)
library(plot3D)
library(vegan)
library(ptw)


#Compare two alkanes chromatograms to see if they align


alk1 <- read.delim("alkanes010722.txt")
alk2 <- read.delim("alkanes200522.txt")


alk1df <- as.data.frame(alk1)
alk2df <- as.data.frame(alk2)

alklist <- list(alk1df, alk2df)

names(alklist) <- c("alk1", "alk2")
alklist

check_input(alklist, plot =T)

peak_interspace(data = alklist, rt_col_name = "RT", quantile_range = c(0, 0.8), quantiles = 0.05)

peak_alklist_aligned <- align_chromatograms(data = alklist, # input data
                                            rt_col_name = "RT", # retention time variable name 
                                            rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                            rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                            reference = NULL, # choose automatically 
                                            max_linear_shift = 0.05, # max. shift for linear corrections
                                            max_diff_peak2mean = 0.3, # max. distance of a peak to the mean across samples
                                            min_diff_peak2peak = 0.3, # min. expected distance between peaks
                                            delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                            write_output = NULL) # add variable names to write aligned data to text files
#The aligned data matrices are now stored in data frames which can be accessed as follows:


gc_heatmap(peak_alklist_aligned)

plot(peak_alklist_aligned, which_plot = "all")

print(peak_alklist_aligned)
peak_alklist_aligned

norm_alklist <- norm_peaks(peak_alklist_aligned, conc_col_name = "Area", rt_col_name = "RT", out = "data.frame")
norm_alklist <- log(norm_alklist + 1)
norm_alklist

norm_alklist <- vegan::metaMDS(comm = norm_apcon, distance = "bray")
norm_nmds

######################BELOWGROUND##########################################
#First, upload the data




ba14 <- read.delim("Below_14.txt")
b14a <- as.data.frame(ba14)

ba2 <- read.delim("NBelow_2A.txt")
b2a <- as.data.frame(ba2)

ba5 <- read.delim("Below_A5.txt")
b5a <- as.data.frame(ba5)

ba6 <- read.delim("Below_A6.txt")
b6a <- as.data.frame(ba6)

ba9 <- read.delim("Below_A9.txt")
b9a <- as.data.frame(ba9)

belap <- list(b14a, b2a, b5a, b6a, b9a)
belap

names(belap) <- c("b14a", "b2a", "b5a", "b6a", "b9a")
belap

check_input(belap, plot =T)


peak_interspace(data = belap, rt_col_name = "RT", quantile_range = c(0, 0.8), quantiles = 0.05)

peak_belap_aligned <- align_chromatograms(data = belap, # input data
                                          rt_col_name = "RT", # retention time variable name 
                                          rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                          rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                          reference = NULL, # choose automatically 
                                          max_linear_shift = 0.05, # max. shift for linear corrections
                                          max_diff_peak2mean = 0.3, # max. distance of a peak to the mean across samples
                                          min_diff_peak2peak = 0.3, # min. expected distance between peaks
                                          delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                          write_output = NULL) # add variable names to write aligned data to text files
#The aligned data matrices are now stored in data frames which can be accessed as follows:


gc_heatmap(peak_belap_aligned)

plot(peak_belap_aligned, which_plot = "all")

print(peak_belap_aligned)




##Now, I'll try to change the parameters

peak_belap_aligned <- align_chromatograms(data = belap, # input data
                                          rt_col_name = "RT", # retention time variable name 
                                          rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                          rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                          reference = NULL, # choose automatically 
                                          max_linear_shift = 0.05, # max. shift for linear corrections
                                          max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                          min_diff_peak2peak = 0.03, # min. expected distance between peaks
                                          delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                          write_output = NULL) # add va





gc_heatmap(peak_belap_aligned)

plot(peak_belap_aligned, which_plot = "all")

print(peak_belap_aligned)

##I will leave the max_diff_peak2mean in 0.03 


#Normalise the data


norm_belap <- norm_peaks(peak_belap_aligned, conc_col_name = "Area", rt_col_name = "RT", out = "data.frame")
norm_belap <- log(norm_belap + 1)
norm_belap

##Include the metadata table
metabelow <- belowmeta
metabelow



norm_belap <- norm_belap[match(row.names(metabelow), row.names(norm_belap)), ] 
norm_belap


###############Including the control


ba14 <- read.delim("Below_14.txt")
b14a <- as.data.frame(ba14)

ba2 <- read.delim("NBelow_2A.txt")
b2a <- as.data.frame(ba2)

ba5 <- read.delim("Below_A5.txt")
b5a <- as.data.frame(ba5)

ba6 <- read.delim("Below_A6.txt")
b6a <- as.data.frame(ba6)

ba9 <- read.delim("Below_A9.txt")
b9a <- as.data.frame(ba9)


bc3 <- read.delim("Below_C3.txt")
b3c <- as.data.frame(bc3)

bc10 <- read.delim("Below_C10.txt")
b10c <- as.data.frame(bc10)


bc12 <- read.delim("Below_C12.txt")
b12c <- as.data.frame(bc12)

bc14 <- read.delim("Below_C14.txt")
b14c <- as.data.frame(bc14)

bc15 <- read.delim("Below_C15.txt")
b15c <- as.data.frame(bc15)



metadata <- Bmetadata



belac <- list(b14a, b2a, b5a, b6a, b9a, b3c, b10c, b12c, b14c, b15c)
belac

names(belac) <- c("b14a", "b2a", "b5a", "b6a", "b9a", "b3c", "b10c", "b12c", "b14c", "b15c")
belac

check_input(belac, plot =T)


peak_interspace(data = belac, rt_col_name = "RT", quantile_range = c(0, 0.8), quantiles = 0.05)

peak_belac_aligned <- align_chromatograms(data = belac, # input data
                                          rt_col_name = "RT", # retention time variable name 
                                          rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                          rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                          reference = NULL, # choose automatically 
                                          max_linear_shift = 0.05, # max. shift for linear corrections
                                          max_diff_peak2mean = 0.3, # max. distance of a peak to the mean across samples
                                          min_diff_peak2peak = 0.3, # min. expected distance between peaks
                                          delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                          write_output = NULL) # add variable names to write aligned data to text files
#The aligned data matrices are now stored in data frames which can be accessed as follows:


gc_heatmap(peak_belac_aligned)

plot(peak_belap_aligned, which_plot = "all")

print(peak_belap_aligned)




##Now, I'll try to change the parameters

peak_belac_aligned <- align_chromatograms(data = belac, # input data
                                          rt_col_name = "RT", # retention time variable name 
                                          rt_cutoff_low = 5, # remove peaks below 15 Minutes
                                          rt_cutoff_high = 40, # remove peaks exceeding 45 Minutes
                                          reference = NULL, # choose automatically 
                                          max_linear_shift = 0.05, # max. shift for linear corrections
                                          max_diff_peak2mean = 0.03, # max. distance of a peak to the mean across samples
                                          min_diff_peak2peak = 0.03, # min. expected distance between peaks
                                          delete_single_peak = TRUE, # delete peaks that are present in just one sample 
                                          write_output = NULL) # add va





gc_heatmap(peak_belac_aligned)

plot(peak_belac_aligned, which_plot = "all")

print(peak_belac_aligned)

PEAKDF <- as.data.frame(peak_belac_aligned)
PEAKDF

write.csv(PEAKDF, "peakdf.csv")

#Normalise the data


norm_belac <- norm_peaks(peak_belac_aligned, conc_col_name = "Area", rt_col_name = "RT", out = "data.frame")
norm_belac <- log(norm_belac + 1)
norm_belac


#Change the number of peaks to categorical names so the heatmap gets easier to read
write.csv (norm_belac , "norm_belac.csv")

norm_belac2 <- norm_belac
norm_belac
########Multivariate analysis
adonis <- vegan::adonis2(norm_belac ~  metadata$Treatment, permutations = 999)
adonis


####NMDS

##Keep order of rows consistent###

norm_belac <- norm_belac[match(row.names(metadata), row.names(norm_belac)), ] 
norm_belac


## NMDS using Bray-Curtis dissimilarities
norm_nmds <- vegan::metaMDS(comm = norm_belac, distance = "bray")
norm_nmds

##I'm getting a warning, so I'll explore the stress plot
vegan::stressplot(norm_nmds)


## get x and y coordinates
norm_nmds <- as.data.frame(norm_nmds[["points"]])  
## add the colony as a factor to each sample
norm_nmds <- cbind(norm_nmds,Treatment = metadata[["Treatment"]])


## ordiplot with ggplot2
##This is not working, I need to use a different analysis
library(ggplot2)
ggplot(data = norm_nmds,aes(MDS1,MDS2)) +
  geom_point() + 
  theme_void() + 
  scale_color_manual(values = c("blue","red")) +
  theme(panel.background = element_rect(colour = "black", size = 1.25,
                                        fill = NA), aspect.ratio = 1, legend.position = "none")


###The NMDS analysis did not work, so I will try the heatmap



library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap


norm <-  read_excel ("norm_apcon.xlsx")

norm_a <- as.data.frame(norm)

norm_ap <- norm_belac[,-1]
rownames(norm_ap) <-  norm_belac[,1]
norm_ap

df <-  norm_ap
matdf <- as.matrix(df)

pheatmap(matdf, scale = "column")

#create data frame for annotations
dfh<-data.frame(sample=as.character(rownames(matdf)),dex="Treatment")

dfh$dex<-ifelse(rownames(dfh) %in% c("b14a", "b2a", "b5a", "b6a", "b9a", "b3c", "b10c", "b12c", "b14c", "b15c")) 

pheatmap(matdf, annotation_row = dfh, scale = "column")



#create data frame for annotations
dfh<-data.frame(sample=as.character(rownames(matdf)),dex=meta$Treatment)%>%
  column_to_rownames("sample")

dfh

dfh$dex<-ifelse(rownames(dfh) %in% c("b14a", "b2a", "b5a", "b6a", "b9a", "b3c", "b10c", "b12c", "b14c", "b15c"))

pheatmap(matdf, annotation_row = dfh, scale = "column")


#############################################################################
########################Heatmap with only the selected names#################



library(pheatmap) ## for heatmap generation
library(tidyverse) ## for data wrangling
library(ggplotify) ## to convert pheatmap to ggplot2
library(heatmaply) ## for constructing interactive heatmap


norm <-  read_excel ("norm_apcon_names_2.xlsx", sheet = "selected")
meta <- read_excel ("norm_apcon_names_2.xlsx", sheet = "Meta")
meta


norm_a <- as.data.frame(norm)
norm_ap <- norm_a[,-1]
rownames(norm_ap) <-  norm_a[,1]
norm_ap

df <-  norm_ap
matdf <- as.matrix(df)

pheatmap(matdf, scale = "column")


#create data frame for annotations
dfh<-data.frame(sample=as.character(rownames(matdf)),dex=meta$Treatment)%>%
  column_to_rownames("sample")

dfh

dfh$dex<-ifelse(rownames(dfh) %in% c("b14a", "b2a", "b5a", "b6a", "b9a", "b3c", "b10c", "b12c", "b14c", "b15c"))

pheatmap(matdf, annotation_row = dfh, scale = "column", border_color = "NA", show_rownames = FALSE, angle_col = 45, fontsize_col = 11)



