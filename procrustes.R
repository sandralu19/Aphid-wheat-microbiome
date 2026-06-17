

#procrustes

library(vegan)

setwd("~/Man/Phylo_M_filtered-20251013T190827Z-1-001/Phylo_M_filtered")

#Remember this is in the Phylo_M_filtered folder

# Make ordinations
ord_micro <- read.csv("pcoa_distances_procrusts.csv", row.names = 1)
ord_vol <- read.csv("pca_distances_procrusts_above.csv", row.names = 1)
ord_eco <- read.csv("pca_distances_procrusts_ecoplates.csv", row.names = 1)

#2 weeks
ord_below <- read.csv("pca_distances_procrusts_below_2weeks.csv", row.names = 1)
ord_above <-read.csv("pca_distances_procrusts_above2weeks.csv", row.names = 1)



#Check for name consistency
all(rownames(ord_micro) == rownames(ord_vol))

#Make the rownames consistent
ord_vol <- ord_vol[rownames(ord_micro), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2

#1.Microbiome and Volatiles#######

proc <- procrustes(as.matrix(ord_vol), as.matrix(ord_micro))
res_protest <- protest(as.matrix(ord_vol), as.matrix(ord_micro), permutations = 9999)

plot(proc)  # arrows show how samples shift between spaces
summary(res_protest)
res_protest

summary(proc)
plot(proc, type = "text")  # Adds sample labels


#metadata
metadata <- read.csv("meta_proc_all.csv", row.names = 1)

# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_vol), as.matrix(ord_micro))

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "salmon", "steelblue", "skyblue")[as.factor(metadata$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # volatile
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # micro
# Optional legend
legend("topleft", legend = unique(metadata$Treatment), col = unique(cols), pch = 19, bty = "n")

#Jacknife

set.seed(123)
n_samples <- nrow(ord_vol)
n_boot <- 10000
boot_r <- numeric(n_boot)

for (i in 1:n_boot) {
  # Sample indices with replacement
  idx <- sample(1:n_samples, replace = TRUE)
  
  # Subset ordinations
  micro_sub <- ord_micro[idx, ]
  vol_sub <- ord_vol[idx, ]
  
  # Procrustes on subset
  res_sub <- procrustes(micro_sub, vol_sub)
  
  # Store correlation
  boot_r[i] <- cor(res_sub$Yrot, res_sub$X, method = "pearson")
}

# Summary
mean_r <- mean(boot_r)
ci <- quantile(boot_r, probs = c(0.025, 0.975))
cat("Bootstrap mean Procrustes correlation:", round(mean_r, 3), "\n")
cat("95% CI:", round(ci[1], 3), "-", round(ci[2], 3), "\n")

#plot procrustes
cols <- c("tomato","steelblue")  # example colors for groups

plot(proc, type = "n")  # set up the plot
points(proc$X[,1], res_pro$X[,2], col = cols, pch = 16, cex = 1.5)  # target
points(proc$Yrot[,1], proc$Yrot[,2], col = cols, pch = 1, cex = 1.5)  # rotated

# Draw arrows
arrows(proc$Yrot[,1], proc$Yrot[,2], proc$X[,1], proc$X[,2], length = 0.1, col = "gray50")
legend("topright", legend = c("Microbiome", "Volatiles"), pch = c(16,1), col = c("tomato","steelblue"))


#histogram
hist(boot_r, breaks = 50, col = "skyblue", main = "Bootstrap Procrustes correlation", xlab = "Correlation r")
abline(v = mean_r, col = "red", lwd = 2, lty = 2)
abline(v = ci[1], col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci[2], col = "darkgreen", lwd = 2, lty = 3)
legend("topright", legend = c("Mean", "95% CI"), col = c("red","darkgreen"), lty = c(2,3), lwd = 2)


hist(boot_r, breaks = 50, col = "skyblue", main = "Bootstrap Procrustes correlation", xlab = "Correlation r")
abline(v = mean_r, col = "red", lwd = 2, lty = 2)
abline(v = ci[1], col = "darkgreen", lwd = 2, lty = 3)
abline(v = ci[2], col = "darkgreen", lwd = 2, lty = 3)
legend("topright", legend = c("Mean", "95% CI"), col = c("red","darkgreen"), lty = c(2,3), lwd = 2)

library(vegan)
library(ggplot2)
library(patchwork)  # for combining plots


# Define colors for groups (adjust to your data)
cols <- c("tomato","steelblue","gold","forestgreen")  # example for 4 groups





#ggplot
library(ggplot2)

df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = metadata$Treatment
)

mivo <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59","tomato",  "#91bfdb", "steelblue"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")

mivo

ggsave("above_micro_updated.pdf", plot = mivo, 
       width = 10,   # adjust width as needed
       height = 5,   # adjust height as needed
       device = "pdf")


library(vegan)
mantel(dist(ord_vol), dist(ord_micro), method = "spearman", permutations = 9999)

#####Jackknige leave-one-out

library(vegan)
jack_r <- sapply(1:4, function(i){
  idx <- setdiff(1:4, i)
  mantel(vegdist(ord_vol[idx,]), vegdist(ord_micro[idx,]), permutations=999)$statistic
})
range(jack_r); mean(jack_r)

jack_r <- sapply(1:16, function(i){
  idx <- setdiff(1:16, i)  # <-- this is where idx is defined
  mantel(dist(ord_vol[idx, ], method = "euclidean"),
         dist(ord_micro[idx, ], method = "euclidean"),
         permutations = 999)$statistic
})


jack_r
mean(jack_r)
range(jack_r)

##Bootstrapping
# Load required package
library(vegan)

# Example:
# ord_micro and ord_vol are your distance matrices (e.g., Bray–Curtis, Euclidean, etc.)

# Number of bootstrap resamples
n_boot <- 10000
boot_r <- numeric(n_boot)

set.seed(123)  # for reproducibility

for (i in 1:n_boot) {
  # Resample indices with replacement
  idx <- sample(1:nrow(ord_micro), replace = TRUE)
  
  # Run Mantel on resampled data
  boot_r[i] <- mantel(
    dist(ord_micro[idx, idx], method = "euclidean"),
    dist(ord_vol[idx, idx], method = "euclidean"),
    permutations = 0   # no permutation test inside bootstrap
  )$statistic
}

# Compute summary statistics
mean_r <- mean(boot_r)
ci <- quantile(boot_r, probs = c(0.025, 0.975))

cat("Bootstrap mean Mantel r:", round(mean_r, 3), "\n")
cat("95% CI:", round(ci[1], 3), "-", round(ci[2], 3), "\n")


#2.Ecoplates and microbiome#########


#Check for name consistency
all(rownames(ord_micro) == rownames(ord_eco))

#Make the rownames consistent
ord_eco <- ord_eco[rownames(ord_micro), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2


#metadata
metadata <- read.csv("meta_proc_all.csv", row.names = 1)

# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_micro), as.matrix(ord_eco))
res_protest <- protest(as.matrix(ord_micro), as.matrix(ord_eco), permutations = 999)
summary(res_protest)
res_protest

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "salmon", "steelblue", "skyblue")[as.factor(metadata$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # volatile
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # micro
# Optional legend
legend("topleft", legend = unique(metadata$Treatment), col = unique(cols), pch = 19, bty = "n")


#ggplot
library(ggplot2)

df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = metadata$Treatment
)

mieco <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59","tomato",  "#91bfdb", "steelblue"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")

mieco

#3.Eco and Vol#############


#Check for name consistency
all(rownames(ord_vol) == rownames(ord_eco))

#Make the rownames consistent
ord_eco <- ord_eco[rownames(ord_vol), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2

proc <- procrustes(as.matrix(ord_vol), as.matrix(ord_eco))
res_protest <- protest(as.matrix(ord_vol), as.matrix(ord_eco), permutations = 999)
res_protest

plot(proc)  # arrows show how samples shift between spaces
summary(res_protest)

plot(proc, type = "text")  # Adds sample labels

#Check for name consistency
all(rownames(ord_micro) == rownames(ord_eco))

#Make the rownames consistent
ord_eco <- ord_eco[rownames(ord_micro), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2


#metadata
metadata <- read.csv("meta_proc_all.csv", row.names = 1)

# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_vol), as.matrix(ord_eco))
res_protest <- protest(as.matrix(ord_vol), as.matrix(ord_eco), permutations = 999)
summary(res_protest)
res_protest

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "salmon", "steelblue", "skyblue")[as.factor(metadata$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # volatile
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # micro
# Optional legend
legend("topleft", legend = unique(metadata$Treatment), col = unique(cols), pch = 19, bty = "n")


#ggplot
library(ggplot2)

df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = metadata$Treatment
)

ecovol <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59","tomato",  "#91bfdb", "steelblue"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")
ecovol

library(vegan)
mantel(dist(ord_vol), dist(ord_eco), method = "spearman", permutations = 999)


#####Join figures

######################join two figures##################################################
library(patchwork)

# Combine side by side
combined <- mivo + ecovol    # side by side
# OR stack vertically
# combined <- p1 / p2

# Display combined plot
combined

ggsave("combined_procrustes_above_micro_eco.pdf", plot = combined, 
       width = 10,   # adjust width as needed
       height = 5,   # adjust height as needed
       device = "pdf")




##########################2 weeks##############################################################

#4. Above and below volatiles########

#Check for name consistency
all(rownames(ord_above) == rownames(ord_below))

#Make the rownames consistent
ord_above <- ord_above[rownames(ord_below), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2

#Microbiome and Volatiles

proc <- procrustes(as.matrix(ord_above), as.matrix(ord_below))
res_protest <- protest(as.matrix(ord_above), as.matrix(ord_below), permutations = 999)

plot(proc)  # arrows show how samples shift between spaces
summary(res_protest)
res_protest



groups <- meta2$Treatment
cols <- c("steelblue", "tomato")[as.factor(groups)]

plot(res_protest,type = "points", point.col = cols)

# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_above), as.matrix(ord_below))

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "steelblue")[as.factor(meta2$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # aboveground
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # belowground

# Optional legend
legend("topleft", legend = unique(meta2$Treatment), col = unique(cols), pch = 19, bty = "n")


#ggplot
library(ggplot2)

df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = meta2$Treatment
)

abolow <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59", "#91bfdb"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")

abolow



#5. Above volatiles and below non-volatiles####################

#Non-metabolites

df <- read.csv("data_norm_metabolites.csv", row.names = 1)
norm_apcon <- t(df)

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

#GET DISTANCE MATRUX FOR PROCRUSTS
ind <- get_pca_ind(res.pca)  # Extract PCA results for individuals (samples)
head(ind$coord)              # PCA coordinates for each sample

metadata <- meta2

pca_scores <- ind$coord[, 1:2]

write.csv(pca_scores, "pca_distances_procrusts_metabolites.csv")


ord_meta <- pca_scores


#Procrustes
#Check for name consistency
all(rownames(ord_above) == rownames(ord_meta))

#Make the rownames consistent
ord_above <- ord_above[rownames(ord_meta), ]
# Suppose:
# pcoa1 has columns Dim1, Dim2
# pcoa2 has columns Axis1, Axis2

#Microbiome and Volatiles

proc <- procrustes(as.matrix(ord_above), as.matrix(ord_meta))
res_protest <- protest(as.matrix(ord_above), as.matrix(ord_meta), permutations = 999)

plot(proc)  # arrows show how samples shift between spaces
summary(res_protest)
res_protest



groups <- meta2$Treatment
cols <- c("steelblue", "tomato")[as.factor(groups)]

plot(res_protest,type = "points", point.col = cols)

# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_above), as.matrix(ord_meta))

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "steelblue")[as.factor(meta2$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # aboveground
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # belowground

# Optional legend
legend("topleft", legend = unique(meta2$Treatment), col = unique(cols), pch = 19, bty = "n")

#ggplot
library(ggplot2)

df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = meta2$Treatment
)

abolites <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59", "#91bfdb"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")

abolites


#6.Above and microbiome two weeks#########


#micro pcoa two weeks

ord_micro2 <- read.csv("pcoa_distances_2weeks_procrusts.csv", row.names = 1)



proc <- procrustes(as.matrix(ord_above), as.matrix(ord_micro2))
res_protest <- protest(as.matrix(ord_above), as.matrix(ord_micro2), permutations = 999)

plot(proc)  # arrows show how samples shift between spaces
summary(res_protest)
res_protest




# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_above), as.matrix(ord_micro2))

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "steelblue")[as.factor(meta2$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # aboveground
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # belowground

# Optional legend
legend("topleft", legend = unique(meta2$Treatment), col = unique(cols), pch = 19, bty = "n")

######################join two figures##################################################
library(patchwork)

# Combine side by side
combined <- abolow + abolites    # side by side
# OR stack vertically
# combined <- p1 / p2

# Display combined plot
combined

ggsave("combined_procrustes_metabo_2weeks_above_below2.pdf", plot = combined, 
       width = 10,   # adjust width as needed
       height = 5,   # adjust height as needed
       device = "pdf")


#7. Below volatiles and microbes


# Step 1: fit the Procrustes object
res_pro <- procrustes(as.matrix(ord_below), as.matrix(ord_micro2))
res_protest <- protest(as.matrix(ord_below), as.matrix(ord_micro2), permutations = 999)
summary(res_protest)
res_protest

# Target (X) and rotated (Y) coordinates
X_coords <- res_pro$X       # aboveground (target)
Y_coords <- res_pro$Yrot    # belowground (rotated)

# Example: create color vector from groups
cols <- c("tomato", "salmon", "steelblue", "skyblue")[as.factor(metadata$Treatment)]

# Base plot of arrows
plot(res_pro, kind = 1)  # draws arrows

# Add colored points
points(X_coords, pch = 19, col = cols, cex = 1.5)   # volatile
points(Y_coords, pch = 1, col = cols, cex = 1.5)    # micro
# Optional legend
legend("topleft", legend = unique(metadata$Treatment), col = unique(cols), pch = 19, bty = "n")


#ggplot
library(ggplot2)

metadata <- read.csv("metadata_procrusts.csv", row.names = 1)


df <- data.frame(
  X1 = X_coords[,1], X2 = X_coords[,2],
  Y1 = Y_coords[,1], Y2 = Y_coords[,2],
  Group = metadata$Treatment
)

mivol <- ggplot(df) +
  geom_segment(aes(x = X1, y = X2, xend = Y1, yend = Y2, color = Group),
               arrow = arrow(length = unit(0.2, "cm")), alpha = 0.7) +
  geom_point(aes(x = X1, y = X2, color = Group), size = 3, shape = 19) +
  geom_point(aes(x = Y1, y = Y2, color = Group), size = 3, shape = 1) +
  theme_minimal() +
  coord_equal() +
  scale_colour_manual(values = c("#fc8d59" , "#91bfdb"))+
  labs(x = "Dimension 1", y = "Dimension 2") +
  theme(legend.position = "top")

mivol



ggsave("belowvol_micro_weeek2.pdf", plot = mivol, 
       width = 10,   # adjust width as needed
       height = 5,   # adjust height as needed
       device = "pdf")

#Mantel to confirm correlation between distance matrices

library(vegan)
mantel(dist(ord_below), dist(ord_micro2), method = "spearman", permutations = 999)

