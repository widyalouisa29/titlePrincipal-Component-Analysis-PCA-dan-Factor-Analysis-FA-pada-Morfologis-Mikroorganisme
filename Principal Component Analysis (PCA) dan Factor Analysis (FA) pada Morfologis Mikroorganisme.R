data <- read.csv("C:/Semester 4/Analisis multivariat/Dataset/microbes.csv")
data_numeric <- data[, sapply(data, is.numeric)]
head(data_numeric)

selected_features <- c("Solidity", "Eccentricity", "Extent", "EulerNumber", 
                       "BoundingBox2", "BoundingBox4", "MinorAxisLength", 
                       "Perimeter", "ConvexArea", "raddi")
data_selected <- data[, selected_features]
head(data_selected)

#pre-processing
sum(is.na(data_selected))
p <- ncol(data_selected)
summary(data_selected)

#KMO & Bartlett Test
library(psych)
r <- cor(data_selected)
KMO(r)
cortest.bartlett(r, n = nrow(data_selected))  # Gunakan ini untuk menggantikan bartlett.test()

################################ Principal Component Ananlysis ########################

scale_data <- scale(data_selected)
r <- cov(scale_data)
head(scale_data)
head(r)
######## EigenValue and EigenVector ########################
pc <- eigen(r)
head(pc)
pc$values

####### Proporsi Varians and Kumulatif
library(dplyr)
sumvar <- sum(pc$values)
propvar <- sapply(pc$values, function(x) x/sumvar)*100
cumvar <- data.frame(cbind(pc$values, propvar)) %>% mutate(cum = cumsum(propvar))
colnames(cumvar)[1] <- "value"
rownames(cumvar) <- paste0("PC", seq_len(nrow(cumvar)))  # fix rowname length
head(cumvar)

#####PCA Scores ######
pc$vectors
scores <- as.matrix(scale_data) %*% pc$vectors
head(scores)

PCA.mod <- prcomp(scale_data)
summary(PCA.mod)        # vs t(cumvar)
PCA.mod$rotation        # vs pc$vectors
head(PCA.mod$x)         # vs head(scores)

# PCA via FactoMineR
library(FactoMineR)
pca_result <- PCA(scale_data, scale.unit = TRUE, graph = FALSE, ncp = p)

# Ringkasan hasil PCA
pca_result$eig          # vs cumvar
pca_result$svd$V        # vs pc$vectors
pca_result$ind$coord    # vs head(scores)

# Visualisasi
library(factoextra)

# Scree plot
fviz_eig(pca_result, 
         addlabels = TRUE, 
         ncp = p, 
         barfill = "skyblue", 
         barcolor = "darkblue", 
         linecolor = "red")

# Biplot
fviz_pca_biplot(pca_result, 
                geom.ind = "point", 
                addEllipses = TRUE)

# Correlation circle
contrib_circle <- fviz_pca_var(pca_result, col.var = "contrib",
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                               repel = TRUE) + ggtitle("Kontribusi Variabel")
plot(contrib_circle)

# Contribution to PCs
contrib_v_PC1 <- fviz_contrib(pca_result, choice = "var", axes = 1, top = 5) + ggtitle("PC1")
contrib_v_PC2 <- fviz_contrib(pca_result, choice = "var", axes = 2, top = 5) + ggtitle("PC2")
contrib_v_PC3 <- fviz_contrib(pca_result, choice = "var", axes = 3, top = 5) + ggtitle("PC3")
plot(contrib_v_PC1)
plot(contrib_v_PC2)
plot(contrib_v_PC3)

# --------- Factor Analysis -----------
varcov <- cov(scale_data)
pc <- eigen(varcov)

# Eigenvalues dan vectors
pc$values
pc$vectors

# Loadings (dari 3 faktor)
L1 <- sqrt(pc$values[1]) * pc$vectors[, 1]
L2 <- sqrt(pc$values[2]) * pc$vectors[, 2]
L3 <- sqrt(pc$values[3]) * pc$vectors[, 3]

L <- cbind(L1, L2, L3)
print(L)

# FA dengan rotasi varimax
fa_result <- fa(r = scale_data, covar = TRUE, nfactors = 3, rotate = "varimax")
load <- fa_result$loadings

# Plot loadings faktor
plot(load[, c(1, 3)], type = "n")
text(load[, c(1, 3)], labels = names(data_selected), cex = 0.7)

# Diagram faktor
fa.diagram(load)