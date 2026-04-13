
# 
#  CLUSTER ANALYSIS SU DATI DI ESPRESSIONE GENICA
#  Stesso dataset dello Script 1: confronta i risultati!
# 

set.seed(123)

# 1. Stessa matrice di espressione (copia da Script 1) 
n_geni   <- 200
campioni <- c(paste0("CTRL_",   1:20),
              paste0("TRATT_A_", 1:20),
              paste0("TRATT_B_", 1:20))
gruppo   <- rep(c("CTRL", "TRATT_A", "TRATT_B"), each = 20)

expr <- matrix(rnorm(n_geni * 60, mean = 6, sd = 1),
               nrow = n_geni, ncol = 60)
rownames(expr) <- paste0("Gene", 1:n_geni)
colnames(expr) <- campioni
expr[1:20,  21:40] <- expr[1:20,  21:40] + 2.5
expr[21:40, 41:60] <- expr[21:40, 41:60] - 2.0

#  2. Clustering gerarchico sui campioni 

#  Distanza: 1 - correlazione di Pearson (tipica per RNA-seq)
#  Misura la similarità del profilo, non della magnitudo
dist_camp <- as.dist(1 - cor(t(t(expr))))
#  Equivalente più leggibile: cor(expr) = correlazione tra campioni
dist_camp <- as.dist(1 - cor(expr))

hc_camp <- hclust(dist_camp, method = "complete")
#  method = "complete": linkage completo (massima distanza tra cluster)
#  Alternative: "ward.D2" (minimizza varianza intra), "average" (UPGMA)

# Dendrogramma campioni — i rami raggruppano CTRL, TRATT_A, TRATT_B?
colori_lab <- c(CTRL    = "#888780",
                TRATT_A = "#1D9E75",
                TRATT_B = "#7F77DD")
col_lab <- colori_lab[gruppo][hc_camp$order]

plot(hc_camp, labels = FALSE,
     main = "Dendrogramma campioni (dist. correlazione)",
     xlab = "", sub = "")
# Coloriamo le etichette manualmente (base R trick)
text(x = 1:60, y = -0.005,
     labels = gruppo[hc_camp$order],
     col    = col_lab,
     srt    = 90, adj = 1, cex = .6, xpd = TRUE)

#  3. Taglio del dendrogramma con K = 3 cluster 
K      <- 3
etich  <- cutree(hc_camp, k = K)
table(etich, gruppo)  # quanti campioni per gruppo finiscono in ogni cluster?

#  4. CLUSTERING GERARCHICO sui GENI (bicluster) 
#  Distanza euclidea su z-score per riga (ogni gene standardizzato)
expr_z   <- t(scale(t(expr)))   # z-score per gene (riga)
dist_gen <- dist(expr_z[1:50, ]) # usiamo i primi 50 geni per velocità
hc_gen   <- hclust(dist_gen, method = "ward.D2")

#  5. Heatmap bicluster (base R) 
#  Riordina righe e colonne secondo i dendrogrammi
heatmap(expr_z[1:50, ],
        Rowv   = as.dendrogram(hc_gen),
        Colv   = as.dendrogram(hc_camp),
        col    = colorRampPalette(c("#534AB7","white","#1D9E75"))(50),
        scale  = "none",        # già fatto lo z-score sopra
        labCol = gruppo,
        cexCol = .6,
        main   = "Heatmap bicluster — z-score espressione")

# 6. k-means sui campioni (alternativa al gerarchico) 
#  Usiamo le prime 5 PC come spazio ridotto (meno rumore)
pca_km  <- prcomp(t(expr), center = TRUE, scale. = TRUE)
spazio  <- pca_km$x[, 1:5]   # campioni × 5 PC

km      <- kmeans(spazio, centers = 3, nstart = 25) # NB 25 (slide 19)
#  nstart = 25: prova 25 inizializzazioni casuali, prende la migliore
table(km$cluster, gruppo)

# Silhouette score per valutare la qualità dei cluster
library(cluster)
sil <- silhouette(km$cluster, dist(spazio))
cat("Silhouette medio:", round(mean(sil[, 3]), 3), "\n")
plot(summary(sil)$clus.avg.widths,
     type = "b", pch = 19,
     xlab = "Cluster", ylab = "Silhouette medio",
     main = "Qualità dei cluster (k-means, K=3)",
     col  = "#7F77DD")
# dendrogramma per verificare se i campioni si raggruppano naturalmente per condizione; 
# heatmap bicluster per visualizzare pattern di espressione; k-means nello spazio PCA ridotto; 
# silhouette per valutare la coerenza interna dei cluster.

# Per scegliere k
#loop su K, non sui cluster a K fisso
spazio <- pca_km$x[, 1:5]   # stessa matrice dello script 2
K_max  <- 8
sil_medio <- numeric(K_max)

for (k in 2:K_max) {
  km_k <- kmeans(spazio, centers = k, nstart = 25)
  sil  <- silhouette(km_k$cluster, dist(spazio))
  sil_medio[k] <- mean(sil[, 3])
}

plot(2:K_max, sil_medio[2:K_max],
     type = "b", pch = 19, col = "#7F77DD",
     xlab = "Numero di cluster K",
     ylab = "Silhouette medio",
     main = "Scelta di K — silhouette al variare di K")
abline(v = which.max(sil_medio[2:K_max]) + 1,
       lty = 2, col = "firebrick")
which.max(sil_medio[2:K_max])+1

# Alternativa con factoextra (un comando solo):
# library(factoextra)
# fviz_nbclust(spazio, kmeans, method = "silhouette", k.max = 8)

