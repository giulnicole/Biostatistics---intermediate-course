# 
#  PCA SU DATI DI ESPRESSIONE GENICA
#  Struttura: 200 geni × 60 campioni (20 CTRL, 20 TRATT_A, 20 TRATT_B)
# 

set.seed(123)

# NB! Creazione dati: non fa parte dell'esercizio, ma questo serve per avere dati congrui
# se avete dei dati non serve

#  1. Simuliamo una matrice di espressione (log2-CPM) 
#        Nella pratica: readRDS("matrice_vst.rds") o simile
n_geni    <- 200
campioni  <- c(paste0("CTRL_",   1:20),
               paste0("TRATT_A_", 1:20),
               paste0("TRATT_B_", 1:20))
gruppo    <- rep(c("CTRL", "TRATT_A", "TRATT_B"), each = 20)

# Matrice base: rumore gaussiano (geni non differenziali)
expr <- matrix(rnorm(n_geni * 60, mean = 6, sd = 1),
               nrow = n_geni, ncol = 60)
rownames(expr) <- paste0("Gene", 1:n_geni)
colnames(expr) <- campioni

# Aggiunge segnale biologico reale ai primi 40 geni
expr[1:20,  21:40] <- expr[1:20,  21:40] + 2.5  # up in TRATT_A
expr[21:40, 41:60] <- expr[21:40, 41:60] - 2.0  # down in TRATT_B

#  2. PCA 
#  prcomp vuole campioni sulle righe, trasponiamo
#  scale. = TRUE: standardizza ogni gene (media 0, sd 1)
#  senza scale: geni con varianza alta dominerebbero
pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)

# Varianza spiegata da ogni PC
var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
print(round(var_exp[1], 1))
print(round(var_exp[2], 1))

# 3. Scree plot (varianza spiegata) 
barplot(var_exp[1:10],
        names.arg = paste0("PC", 1:10),
        ylab = "Varianza spiegata (%)",
        main = "Scree plot",
        col  = "#5DCAA5")

# 4. Score plot: PC1 vs PC2 
#  pca$x contiene le coordinate dei campioni nello spazio PC
df <- data.frame(
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  gruppo = gruppo
)
colori <- c(CTRL    = "#888780",
            TRATT_A = "#1D9E75",
            TRATT_B = "#7F77DD")

plot(df$PC1, df$PC2,
     col  = colori[df$gruppo],
     pch  = 19,
     cex  = 1.4,
     xlab = paste0("PC1 (", round(var_exp[1],1), "%)"),
     ylab = paste0("PC2 (", round(var_exp[2],1), "%)"),
     main = "Score plot PCA (campioni)")
legend("topright", legend = names(colori),
       col = colori, pch = 19, bty = "n")

# 5. Loading plot: quali geni guidano PC1? 
#  pca$rotation: contributo di ogni gene alle PC (eigengene)
loadings_pc1 <- pca$rotation[, 1]
top10        <- sort(abs(loadings_pc1), decreasing = TRUE)[1:10]

barplot(loadings_pc1[names(top10)],
        horiz = TRUE, las = 1,
        xlab  = "Loading su PC1",
        main  = "Top 10 geni — Loading plot",
        col   = ifelse(loadings_pc1[names(top10)] > 0,
                       "#1D9E75", "#7F77DD"))
abline(v = 0, lty = 2, col = "gray50")

# scree plot per scegliere quante PC conservare; 
# score plot per vedere se i campioni si separano per condizione; 
# loading plot per capire quali geni guidano la separazione. 
# Nessuna etichetta discreta assegnata: 
# la struttura è continua.


