


#install.packages(c("ggplot2","dplyr","tidyr","patchwork"))
library(ggplot2); library(dplyr)
library(tidyr);   library(patchwork)


 
# (1) Generazione dei dati 

set.seed(530)

n_eur  <- 500          # individui popolazione EUR
n_afr  <- 500          # individui popolazione AFR
n      <- n_eur + n_afr # totale
M      <- 5000         # SNP da testare
M_causal <- 20        # SNP causali (veri positivi)

# Frequenze alleliche (MAF) per ogni SNP nelle due popolazioni
# La differenza nelle frequenze simula la stratificazione
maf_eur <- runif(M, 0.05, 0.50)
delta   <- rnorm(M, mean = 0, sd = 0.15) # diff. differenziazione
maf_afr <- pmin(pmax(maf_eur + delta, 0.01), 0.99)


# matrice genotipica
# Genotipo: 0,1,2 copie dell'allele di rischio (Hardy-Weinberg)
# EUR: ciascun individuo ha 2 draw bernoulliani per ogni SNP
geno_eur <- matrix(
  rbinom(n_eur * M, size = 2, prob = rep(maf_eur, each = n_eur)),
  nrow = n_eur, ncol = M
)
geno_afr <- matrix(
  rbinom(n_afr * M, size = 2, prob = rep(maf_afr, each = n_afr)),
  nrow = n_afr, ncol = M
)
geno <- rbind(geno_eur, geno_afr)  # matrice n × M

# Etichetta di popolazione (covariate "nascosta")
pop <- c(rep("EUR", n_eur), rep("AFR", n_afr))


# simulazione del fenotipo di stratificazione

# Fenotipo binario: la malattia è più comune in AFR (struttura!
# Questo NON è un effetto genetico — è confounding per popolazione)
prob_eur <- 0.30   # prevalenza in EUR
prob_afr <- 0.55   # prevalenza in AFR (confounding)

y_eur <- rbinom(n_eur, 1, prob_eur)
y_afr <- rbinom(n_afr, 1, prob_afr)
y     <- c(y_eur, y_afr)

# Aggiungi effetto reale per M_causal SNP scelti a caso
causal_idx <- sample(1:M, M_causal)
for (i in causal_idx) {
  y <- y + rbinom(n, 1, prob = geno[, i] * 0.08)
}
y <- as.integer(y > 0)  # rimette a binario


# Plot del background genetico


# Test di associazione: regressione logistica semplice, un SNP alla volta
# In pratica equivale a un chi-quadro 1 gdl per grandi n
pvals_naive <- numeric(M)
for (j in 1:M) {
  fit <- glm(y ~ geno[, j], family = binomial())
  pvals_naive[j] <- coef(summary(fit))[2, "Pr(>|z|)"]
}


# . . . 

#  (2) Calcolo λ genomico e QQplot

# Formula: mediana(chi2_obs) / mediana(chi2 atteso sotto H0)
chi2_naive  <- qchisq(1 - pvals_naive, df = 1)
lambda_naive <- median(chi2_naive) / qchisq(0.5, df = 1)

# atteso >> 1 — inflazione da stratificazione

qq_plot <- function(pvals, title = "", lambda = NULL) {
  n_snp   <- length(pvals)
  obs     <- sort(-log10(pvals))
  exp     <- -log10((1:n_snp) / (n_snp + 1))
  df      <- data.frame(exp = rev(exp), obs = obs)
  lbl     <- if (!is.null(lambda))
    paste0(title, "\nλ = ", round(lambda, 3)) else title
  
  ggplot(df, aes(exp, obs)) +
    geom_abline(slope = 1, intercept = 0,
                color = "firebrick", lty = 2, linewidth = .8) +
    geom_point(alpha = .4, size = .7, color = "steelblue") +
    labs(x = "Expected -log10(p)",
         y = "Observed -log10(p)", title = lbl) +
    theme_bw(base_size = 11)
}

p1 <- qq_plot(pvals_naive,
              title = "Senza correzione",
              lambda = lambda_naive)
print(p1)



# . . . 

# (3) PCA e correzione


# Normalizzazione di Price et al. 2006
# x*_ij = (x_ij - 2*p_j) / sqrt(2 * p_j * (1 - p_j))
maf_obs <- colMeans(geno) / 2         # frequenza stimata per ogni SNP
sd_snp  <- sqrt(2 * maf_obs * (1 - maf_obs))

# Rimuovi SNP monomorfi (sd = 0) per evitare divisioni per zero
keep    <- sd_snp > 0
geno_n  <- scale(geno[, keep],
                 center = 2 * maf_obs[keep],
                 scale  = sd_snp[keep])

# PCA via SVD (più efficiente di prcomp per grandi matrici)
# Usiamo solo le prime K componenti — via svds dal pacchetto RSpectra
# Alternativa base R (più lenta ma sempre disponibile):
pca_res <- prcomp(geno_n, center = FALSE, scale. = FALSE, rank. = 10)
PCs     <- pca_res$x   # matrice n × 10

# Plot
df_pca <- data.frame(PC1 = PCs[,1], PC2 = PCs[,2], pop = pop)

ggplot(df_pca, aes(PC1, PC2, color = pop)) +
  geom_point(alpha = .5, size = 1) +
  scale_color_manual(values = c(EUR = "#2196F3", AFR = "#FF5722")) +
  labs(title = "PC1 vs PC2: struttura ancestrale",
       subtitle = "Separazione attesa tra EUR e AFR") +
  theme_bw(base_size = 11)

# Varianza spiegata
var_exp <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
cat("PC1:", round(var_exp[1], 1), "% varianza\n")
cat("PC2:", round(var_exp[2], 1), "% varianza\n")

# Test corretto per le prime k covariate

K <- 2  # tipicamente 10; qui 2 bastano (solo 2 pop)
covars  <- PCs[, 1:K, drop = FALSE]

pvals_corrected <- numeric(M)
for (j in 1:M) {
  fit <- glm(y ~ geno[, j] + covars, family = binomial())
  pvals_corrected[j] <- coef(summary(fit))[2, "Pr(>|z|)"]
}

# Ricalcola λ dopo correzione
chi2_corr     <- qchisq(1 - pvals_corrected, df = 1)
lambda_corr   <- median(chi2_corr) / qchisq(0.5, df = 1)
cat("λ corretto:", round(lambda_corr, 3), "\n")

# Plot a confronto (patchwork)
p1 <- qq_plot(pvals_naive,     "Naive",    lambda_naive)
p2 <- qq_plot(pvals_corrected, "Corretto", lambda_corr)
p1 + p2 + plot_annotation(
  title = "QQ-plot prima e dopo correzione con PCA")


