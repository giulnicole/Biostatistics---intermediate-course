
# . . . 
# Author: Giulia N Baldrighi
# R EXERCISES: Genetic models and multiple testing correction
#  
# . . . 
# Packages needed:
#   install.packages(c("tidyverse", "ggplot2", "gap", "qvalue"))
#   BiocManager::install("qvalue")  # if using Bioconductor
# . . . 
#
# EXERCISE 1 Simulating genotypes & HWE test

# LEARNING OBJECTIVES:
#   - Understand SNP genotypes and allele frequencies
#   - Simulate genotype data under HWE
#   - Perform HWE chi-square test
#   - Visualise genotype distribution
# . . . 


library(ggplot2)
library(dplyr)
library(tidyverse)

n    <- 1000   # sample size
maf  <- 0.30   # minor allele frequency (T allele)
p    <- 1 - maf  # major allele freq (C)
q    <- maf

# Simulate genotypes under HWE using multinomial
# P(CC) = p², P(CT) = 2pq, P(TT) = q²
probs     <- c(p^2, 2*p*q, q^2)
genotypes <- sample(c("CC","CT","TT"), size = n, replace = TRUE, prob = probs)


obs <- table(genotypes)
print(obs)

# Expected counts
exp_counts <- c(CC = p^2 * n, CT = 2*p*q * n, TT = q^2 * n)
print(round(exp_counts, 1))

# HWE chi-square test (1 df)
chi2_hwe <- sum((obs[c("CC","CT","TT")] - exp_counts)^2 / exp_counts)
p_hwe    <- pchisq(chi2_hwe, df = 1, lower.tail = FALSE)


# if (p_hwe > 0.05): NO significant deviation from HWE (as expected — we simulated under HWE).\n")

# Now we Simulate HWE violation (genotyping error: excess heterozygotes)
bad_genos <- sample(c("CC","CT","TT"), size = n, replace = TRUE,
                    prob = c(0.15, 0.70, 0.15))  # too many CT
obs_bad   <- table(bad_genos)

# Estimate p from allele counts
p_est     <- (2*obs_bad["CC"] + obs_bad["CT"]) / (2*n)
q_est     <- 1 - p_est
exp_bad   <- c(CC = p_est^2 * n, CT = 2*p_est*q_est * n, TT = q_est^2 * n)
chi2_bad  <- sum((obs_bad[c("CC","CT","TT")] - exp_bad)^2 / exp_bad)
p_bad     <- pchisq(chi2_bad, df = 1, lower.tail = FALSE)


# This SNP would be EXCLUDED in QC (p < 1e-6 threshold in large studies)

# Visualisation
df_geno <- data.frame(
  Genotype  = c("CC","CT","TT"),
  Observed  = as.numeric(obs[c("CC","CT","TT")]),
  Expected  = round(exp_counts, 1)
) %>% pivot_longer(cols = c(Observed, Expected), names_to = "Type", values_to = "Count")

ggplot(df_geno, aes(x = Genotype, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Observed" = "#0D7377", "Expected" = "#E67E22")) +
  labs(title = "Exercise 1: Observed vs Expected Genotype Counts (HWE test)",
       subtitle = paste("χ² =", round(chi2_hwe, 3), "  p =", round(p_hwe, 3)),
       y = "Count", x = "Genotype") +
  theme_minimal(base_size = 13)


# . . . 
# EXERCISE 2: Testing genetic models 

# LEARNING OBJECTIVES:
#   - Simulate a SNP with a known additive effect on disease
#   - Fit logistic regression under 4 genetic models
#   - Compare OR estimates across models
#   - Understand which model best captures the truth
# . .  .

n2  <- 2000
maf2 <- 0.25
p2   <- 1 - maf2
q2   <- maf2

# Simulate SNP dosage (0, 1, 2 copies of risk allele)
dosage <- rbinom(n2, 2, q2)

# TRUE model is ADDITIVE: each allele adds log(1.5) to log-odds
# Base disease prevalence ~20%
log_odds <- -1.5 + log(1.5) * dosage
prob_dis  <- plogis(log_odds)
disease   <- rbinom(n2, 1, prob_dis)

df2 <- data.frame(
  dosage  = dosage,
  dom     = as.integer(dosage >= 1),   # dominant:  Aa or aa = 1
  rec     = as.integer(dosage == 2),   # recessive: only aa = 1
  Aa      = as.integer(dosage == 1),   # co-dominant dummy 1
  aa      = as.integer(dosage == 2),   # co-dominant dummy 2
  disease = disease
)

# Genotype distribution
print(table(dosage))

#  Fit all four models
fit_add  <- glm(disease ~ dosage, data = df2, family = binomial)
fit_dom  <- glm(disease ~ dom,    data = df2, family = binomial)
fit_rec  <- glm(disease ~ rec,    data = df2, family = binomial)
fit_cod  <- glm(disease ~ Aa + aa, data = df2, family = binomial)


summary(fit_add)

print(tidy(fit_add, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))



summary(fit_dom)

print(tidy(fit_dom, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))

df2$dom


print(tidy(fit_rec, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))


df2$rec


print(tidy(fit_cod, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))


# Extract key stats
extract_model <- function(fit, model_name, term) {
  b  <- coef(fit)[term]
  se <- summary(fit)$coefficients[term, "Std. Error"]
  pv <- summary(fit)$coefficients[term, "Pr(>|z|)"]
  data.frame(Model = model_name, OR = round(exp(b), 3),
             CI_low = round(exp(b - 1.96*se), 3),
             CI_high = round(exp(b + 1.96*se), 3),
             p_value = round(pv, 4),
             AIC = round(AIC(fit), 1))
}

results <- bind_rows(
  extract_model(fit_add, "Additive",   "dosage"),
  extract_model(fit_dom, "Dominant",   "dom"),
  extract_model(fit_rec, "Recessive",  "rec"),
  extract_model(fit_cod, "Co-dom Aa",  "Aa"),
  extract_model(fit_cod, "Co-dom aa",  "aa")
)


print(results)


aic_df <- data.frame(
  Model = c("Additive","Dominant","Recessive","Co-dominant"),
  AIC   = round(c(AIC(fit_add), AIC(fit_dom), AIC(fit_rec), AIC(fit_cod)), 1)
)
print(aic_df)

# AIC has to be the lowest

# Forest plot of ORs
or_plot <- results %>% filter(Model %in% c("Additive","Dominant","Recessive","Co-dom aa"))

ggplot(or_plot, aes(x = OR, y = Model, color = Model)) +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("#0D7377","#E67E22","#C0392B","#5C35A0")) +
  labs(title = "Exercise 2: OR under the four genetic models",
       subtitle = "True model = additive.",
       x = "Odds Ratio (95% CI)", y = "") +
  theme_minimal(base_size = 13) + theme(legend.position = "none")


# . . . 
# EXERCISE 3: GWAS simulation and multiple testing
# 
# LEARNING OBJECTIVES:
#   - Simulate a small GWAS (1000 SNPs, 1 truly associated)
#   - Apply FDR (Benjamini-Hochberg) correction
#   - Compare uncorrected vs FDR-corrected results
#   - Compute genomic inflation factor λ
#   - Draw QQ plot
# . . . 

n3    <- 500     # individuals
m     <- 1000    # number of SNPs to test
n_causal <- 3   # truly associated SNPs

# Simulate phenotype (quantitative trait, e.g. BMI)
phenotype <- rnorm(n3, mean = 25, sd = 4)

# Simulate m SNPs — all null (no association)
maf_vec <- runif(m, 0.05, 0.45)
geno_mat <- sapply(maf_vec, function(q) rbinom(n3, 2, q))  # n3 x m matrix

# Inject TRUE signal into n_causal randomly chosen SNPs
causal_idx <- sample(1:m, n_causal)
true_betas <- c(1.5, -1.2, 2.0)  # effect sizes (in BMI units per allele)
for (i in seq_along(causal_idx)) {
  phenotype <- phenotype + true_betas[i] * geno_mat[, causal_idx[i]]
}



# Run association tests (linear regression per SNP)
p_vals <- numeric(m)
for (j in 1:m) {
  fit_j   <- lm(phenotype ~ geno_mat[, j])
  p_vals[j] <- summary(fit_j)$coefficients[2, 4]
}

# Multiple testing correction: FDR (Benjamini-Hochberg)
p_adj_bh  <- p.adjust(p_vals, method = "BH")


# Summary of discoveries
alpha     <- 0.05
raw_sig   <- sum(p_vals < alpha)
fdr_sig   <- sum(p_adj_bh < alpha)
causal_detected_raw <- sum(p_vals[causal_idx] < alpha)
causal_detected_fdr <- sum(p_adj_bh[causal_idx] < alpha)

# unadjusted: raw_sig
# FDR-corrected: fdr_sig
# True causal SNPs detected (raw): is the causal_detected_raw
# True causal SNPs detected (FDR): is the causal_detected_fdr

#  Genomic inflation factor λ
chi2_obs   <- qchisq(p_vals, df = 1, lower.tail = FALSE)
lambda_gc  <- median(chi2_obs) / qchisq(0.5, df = 1, lower.tail = FALSE)


# QQ Plot
expected_p <- (1:m) / (m + 1)
qq_df <- data.frame(
  expected = -log10(sort(expected_p)),
  observed = -log10(sort(p_vals)),
  rank     = 1:m
)

ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.5, size = 1, color = "#0D7377") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
  annotate("text", x = 0.5, y = max(qq_df$observed) * 0.95,
           label = paste0("λ = ", round(lambda_gc, 3)), size = 5, color = "#1B2A4A") +
  labs(title = "Exercise 3: QQ Plot — GWAS p-values",
       subtitle = "Points should follow diagonal (red). Departures = inflation or true signal.",
       x = "Expected -log₁₀(p)",
       y = "Observed -log₁₀(p)") +
  theme_minimal(base_size = 13)

# Manhattan-style plot (simplified)
snp_df <- data.frame(
  snp      = 1:m,
  neg_logp = -log10(p_vals),
  adj_p    = p_adj_bh,
  causal   = 1:m %in% causal_idx
)

ggplot(snp_df, aes(x = snp, y = neg_logp, color = causal)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = -log10(0.05/m), linetype = "dashed",
             color = "red", linewidth = 0.8) +           # Bonferroni line
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
             color = "orange", linewidth = 0.8) +         # nominal alpha
  scale_color_manual(values = c("FALSE" = "#AACFE4", "TRUE" = "#C0392B"),
                     labels = c("Null SNP", "Causal SNP")) +
  labs(title = "Exercise 3: Manhattan Plot (simplified GWAS)",
       subtitle = "Red dashed = Bonferroni threshold | Orange dotted = nominal α = 0.05",
       x = "SNP index", y = "-log_10(p)", color = "") +
  theme_minimal(base_size = 13)




# . . . 
# EXERCISE 4: permutation testing
# 
# LEARNING OBJECTIVES:
#   - Understand permutation testing from scratch
#   - Implement family-wise error rate (FWER) control
#   - Compare permutation p-values to Bonferroni
#   - Understand why permutation handles LD better
# . . . 

# smaller problem for speed: 200 SNPs, 300 individuals
n4     <- 300
m4     <- 200
n_perm <- 1000   # number of permutations (use 10,000 in real analyses)

# Simulate genotypes
maf4   <- runif(m4, 0.10, 0.40)
geno4  <- sapply(maf4, function(q) rbinom(n4, 2, q))

# Simulate phenotype — 2 truly associated SNPs
pheno4 <- rnorm(n4)
causal4 <- c(42, 137)
pheno4 <- pheno4 + 1.8 * geno4[, causal4[1]] - 1.5 * geno4[, causal4[2]]

# Step 1: Observed test statistics
obs_stats <- sapply(1:m4, function(j) {
  summary(lm(pheno4 ~ geno4[, j]))$coefficients[2, 3]  # t-statistic
})
obs_pvals <- 2 * pt(-abs(obs_stats), df = n4 - 2)

# Step 2: Permutation distribution of MAX test statistic
perm_max <- numeric(n_perm)

for (b in 1:n_perm) {
  perm_pheno <- sample(pheno4)   # shuffle phenotype labels
  perm_stats <- sapply(1:m4, function(j) {
    abs(summary(lm(perm_pheno ~ geno4[, j]))$coefficients[2, 3])
  })
  perm_max[b] <- max(perm_stats)  # record maximum across all SNPs
}

# Step 3: Permutation p-values (compare each obs stat to null max dist)
perm_pvals <- sapply(abs(obs_stats), function(t) mean(perm_max >= t))

# Step 4: Bonferroni p-values for comparison
bonf_pvals <- pmin(obs_pvals * m4, 1)

result4 <- data.frame(
  SNP         = 1:m4,
  obs_t       = round(obs_stats, 3),
  raw_p       = round(obs_pvals, 4),
  bonf_p      = round(bonf_pvals, 4),
  perm_p      = round(perm_pvals, 4),
  is_causal   = 1:m4 %in% causal4
)

# 10 snps by raw p
print(head(result4[order(result4$raw_p), ], 10))

print(result4[causal4, ])


# false positives
null_idx <- setdiff(1:m4, causal4)
cat("Bonferroni false positives:", sum(result4$bonf_p[null_idx] < 0.05), "\n")
cat("Permutation false positives:", sum(result4$perm_p[null_idx] < 0.05), "\n")

# Visualisation: compare Bonferroni vs Permutation adjusted p-values
compare_df <- data.frame(
  bonf = -log10(pmax(bonf_pvals, 1e-10)),
  perm = -log10(pmax(perm_pvals, 1e-10)),
  causal = 1:m4 %in% causal4
)

ggplot(compare_df, aes(x = bonf, y = perm, color = causal)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "red") +
  scale_color_manual(values = c("FALSE" = "#AACFE4", "TRUE" = "#C0392B"),
                     labels = c("Null SNP", "Causal SNP")) +
  labs(
    title    = "Exercise 4: Bonferroni vs Permutation-Corrected p-values",
    x        = "Bonferroni -log₁₀(p)",
    y        = "Permutation -log₁₀(p)",
    color    = ""
  ) +
  theme_minimal(base_size = 13)

# Distribution of permutation max statistics
perm_df <- data.frame(max_t = perm_max)
causal_t <- max(abs(obs_stats[causal4]))

ggplot(perm_df, aes(x = max_t)) +
  geom_histogram(bins = 50, fill = "#0D7377", alpha = 0.8) +
  geom_vline(xintercept = causal_t, color = "#C0392B", linewidth = 1.2) +
  annotate("text", x = causal_t + 0.1, y = n_perm * 0.03,
           label = paste("Strongest causal\nSNP t-stat =", round(causal_t, 2)),
           color = "#C0392B", hjust = 0, size = 4) +
  labs(title = "Exercise 4: Permutation Null Distribution of Max |t|",
       subtitle = "Red line = observed max t-stat for causal SNP",
       x = "Max |t| across all SNPs (under null)",
       y = "Frequency") +
  theme_minimal(base_size = 13)

