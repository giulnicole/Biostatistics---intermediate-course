# . . . 
# Author: Giulia N Baldrighi
# R EXERCISES: Confounding, interaction and effect modification
#  
# . . . 
# Packages needed:
#   install.packages(c("tidyverse", "epitools", "epiR", 
#                      "survival", "broom", "ggplot2"))
# 
# LEARNING OBJECTIVES:
#   - Check models' assumptions
#   - Understand linear model and logistic model interpretation
#   - Considering confounding and interaction
#   - Understand differences between regressors' effects
#   - Vizualise the effects
# . . . 

library(tidyverse)
library(epitools)
library(broom)
library(ggplot2)

set.seed(142)  # for reproducibility

# Pre exercise . . . . .. . . 

# Recap: linear model and logistic model 
n<- 1000


# Dataset from last time

# Creating treatment groups
# groups
group <- sample(c("Control", "Drug A", "Drug B"), n,
                replace = TRUE, prob = c(0.33, 0.33, 0.34))

# Demographics variables
age        <- round(rnorm(n, mean = 58, sd = 12))
age        <- pmax(pmin(age, 90), 25)   # clamp 25-90

sex        <- sample(c("Male", "Female"), n, replace = TRUE)

bmi        <- round(rnorm(n, mean = 27.5, sd = 4.5), 1)
bmi        <- pmax(pmin(bmi, 45), 16)

smoker     <- sample(c("Yes", "No", "Ex-smoker"), n,
                     replace = TRUE, prob = c(0.25, 0.55, 0.20))

# Comorbidities
diabetes   <- rbinom(n, 1, prob = ifelse(bmi > 30, 0.40, 0.15))
hypertension <- rbinom(n, 1,
                       prob = ifelse(age > 60, 0.55,
                                     ifelse(age > 45, 0.35, 0.15)))


# Continuous clinical variables (with group effect)
group_effect <- ifelse(group == "Drug A", -8,
                       ifelse(group == "Drug B", -14, 0))

systolic_bp  <- round(130 + group_effect * 0.6 +
                        (age - 58) * 0.4 +
                        diabetes * 8 +
                        rnorm(n, 0, 12))

diastolic_bp <- round(80 + group_effect * 0.3 +
                        rnorm(n, 0, 8))

ldl          <- round(130 + group_effect * 1.2 +
                        rnorm(n, 0, 25), 1)

hdl          <- round(50 - group_effect * 0.3 +
                        rnorm(n, 0, 10), 1)

hba1c        <- round(6.2 + diabetes * 1.8 +
                        (bmi - 27.5) * 0.05 +
                        rnorm(n, 0, 0.6), 1)

creatinine   <- round(abs(rnorm(n, 1.1, 0.25)), 2)

# Outcome
hospitalized <- rbinom(n, 1,
                       prob = plogis(-2 + 0.03 * (age - 58) +
                                       0.05 * diabetes +
                                       0.01 * (systolic_bp - 130) +
                                       ifelse(group == "Drug B", -0.5, 0)))

# Quality of life score (0-100)
qol_score <- round(70 + group_effect * 0.8 +
                     -0.3 * (age - 58) +
                     -5  * diabetes +
                     rnorm(n, 0, 10))
qol_score <- pmax(pmin(qol_score, 100), 0)


# here we have generated all the variables, now we can put all together in a dataset:

# Assemble dataset
df <- tibble(
  id            = 1:n,
  group         = factor(group, levels = c("Control", "Drug A", "Drug B")),
  age           = age,
  sex           = factor(sex),
  bmi           = bmi,
  smoker        = factor(smoker, levels = c("No", "Ex-smoker", "Yes")),
  diabetes      = factor(diabetes, labels = c("No", "Yes")),
  hypertension  = factor(hypertension, labels = c("No", "Yes")),
  systolic_bp   = systolic_bp,
  diastolic_bp  = diastolic_bp,
  ldl           = ldl,
  hdl           = hdl,
  hba1c         = hba1c,
  creatinine    = creatinine,
  qol_score     = qol_score,
  hospitalized  = factor(hospitalized, labels = c("No", "Yes"))
)

View(df)


# Recap
# modello lineare (var. outcome continua)


# y = systolic_bp
# x = bmi/age

# precheck sulla var. indipendente# per vedere se adotatre misure di correzione per deviazioni della distrib dalla normalità 
# o presenza di outlier

plot(df$bmi, df$systolic_bp)
plot(df$age, df$systolic_bp)

# normality check outcome (shapito test)
shapiro.test(df$systolic_bp)

# modello lineare
mod1 <- lm(systolic_bp ~ age, data=df)
summary(mod1)

tidy(mod1) |>
  mutate(across(where(is.numeric), \(x) round(x, 1000))) #|>
  #write.csv("mod1_lm_table.csv", row.names = FALSE)


# testiamo ora un modello logistico (var. binaria)
# attenzione alla categoria di riferiemnto

df$sex
df$diabetes

mod2 <- glm(diabetes ~ sex, data=df, family = binomial)
summary(mod2)

tidy(mod2, exponentiate = TRUE, conf.int = TRUE) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) #|>
 # write.csv("mod2_glm_table.csv", row.names = FALSE)


# . . . 
# EXERCISE 1: Detecting Confounding
# . . . 
# SCENARIO:
#   We study whether coffee drinking (E) is associated with
#   pancreatic cancer (O). Smoking (C) is a confounder:
#   smokers drink more coffee AND have higher cancer risk.
#
# We will:
#   1. Simulate data with known confounding structure
#   2. Compute crude OR (unadjusted)
#   3. Compute adjusted OR (controlling for smoking)
#   4. Apply the 10% rule to confirm confounding
# . . . 


n <- 2000

# Simulate smoking (confounder) - 30% prevalence
smoking <- rbinom(n, 1, 0.30)

# Coffee drinking (exposure): smokers drink more coffee
# P(coffee | smoker) = 0.70,  P(coffee | non-smoker) = 0.30
p_coffee <- ifelse(smoking == 1, 0.70, 0.30)
coffee <- rbinom(n, 1, p_coffee)

# Pancreatic cancer (outcome):
# Smoking increases risk (OR ~5), coffee has NO TRUE effect
# log-odds base: -4 (rare disease ~1.8%)
log_odds_cancer <- -4 + log(5) * smoking
p_cancer <- plogis(log_odds_cancer)
cancer <- rbinom(n, 1, p_cancer)

df1 <- data.frame(coffee, smoking, cancer)
View(df1)

#  CRUDE MODEL (unadjusted) 
# H0: X (=coffee) non è associato al tumore/ non ha effetto sul tumore al pancreas -> OR=1
# H1: X (=coffee) è associato al tumore / ha effetto significativo -> OR!=1 
crude_model <- glm(cancer ~ coffee, data = df1, family = binomial)
summary(crude_model)
crude_or <- exp(coef(crude_model)["coffee"])


# ADJUSTED MODEL (controlling for smoking) 

adj_model <- glm(cancer ~ coffee + smoking, data = df1, family = binomial)

summary(adj_model)
summary(crude_model)
adj_or <- exp(coef(adj_model)["coffee"])


#  10% RULE 
pct_change <- abs(crude_or - adj_or) / crude_or * 100


print(tidy(adj_model, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))


# Bar plot: crude vs adjusted OR
or_data <- data.frame(
  Model = c("Crude", "Adjusted"),
  OR    = c(crude_or, adj_or)
)

p1 <- ggplot(or_data, aes(x = Model, y = OR, fill = Model)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("#0D7377", "#1B2A4A")) +
  labs(
    title    = "Exercise 1: crude vs adjusted OR for coffee vs. cancer",
    y        = "Odds Ratio",
    x        = "Model"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p1)



# . . . 
# EXERCISE 2: Interaction in logistic regression
# . . . 
# SCENARIO:
#   Drug A (E) reduces mortality in a patient population.
#   Patients are either high-risk (C=1) or low-risk (C=0).
#   We hypothesise that Drug A works better in high-risk patients
#   -> i.e., there may be a positive INTERACTION (E × C).
#
# We will:
#   1. Simulate data with a built-in interaction
#   2. Fit a model WITHOUT interaction term
#   3. Fit a model WITH interaction term
#   4. Compare models and interpret β₃
#   5. Compute interaction on both scales (additive, multiplicative)
# . . . 

n2 <- 1500

# High-risk group (C): 40% of patients
high_risk <- rbinom(n2, 1, 0.40)

# Drug A (E): randomly assigned (50/50, like an RCT)
drug_A <- rbinom(n2, 1, 0.50)

# Mortality (outcome):
# Base log-odds: -2 (low risk, no drug)
# High risk adds log(3) to odds
# Drug A alone reduces by log(0.7) = small benefit for low risk
# Drug A * high risk adds an extra benefit: log(0.4) (interaction)
log_odds_death <- -2 +
  log(3.0) * high_risk +
  log(0.7) * drug_A +
  log(0.40) * drug_A * high_risk  # <-- TRUE interaction term

p_death <- plogis(log_odds_death)
death <- rbinom(n2, 1, p_death)

df2 <- data.frame(drug_A, high_risk, death)

# MODEL WITHOUT INTERACTION
model_no_int <- glm(death ~ drug_A + high_risk, data = df2, family = binomial)

print(tidy(model_no_int, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))


# MODEL WITH INTERACTION 
model_int <- glm(death ~ drug_A * high_risk, data = df2, family = binomial)

print(tidy(model_int, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))

# Extract interaction term
b3      <- coef(model_int)["drug_A:high_risk"]
b3_or   <- exp(b3)
b3_pval <- summary(model_int)$coefficients["drug_A:high_risk", "Pr(>|z|)"]


# STRATUM-SPECIFIC EFFECTS 
low_df  <- df2 %>% filter(high_risk == 0)
high_df <- df2 %>% filter(high_risk == 1)


or_low  <- exp(coef(glm(death ~ drug_A, data = low_df,  family = binomial))["drug_A"])
or_high <- exp(coef(glm(death ~ drug_A, data = high_df, family = binomial))["drug_A"])

#  RERI (Additive interaction) 
# Need RR estimates from risk ratios rather than OR
# Use log-binomial model for additive scale
rr_model <- glm(death ~ drug_A * high_risk, data = df2, family = poisson(link = "log"))

RR <- exp(coef(rr_model))

RR_10 <- exp(coef(rr_model)["drug_A"])
RR_01 <- exp(coef(rr_model)["high_risk"])
RR_11 <- exp(coef(rr_model)["(Intercept)"] + coef(rr_model)["drug_A"] +
              coef(rr_model)["high_risk"]   + coef(rr_model)["drug_A:high_risk"]) /
         exp(coef(rr_model)["(Intercept)"])

RERI <- RR_11 - RR_10 - RR_01 + 1
RERI

#  Plot for visualizing interaction
pred_data <- expand.grid(drug_A = c(0, 1), high_risk = c(0, 1))
pred_data$pred_prob <- predict(model_int, newdata = pred_data, type = "response")
pred_data$RiskGroup <- ifelse(pred_data$high_risk == 1, "High Risk", "Low Risk")
pred_data$Treatment  <- ifelse(pred_data$drug_A  == 1, "Drug A", "No Drug")

p2 <- ggplot(pred_data, aes(x = Treatment, y = pred_prob,
                             color = RiskGroup, group = RiskGroup)) +
  geom_line(size = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("High Risk" = "#C0392B", "Low Risk" = "#0D7377")) +
  labs(
    title    = "Exercise 2: interaction Drug A x risk group",
    subtitle = "Non-parallel lines indicate interaction",
    x        = "Treatment",
    y        = "Predicted probability of death",
    color    = "Risk group"
  ) +
  theme_minimal(base_size = 13)

print(p2)


# . . . 
# EXERCISE 3: Effect modification by stratification
# . . . 
# SCENARIO:
#   Aspirin use (E) and GI bleeding (O).
#   Sex is a suspected effect modifier: aspirin may carry
#   higher bleeding risk in females than in males.
#
# We will:
#   1. Simulate data with sex as an effect modifier
#   2. Compute overall (crude) OR
#   3. Compute sex-stratified ORs
#   4. Mantel-Haenszel pooled estimate
#   5. Woolf test of homogeneity
# . . . 


n3 <- 3000
sex     <- rbinom(n3, 1, 0.50)   # 0 = male, 1 = female
aspirin <- rbinom(n3, 1, 0.45)   # ~45% on aspirin

# GI bleeding:
# Base log-odds: -3
# Aspirin effect differs by sex:
#   Males:   OR ~1.5  (log = 0.405)
#   Females: OR ~3.5  (log = 1.253)
# This IS effect modification
log_odds_bleed <- -3 +
  log(1.5) * aspirin * (1 - sex) +
  log(3.5) * aspirin * sex

p_bleed <- plogis(log_odds_bleed)
bleed   <- rbinom(n3, 1, p_bleed)

df3 <- data.frame(aspirin, sex, bleed)

crude3 <- glm(bleed ~ aspirin, data = df3, family = binomial)
crude3_or <- exp(coef(crude3)["aspirin"])

#  STRATIFIED ANALYSIS 
male_df   <- df3 %>% filter(sex == 0)
female_df <- df3 %>% filter(sex == 1)

or_male   <- exp(coef(glm(bleed ~ aspirin, data = male_df,   family = binomial))["aspirin"])
or_female <- exp(coef(glm(bleed ~ aspirin, data = female_df, family = binomial))["aspirin"])


#  FORMAL TEST FOR EFFECT MODIFICATION 
model_em <- glm(bleed ~ aspirin * sex, data = df3, family = binomial)
em_pval  <- summary(model_em)$coefficients["aspirin:sex", "Pr(>|z|)"]

# only aspirin
model0 <- glm(bleed ~ aspirin, data = df3, family = binomial)

print(tidy(model0, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))

# 0 female
# 1 male 

# model 2
summary(model_em)
print(tidy(model_em, exponentiate = TRUE, conf.int = TRUE) %>%
        select(term, estimate, conf.low, conf.high, p.value))


#  WOOLF TEST (manual χ² of heterogeneity) 
# log(OR) for each stratum + variance
log_or_m <- log(or_male);   var_m <- sum(1 / table(male_df$aspirin,   male_df$bleed))
log_or_f <- log(or_female); var_f <- sum(1 / table(female_df$aspirin, female_df$bleed))

pooled_log_or <- (log_or_m / var_m + log_or_f / var_f) / (1/var_m + 1/var_f)
Q_stat <- (log_or_m - pooled_log_or)^2 / var_m + (log_or_f - pooled_log_or)^2 / var_f
woolf_p <- pchisq(Q_stat, df = 1, lower.tail = FALSE)


#  stratified ORs 
or_df <- data.frame(
  Group = c("Males", "Females", "Pooled (MH)"),
  OR    = c(or_male, or_female, exp(pooled_log_or)),
  Lower = c(exp(log_or_m - 1.96*sqrt(var_m)),
            exp(log_or_f - 1.96*sqrt(var_f)),
            exp(pooled_log_or - 1.96*sqrt(1/(1/var_m + 1/var_f)))),
  Upper = c(exp(log_or_m + 1.96*sqrt(var_m)),
            exp(log_or_f + 1.96*sqrt(var_f)),
            exp(pooled_log_or + 1.96*sqrt(1/(1/var_m + 1/var_f))))
)

p3 <- ggplot(or_df, aes(x = OR, y = Group)) +
  geom_point(aes(color = Group), size = 4) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Group), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Males" = "#0D7377", "Females" = "#C0392B",
                                "Pooled (MH)" = "#1B2A4A")) +
  labs(
    title    = "Exercise 3: Stratified OR: Aspirin vs. GI bleed by Sex",
    subtitle = "Non-overlapping CIs and interaction p-value confirm effect modification",
    x        = "Odds Ratio (95% CI)",
    y        = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p3)

# . . . 
# Comprehensive visualisation of the three effects
# . . . 
# Using the datasets from Ex 1–3, create three combined plots:
#   A. Coefficient comparison (confounding)
#   B. Interaction plot (Exercise 2)
#   C. Stratified effects summary (Exercise 3)
# These three plots together illustrate how confounding,
# interaction, and effect modification look in practice.
# . . . 

library(gridExtra)

# Panel A: Confounding — crude vs adjusted coefficient
coef_df <- data.frame(
  Model       = factor(c("Crude OR\n(confounded)",
                         "Adjusted OR\n(true effect)"), 
                       levels = c("Crude OR\n(confounded)", "Adjusted OR\n(true effect)")),
  OR          = c(crude_or, adj_or),
  Description = c("Includes confounding\nby smoking",
                  "Smoking controlled\n(true null effect)")
)

pA <- ggplot(coef_df, aes(x = Model, y = OR, fill = Model)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.9) +
  geom_text(aes(label = round(OR, 2)), vjust = -0.4, size = 5, fontface = "bold") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
  scale_fill_manual(values = c("#E67E22", "#0D7377")) +
  ylim(0, max(crude_or, adj_or) + 0.5) +
  labs(title = "A — Confounding", subtitle = "Coffee vs. Cancer OR",
       y = "Odds Ratio", x = "") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")

# Panel B: Interaction — already created as p2 above
pB <- p2 + labs(title = "B — Interaction", subtitle = "Drug A × risk group") +
  theme_minimal(base_size = 11)

# Panel C: Effect modification — already created as p3 above
pC <- p3 + labs(title = "C — Effect modification", subtitle = "Aspirin> GI bleed by Sex") +
  theme_minimal(base_size = 11)

# Combine all three
grid.arrange(pA, pB, pC, ncol = 3,
             top = "Confounding vs interaction vs effect modification")

