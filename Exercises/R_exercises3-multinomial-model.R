# . . . 
# Author: Giulia N Baldrighi
# R EXERCISE: multinomial model
#  


# In R 
# Interpretare un modello multinomiale multinom() significa capire come le variabili esplicative 
# influenzano la probabilità di appartenere a ciascuna categoria rispetto a una categoria di riferimento.


set.seed(123)

n <- 500

# Variabili genetiche (binarie)
VAR1 <- rbinom(n, 1, 0.3)
VAR2 <- rbinom(n, 1, 0.4)
VAR3 <- rbinom(n, 1, 0.2)
VAR4 <- rbinom(n, 1, 0.5)

# Confondenti
eta <- rnorm(n, mean = 60, sd = 10)
sesso <- factor(sample(c("M", "F"), n, replace = TRUE))

# Probabilità per le classi (logit multinomiale)
lp_aggr <- -1 + 0.8*VAR1 + 0.5*VAR2 + 0.03*eta
lp_res  <- -1 + 1.2*VAR3 + 0.7*VAR4 + 0.02*eta

# Trasformazione softmax

# perchè serve questa trasformazione: 
# funzione di attivazione fondamentale utilizzata per trasformare un vettore di numeri reali grezzi 
# (chiamati logits) in una distribuzione di probabilità su più classi mutuamente esclusive

exp_base <- 1
exp_aggr <- exp(lp_aggr)
exp_res  <- exp(lp_res)

prob_base <- exp_base / (exp_base + exp_aggr + exp_res)
prob_aggr <- exp_aggr / (exp_base + exp_aggr + exp_res)
prob_res  <- exp_res  / (exp_base + exp_aggr + exp_res)

# Outcome multinomiale
outcome <- mapply(function(p1, p2, p3) {
  sample(c("Base", "Aggressiva", "Resistente"), 1, prob = c(p1, p2, p3))
}, prob_base, prob_aggr, prob_res)

outcome <- factor(outcome)

# Dataset finale
dati <- data.frame(outcome, VAR1, VAR2, VAR3, VAR4, eta, sesso)


# Ricodifica per avere base come riferimento
dati$outcome <- relevel(dati$outcome, ref = "Base")

library(nnet)
library(car)

modello <- multinom(outcome ~ VAR1 + VAR2 + VAR3 + VAR4 + eta + sesso, data = dati)
summary(modello)

vif(modello)



# test delle varianti (p-valaue)
z <- summary(modello)$coefficients / summary(modello)$standard.errors
p_value <- (1 - pnorm(abs(z), 0, 1)) * 2

p_value

# modello di interazione
# Es interazione tra VAR1 e VAR2
mod_interazione <- multinom(outcome ~ VAR1*VAR2 + VAR3 + VAR4 + eta + sesso, data = dati)

summary(mod_interazione)


# verifica
mod_senza_eta <- multinom(outcome ~ VAR1 + VAR2 + VAR3 + VAR4 + sesso, data = dati)
mod_con_eta   <- multinom(outcome ~ VAR1 + VAR2 + VAR3 + VAR4 + sesso + eta, data = dati)

summary(mod_senza_eta)
summary(mod_con_eta)

