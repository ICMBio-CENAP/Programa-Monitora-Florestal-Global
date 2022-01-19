# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyr)
library(ggplot2)
library(dplyr)
library(R2jags)

# carregar funcoes
source(here("bin", "lpi_icmbio.R"))

# carregar dados
dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2019.rds"))
head(dadosICMBio)

# check species
sort(unique(dadosICMBio$Binomial))


#---------- Parte 1: calcular taxas de encontro ----------

encounter.rate <- function(mydata, taxon) {
  mydata[,taxon] <- factor(mydata[,taxon])
  # passo 1, criar objeto encounter_rate
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata[,taxon]))) )
  #mydata2 <- data.frame(matrix(ncol = (3+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata[,taxon]))) )
  colnames(mydata2) <- c("ID", "taxon", sort(unique(seq(min(mydata$Ano), max(mydata$Ano))))) # cria automaticamente os nomes de colunas de anos
  #colnames(mydata2) <- c("ID", "UC", "taxon", sort(unique(seq(min(mydata$Ano), max(mydata$Ano)))))
  mydata2$ID <- c(1:nrow(mydata2))
  mydata2$taxon <- sort(unique(mydata[,taxon]))
  vetor.Ano <- seq(min(mydata$Ano), max(mydata$Ano))
  # passo 2, preencher objeto mydata2
  for(i in 1:nrow(mydata2))
    for(j in 1:length(vetor.Ano)){
      a <- subset(mydata, mydata[,taxon] == mydata2[i,2])
      b <- subset(a, Ano == vetor.Ano[j]) # extrai o ano automaticamente
      cduc <- unique(a$CDUC)
      c <- subset(mydata, CDUC %in% cduc) # effort must be calculated separately per UC
      
      if ( nrow(subset(c, Ano == vetor.Ano[j])) <= 0)  { mydata2[i,j+2] <- NA } else {
        if ( nrow(b) == 0)  { mydata2[i,j+2] <- 0 
        }
        else {
          mydata2[i,j+2] <- round(nrow(b)/(sum(subset(c, Ano == vetor.Ano[j])$esforço, na.rm=TRUE)/10000), 3)
        }
      }}
  #sitename <- deparse(substitute(mydata))
  #assign(paste("encounter_rate", sitename, sep="_"), mydata2, .GlobalEnv)
  assign("encounter_rate", mydata2, .GlobalEnv)
}

encounter.rate(dadosICMBio, "Binomial")


# selecionar populações com criterios minimos
# (taxa de avistamento >= 0.5) e pelo menos 3 anos de dados
encounter_rate$mean <- rowMeans(encounter_rate[,c(3:8)], na.rm = TRUE)
encounter_rate$max <- apply(encounter_rate[,c(3:8)], 1, max, na.rm=T)
#use.this <- subset(encounter_rate, mean >= 0.5) # usar somente espécies com taxa de avistamento médio > 0.1
use.this <- subset(encounter_rate, mean >= 0.5 & max <= 10) # taxa de avistamento médio > 0.1 e max < 10
use.this <- use.this[rowSums(is.na(use.this)) <= 3, ]
vector.taxon <- as.character(factor(use.this$taxon)) # criar vetor com "Binomials" (populações) a serem incluídas na análise
vector.taxon <- vector.taxon[-c(1:7,73:74)] # remover especies "estranhas"
encounter_rate <- subset(encounter_rate, taxon %in% vector.taxon) # keeping only species in species.list
encounter_rate$taxon <- factor(encounter_rate$taxon) # to remove excluded species from factor levels otherwise they will end up as zeros in the paMatrix
encounter_rate <- encounter_rate[,-c(9:10)] # remover coluna "mean"
head(encounter_rate)


# for the analysis we will need covariates, e.g. UC category
# use CDUC in Binomial name to recover categories
# extract CDUC code from taxon (only works if using Binomial for taxon)
encounter_rate$CDUC <- sub(".*_", "", encounter_rate$taxon)
encounter_rate <- merge(encounter_rate, distinct(dadosICMBio, nome.UC, CDUC, .keep_all=FALSE), by="CDUC", all.y=FALSE)
dim(encounter_rate)
# extract category from nome.UC
encounter_rate$nome.UC <- sub(" .*", "", encounter_rate$nome.UC)
names(encounter_rate)[10] <- "PA_type"
head(encounter_rate)

# create dummy covariates just to run model
encounter_rate$dummy1 <- floor(runif(nrow(encounter_rate), min=1, max=10))
encounter_rate$dummy2 <- floor(runif(nrow(encounter_rate), min=1, max=10))

# reorder columns
encounter_rate <- encounter_rate[,c(2:9,1,10:12)]
head(encounter_rate)

# converter para numero inteiro para rodar com poisson
# para mais tarde: mudar modelo para outra distribuicao
encounter_rate[,3:8] <- round(encounter_rate[,3:8]*100)

# test to see if NAs must be replaced by value to run GLM5 (see below)
#encounter_rate[,3:8][is.na(encounter_rate[,3:8])] <- 0 

# plotar serie temporal de todas as especies
matplot(2014:2019, t(encounter_rate[,3:8]), type = "l", lty = 1, lwd = 2, main = "", las = 1, ylab = "Taxa de encontro", xlab = "Year", ylim = c(0, max(encounter_rate[,3:8], na.rm=T)), frame = FALSE)


#---------- Parte 2: modelo bayesiano ----------

#  (a) Null or intercept-only model
# Specify model in BUGS language
sink(here("experimental", "GLM0.jags"))
cat("
model {

# Prior
alpha ~ dnorm(0, 0.01)    # log(mean count)

# Likelihood
for (i in 1:nyear){
   for (j in 1:nsite){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      log.lambda[i,j] <- alpha
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(C = t(encounter_rate[,3:8]), nsite = nrow(encounter_rate), nyear = ncol(encounter_rate[,3:8]))

# Initial values
inits <- function() list(alpha = runif(1, -10, 10))

# Parameters monitored
params <- c("alpha")

# MCMC settings
ni <- 1200
nt <- 2
nb <- 200
nc <- 3

# Call JAGS from R (BRT < 1 min)
out0 <- jags(win.data, inits, params, here("experimental", "GLM0.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out0, dig = 3)


#----------------------------------------------------------------------
#---------- A fazer: outros modelos do livro BPA with WinBUGS ---------

# (h) The full model 
# Specify model in BUGS language
sink(here("experimental", "GLMM5.jags"))
cat("
model {

# Priors
mu ~ dnorm(0, 0.01)                  # Overall intercept
beta1 ~ dnorm(0, 0.01)               # Overall trend 
#beta2 ~ dnorm(0, 0.01)              # First-year observer effect (for ICMBio let us use dummy1 as the effect of a dummy variable)
                                     # actually dummy1 must be a site or a year covariate
                                     # in the first case it affects first value
                                     # in the second case it should affect beta1? 

for (j in 1:ntaxon){
   alpha[j] ~ dnorm(0, tau.alpha)    # Random site effects (for ICMBio it is Binomial_CDUC or taxon by site)
   }
tau.alpha <- 1/ (sd.alpha * sd.alpha)
sd.alpha ~ dunif(0, 3)

for (i in 1:nyear){
   eps[i] ~ dnorm(0, tau.eps)        # Random year effects
   }
tau.eps <- 1/ (sd.eps * sd.eps)
sd.eps ~ dunif(0, 1)

#for (k in 1:nobs){
#   gamma[k] ~ dnorm(0, tau.gamma)   # Random observer effects
#   }
#tau.gamma <- 1/ (sd.gamma * sd.gamma)
#sd.gamma ~ dunif(0, 1)


# Likelihood
for (i in 1:nyear){
   for (j in 1:ntaxon){
      C[i,j] ~ dpois(lambda[i,j])
      lambda[i,j] <- exp(log.lambda[i,j])
      #log.lambda[i,j] <- mu + beta1 * year[i] + beta2 * first[i,j] + alpha[j] + gamma[newobs[i,j]] + eps[i]
      log.lambda[i,j] <- mu + beta1 * year[i] + alpha[j] + eps[i]
      }  #j
   }  #i
}
",fill = TRUE)
sink()

# Bundle data
#win.data <- list(C = t(C), ntaxon = nrow(C), nyear = ncol(C), nobs = 272, newobs = t(newobs), first = t(first), year = ((1:9)-5) / 4)
win.data <- list(C = t(encounter_rate[,3:8]), ntaxon = nrow(encounter_rate),
                 year = 1:6, nyear = ncol(encounter_rate[,3:8]), 
                 dummy1 = encounter_rate$dummy1)

# Initial values
#inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), beta2 = runif(1, -1, 1), alpha = runif(235, -1, 1), gamma = runif(272, -1, 1), eps = runif(9, -1, 1))
inits <- function() list(mu = runif(1, 0, 4), beta1 = runif(1, -1, 1), 
                         alpha = runif(nrow(encounter_rate), -1, 1), eps = runif(ncol(encounter_rate[,3:8]), -1, 1))

# Parameters monitored
#params <- c("mu", "beta1", "beta2", "alpha", "gamma", "eps", "sd.alpha", "sd.gamma", "sd.eps")
params <- c("mu", "beta1", "alpha", "eps", "sd.alpha", "sd.eps")

# MCMC settings
ni <- 12000
nt <- 6
nb <- 6000
nc <- 3

# Call JAGS from R (BRT 11 min)
out7 <- jags(win.data, inits, params, here("experimental", "GLMM5.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out7, dig = 2)


