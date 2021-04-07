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

# para os testes usar dados de Cazumba-Iracema
cazumba <- subset(dadosICMBio, nome.UC == "Resex Cazumbá-Iracema")
head(cazumba)
dim(cazumba)

# checar estacoes amostrais
sort(unique(cazumba$estacao.amostral))
sort(unique(cazumba$estacao.amostral))

caz1 <- subset(cazumba, estacao.amostral == 1)
caz2 <- subset(cazumba, estacao.amostral == 2)
caz3 <- subset(cazumba, estacao.amostral == 3)

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

encounter.rate(caz1, "Binomial")
ER_caz1 <- encounter_rate
#ER_caz1$trilha <- 1

encounter.rate(caz2, "Binomial")
ER_caz2 <- encounter_rate
#ER_caz2$trilha <- 2

encounter.rate(caz3, "Binomial")
ER_caz3 <- encounter_rate
#ER_caz3$trilha <- 3

cutia <- rbind(ER_caz1[ER_caz1$taxon=="Dasyprocta_232",], ER_caz2[ER_caz2$taxon=="Dasyprocta_232",], ER_caz3[ER_caz3$taxon=="Dasyprocta_232",])
mean(as.matrix(cutia[,3:8])) # taxa media para todas as trilhas e anos

#-------------------------
# 5.4. Real example: House martin population counts in the village of Magden
# Specify model in BUGS language
sink(here("experimental", "ssm.jags"))
cat("
model {
# Priors and constraints
for (j in 1:nsite) {
  alpha[j] ~ dnorm(mu.alpha, tau.alpha)   # random site effects
  }
mu.alpha ~ dnorm(0, 0.01)
tau.alpha <- 1/(sd.alpha*sd.alpha)
sd.alpha ~ dunif(0,5)

logN.est[1] ~ dnorm(1.27, 0.01)       # Prior for initial population size
mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
  for (j in 1:nsite) {
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t+1] <- logN.est[t] + r[t]
   logN.est[t] <- mean(logN.est.Site[j,t])
  } #t
} #j

# Observation process
for (t in 1:T) {
  for (j in 1:nsite) {
   y[j,t] ~ dnorm(logN.est.Site[j,t], tau.obs)
   logN.est.Site[j,t] <- alpha[j]
   } #t
} #j
  
# Population sizes on real scale
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])
   
   } # t
}

",fill = TRUE)
sink()

# House martin population data from Magden
#pyears <- 6 # Number of future years with predictions
#hm <- c(271, 261, 309, 318, 231, 216, 208, 226, 195, 226, 233, 209, 226, 192, 191, 225, 245, 205, 191, 174, rep(NA, pyears))
year <- 2014:2019

# Bundle data
jags.data <- list(y = log(cutia[,3:8]), T = length(year), nsite=nrow(cutia))

# Initial values
inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1), 
                         logN.est = log(rbind(c(rnorm(1, 1.27, 0.1), rep(NA, (length(year)-1))), c(rnorm(1, 1.27, 0.1), rep(NA, (length(year)-1))), c(rnorm(1, 1.27, 0.1), rep(NA, (length(year)-1))) )) )}

# Parameters monitored
parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")

# MCMC settings
ni <- 200000
nt <- 6
nb <- 100000
nc <- 3

# Call JAGS from R (BRT 3 min)
cutia.ssm <- jags(jags.data, inits, parameters, here("experimental", "ssm.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(cutia.ssm, digits = 3)


# Draw figure
fitted <- lower <- upper <- numeric()
year <- 1990:2015
n.years <- length(hm)
for (i in 1:n.years){
  fitted[i] <- mean(hm.ssm$BUGSoutput$sims.list$N.est[,i])
  lower[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
  upper[i] <- quantile(hm.ssm$BUGSoutput$sims.list$N.est[,i], 0.975)}
m1 <- min(c(fitted, hm, lower), na.rm = TRUE)
m2 <- max(c(fitted, hm, upper), na.rm = TRUE)
par(mar = c(4.5, 4, 1, 1))
plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Population size", xlab = "Year", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
axis(2, las = 1)
axis(1, at = 1:n.years, labels = year)
polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
points(hm, type = "l", col = "black", lwd = 2)
points(fitted, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 150, legend = c("Counts", "Estimates"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 1)

# Probability of N(2015) < N(2009)
mean(hm.ssm$BUGSoutput$sims.list$N.est[,26] < hm.ssm$BUGSoutput$mean$N.est[20])


