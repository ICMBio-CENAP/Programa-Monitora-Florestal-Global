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

# para os testes usar dados de Cazumba-Iracema
cazumba <- subset(dadosICMBio, nome.UC == "Resex Cazumbá-Iracema")
head(cazumba)
dim(cazumba)

#---------- Parte 1: calcular taxas de encontro ----------

encounter.rate <- function(mydata, taxon) {
  mydata[,taxon] <- factor(mydata[,taxon])
  # passo 1, criar objeto mydata2
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata[,taxon]))) )
  colnames(mydata2) <- c("ID", "taxon", sort(unique(seq(min(mydata$Ano), max(mydata$Ano))))) # cria automaticamente os nomes de colunas de anos
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

encounter.rate(cazumba, "Genero")
#encounter_rate_cazumba
encounter_rate

#---------- Parte 2: modelo bayesiano ----------

state.space.model <- function(y, n.years) {
  
  # Specify model in BUGS language
  sink(here("experimental", "ssm.jags"))
  cat("
model { 
# Priors and constraints
N.est[1] ~ dunif(0, 2)            # Prior for initial encounter rate
mean.lambda ~ dunif(0, 10)          # Prior for mean growth rate
sigma.proc ~ dunif(0, 10)           # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 10)           # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   lambda[t] ~ dnorm(mean.lambda, tau.proc) 
   N.est[t+1] <- N.est[t] * lambda[t] 
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(N.est[t], tau.obs)
   }
}
",fill = TRUE)
  sink()
  
  # Bundle data
  jags.data <- list(y = y, T = n.years)
  
  # Initial values
  inits <- function(){list(sigma.proc = runif(1, 0, 5), mean.lambda = runif(1, 0.1, 2), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 0, 2), rep(NA, (n.years-1))))} 
  
  # Parameters monitored
  parameters <- c("lambda", "mean.lambda", "sigma2.obs", "sigma2.proc", "N.est")
  
  # MCMC settings
  ni <- 25000
  nt <- 3
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R (BRT <1 min)
  ssm <- jags(jags.data, inits, parameters, here("experimental", "ssm.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # ccheck results
  print(ssm, digits = 2)
  
  # Draw figure
  fitted <- lower <- upper <- numeric()
  year <- 2014:2019
  n.years <- length(3:ncol(encounter_rate))
  
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$BUGSoutput$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.975)
  }
  m1 <- min(c(fitted, y, lower), na.rm = TRUE)
  m2 <- max(c(fitted, y, upper), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  #plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  plot(0, 0, ylim = c(m1-0.5, m2+1), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:n.years, labels = year)
  polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  points(y, type = "l", col = "black", lwd = 1, lty = 2)
  points(fitted, type = "l", col = "blue", lwd = 2)
  legend(x = 4.5, y = m2+0.5, legend = c("Observada", "Estimada"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "blue"), bty = "n", cex = 0.8)
  
  # Probability of N(2019) < N(2014)
  mean(ssm$BUGSoutput$sims.list$N.est[,6] < ssm$BUGSoutput$mean$N.est[1])
}

# definir numero de anos
n.years <- length(3:ncol(encounter_rate))

# Alouatta
y  <- as.numeric(encounter_rate[3, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Alouatta_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()


# Dasyprocta
y  <- as.numeric(encounter_rate[11, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Dasyprocta_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()


# Pecari
y  <- as.numeric(encounter_rate[27, 3:ncol(encounter_rate)])
# save jpeg
jpeg(here("experimental", "Pecari_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()

# Mazama
y  <- as.numeric(encounter_rate[17, 3:ncol(encounter_rate)])
# save jpeg
jpeg(here("experimental", "Mazama_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()

# Sapajus
y  <- as.numeric(encounter_rate[34, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Sapajus_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()

# Cebus
y  <- as.numeric(encounter_rate[8, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Cebus_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()

# Penelope
y  <- as.numeric(encounter_rate[28, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Penelope_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()

# Crypturellus
y  <- as.numeric(encounter_rate[10, 3:ncol(encounter_rate)])
state.space.model(y, n.years)
# save jpeg
jpeg(here("experimental", "Crypturellus_cazumba.jpg"), width=1000, height=600, res=120) # Open jpeg file
state.space.model(y, n.years)
dev.off()


