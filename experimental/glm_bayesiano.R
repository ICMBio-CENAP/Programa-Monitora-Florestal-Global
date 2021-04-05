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
encounter_rate <- encounter_rate[,-9] # remover coluna "mean"
head(encounter_rate)

# converter para numero inteiro para rodar com poisson
# para mais tarde: mudar modelo para outra distribuicao
encounter_rate[,3:8] <- round(encounter_rate[,3:8]*100)

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
#---------- A fazer: outros modelos di libro BPA with WinBUGS ---------

