# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyr)
library(stringr)
library(htmlTable)
library(ggplot2)
library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
source(here("royle-nichols", "population-functions-royle-nichols.R"))

# para os testes usar dados de Cazumba-Iracema
cazumba <- read.csv(here("experimental", "cazumba_2014-2020.csv"))
head(cazumba)

# adicionar ano e corrigir formato da data
cazumba$Ano <- as.numeric(str_sub(cazumba$data.da.amostragem,-4,-1))
head(cazumba)

cazumba$data.da.amostragem <- as.Date(cazumba$data.da.amostragem, format = "%d/%m/%Y")

# adicionar esforco de 5km por sitio-dia
cazumba$site.day <- paste(cazumba$estacao.amostral, cazumba$data.da.amostragem, sep="_")
for(i in 1:length(unique(cazumba$site.day))) {
  cazumba[cazumba$site.day == unique(cazumba$site.day)[i],][1,]$esforco <- 5000
}
cazumba$esforço <- cazumba$esforco
cazumba # check

# como esses dados ainda nao foram validados, ha spp com nomes errados etc
# portanto usar somente genero nesse teste
cazumba$Genero <- gsub( " .*$", "", cazumba$Especie)
head(cazumba)

# ... e somente mamiferos 
#cazumba <- subset(cazumba, Classe == "Mammalia")
#cazumba[cazumba$Especie == "Dasyprocta cristata", "Especie"] <- "Dasyprocta fuliginosa"
#cazumba[cazumba$Especie == "Dasyprocta leporina", "Especie"] <- "Dasyprocta fuliginosa"
#cazumba[cazumba$Especie == "Dasyprocta sp.", "Especie"] <- "Dasyprocta fuliginosa"
#head(cazumba)
#dim(cazumba)

# criar uma coluna especie-trilha
cazumba$spp_trilha <- paste(cazumba$Genero, cazumba$estacao.amostral)


#---------- Parte 1: Historico de deteccao ----------

teste <- detection.history(cazumba)
names(teste)
teste["Dasyprocta"]


# PAREI AQUI!!!
# Ja tem o historico de deteccao para cada especie para cada dia amostrado
# FALTA:
# Separar os anos (estao todos na mesma matriz)
# calcular numero de deteccoes e esforco por estacao e por ano
# rodar o modelo Royle-Nichols :)
# Partes do modelo copiadas abaixo, falta adaptar

#------------------------------------
# ------------------------------------

##----- 1 - Load libraries -----
library(R2jags)
library(here)


##----- 2 - Write models -----
# Specify the model in JAGS language
sink(here("royle-nichols", "RNmodel.txt"))
cat("
   model{
   
   ## Prior distributions
   
   ## Coefficients
   
   # mean value
   # parameter related to abundance
   mu.a0 ~ dnorm(0,0.001)	# intercept
   #mu.a1 ~ dnorm(0,0.001)	# covar

   # parameter related to detectability
   mu.r0	~ dnorm(0,0.001) # intercept
   #mu.r1 ~ dnorm(0,0.001)  # covar

   # standard deviation
   # hyper parameters for abundance
   sigma.a0 ~ dunif(0,10)
   sigma.a1 ~ dunif(0,10)
   sigma.a2 ~ dunif(0,10)
   sigma.a3 ~ dunif(0,10)
   #sigma.a4 ~ dunif(0,10)
   
   # hyper parameters for detectability
   sigma.r0 ~ dunif(0,10)
   sigma.r1 ~ dunif(0,10)
   sigma.r2 ~ dunif(0,10)
   sigma.r3 ~ dunif(0,10)
   #sigma.r4 ~ dunif(0,10)
   
   # create precision
   tau.a0 <- pow(sigma.a0,-2)
   tau.a1 <- pow(sigma.a1,-2)
   tau.a2 <- pow(sigma.a2,-2)
   tau.a3 <- pow(sigma.a3,-2)
   #tau.a4 <- pow(sigma.a4,-2)
   tau.r0 <- pow(sigma.r0,-2)
   tau.r1 <- pow(sigma.r1,-2)
   tau.r2 <- pow(sigma.r2,-2)
   tau.r3 <- pow(sigma.r3,-2)
   #tau.r4 <- pow(sigma.r4,-2)
   
   # parameters related to abundance
   
   a0[i] ~ dnorm(mu.a0,tau.a0) #I(-10,10)
   #a1[i] ~ dnorm(mu.a1,tau.a1) #I(-10,10)

   # parameters related to detectability
   r0[i] ~ dnorm(mu.r0,tau.r0) #I(-10,10)
   #r1[i] ~ dnorm(mu.r1,tau.r1) #I(-10,10)

   ## Likelihood
   
   # loop to estimate the Z matrix (ecological process, abundance for species i at site j)
   for(j in 1:J){
   log(lambda[i,j]) <- a0[i] # + a1[i]*covar1[j]
   
   # detection process model (for species i at site j)
   r[i,j] <- 1/(1+exp(-( r0[i] #+ r1[i]*covar1[j] )))	# from eqn (4.3.1) in Royle and Dorazio, eqn (5) in Yamaura
   p[i,j] <- 1-pow(1-r[i,j])	# from eqn (4.1.1) in Royle and Dorazio, eqn (1) in Yamaura
   X[i,j] ~ dbin(p[i,j],effort[j])	# detection proc in each site
   
    # } # k
    } #j
   
}", fill=TRUE)
sink()


##-----  - Create the necessary arguments to run the jags.model() command -----

# define which data will be used 
jags_data = list( n=n, nzeroes=nzeroes, J=J, X=X,
                  effort=effort, block=block, nblock=nblock, year=year, nyear=nyear,
                  recovery=recovery,  intensity500=intensity500, road500=road500 )


# Specify the parameters to be monitored
params = c('a1', 'a2', 'a3',
           'r1', 'r2', 'r3',
           'A', 'lambda', 'o', 'r', 'N', 'Nsite', 'Nsite2')


inits <- function()list(mu.a0=rnorm(1), mu.a1=rnorm(1), mu.a2=rnorm(1), mu.a3=rnorm(1), #mu.a4=rnorm(1),
                        a.site=matrix(rnorm(1), nrow = n+nzeroes, ncol = nsite),
                        a.year=matrix(rnorm(1), nrow = n+nzeroes, ncol = nyear),
                        mu.r0=rnorm(1), mu.r1=rnorm(1), mu.r2=rnorm(1), mu.r3=rnorm(1), #mu.r4=rnorm(1),
                        sigma.a0=runif(n=1,min=0.5,max=1),
                        sigma.a1=runif(n=1,min=0.5,max=1),
                        sigma.a2=runif(n=1,min=0.5,max=1),
                        sigma.a3=runif(n=1,min=0.5,max=1),
                        #sigma.a4=runif(n=1,min=0.5,max=1),
                        sigma.r0=runif(n=1,min=0.5,max=1), 
                        sigma.r1=runif(n=1,min=0.5,max=1),
                        sigma.r2=runif(n=1,min=0.5,max=1),
                        sigma.r3=runif(n=1,min=0.5,max=1),
                        #sigma.r4=runif(n=1,min=0.5,max=1),
                        #Z = matrix(1, nrow=(n+nzeroes), ncol=J),
                        Z=1+matrix(rpois(n=(n+nzeroes)*T,lambda=0.25),ncol=J,nrow=(n+nzeroes)),
                        psi=runif(1), 
                        w=rbinom((n+nzeroes),1,0.5) )


##----- 4 - Run model and save results -----

# Run the model in jags
fit <- jags(data=jags_data, inits=inits, parameters.to.save=params, n.chains=3, n.iter=500, n.burnin=50, n.thin=10, model.file=here("FORECO", "RNmodel_noK.txt"))
#fit <- jags(data=jags_data, inits=inits, parameters.to.save=params, n.chains=3, n.iter=1000, n.burnin=500, n.thin=30, model.file=here("FORECO", "RNmodel_noK.txt"))

fit <- jags(data=jags_data, inits=inits, parameters.to.save=params, n.chains=3, n.iter=100000, n.burnin=50000, n.thin=100, model.file=here("FORECO","RNmodel_noK.txt"))


#----------------------------
#----------------------------

#encounter.rate(cazumba, "spp_trilha")
#encounter_rate

# selecionar uma especie para teste: Dasyprocta
encounter_rate <- encounter_rate[grep("Dasyprocta", encounter_rate$taxon), ]



# selecionar somente espécies que atendem a criterios minimos
# (taxa de avistamento >= 0.5) e pelo menos 3 anos de dados
names(encounter_rate) # checar numeracao das colunas

# quantas tem NAs?
encounter_rate$na_count <- apply(encounter_rate, 1, function(x) sum(is.na(x)))
encounter_rate <- subset(encounter_rate, na_count < 2)

encounter_rate$mean <- rowMeans(encounter_rate[,c(3:8)], na.rm = TRUE)
#encounter_rate$max <- apply(encounter_rate[,c(3:8)], 1, max, na.rm=T)
#encounter_rate$min <- apply(encounter_rate[,c(3:8)], 1, min, na.rm=T)
#use.this <- subset(encounter_rate, mean >= 0.25 & min >= 0.2) # usar somente espécies com taxa de avistamento médio > 0.1
use.this <- subset(encounter_rate, mean >= 0.25) 
encounter_rate <- use.this
encounter_rate$taxon <- factor(encounter_rate$taxon)
encounter_rate <- encounter_rate[,-c(9:10)]
encounter_rate
#head(encounter_rate)
dim(encounter_rate)


#---------- Parte 2: modelo bayesiano ----------

# check species names for analysis
pops <- encounter_rate$taxon

# criar tabela para receber resultados do loop
tabela3 <- tibble(
  'População' = character(),
  'r' = numeric(),
  IC = character(),
  'Prob. declinio' = character(),
  'Prob. aumento' = character(),
)

# criar lista para receber N.est
N.est_all_pops <- list()

#loop para rodar modelo para todas as spp
for(i in 1:length(pops)) {
  y  <- as.numeric(encounter_rate[pops[i], 3:ncol(encounter_rate)])
  state.space.model(y, n.years)
  # plot trends
  jpeg(here("report", paste(gsub(" ", "_", pops)[i], "cazumba.jpg", sep="_")), width=1000, height=600, res=120) # Open jpeg file
  pop.trends()
  dev.off()
  # save ssm results to be used for making figures later independent of loop and add them to list
  N.est_all_pops[[paste(gsub(" ", "_", pops)[i], "ssm", sep="_")]] <- ssm$BUGSoutput$sims.list$N.est
  # build table 1
  tabela3[i,] <- list(pops[i],
                      round(mean(meanR), 2),
                      paste(round(quantile(meanR, probs = c(0.025)), 2), " a ", round(quantile(meanR, probs = c(0.975)), 2), sep=""),
                      paste((round(length(which(meanR < 0))/length(meanR), 2)*100), "%", sep= ""),
                      paste((round(length(which(meanR > 0))/length(meanR), 2)*100), "%", sep= ""))
}
tabela3
print(tabela3, n=Inf)
#write.csv(tabela3, here("experimental", "tabela3_cazumba.csv"))
class(tabela3)

#library(pander)
#pandoc.table(tabela3, keep.line.breaks = TRUE)

#table2 <- as.data.frame(tabela3)
#table2
#class(table2)

#tabela3.html <- htmlTable(table2, rnames = FALSE)
#tabela3.html
#class(tabela3.html)

# install.packages("xtable")
#library("xtable")
#print(xtable(tabela3), type="html", file=here("experimental", "tabela3_cazumba.html"))
#print(xtable(tabela3, type="html")

