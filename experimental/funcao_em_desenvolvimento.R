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

#---------- Parte 1: calcular taxas de encontro ----------

encounter.rate <- function(mydata, taxon) {
  
  # split data from each sampling station to calculate encounter rates separately
  samplingStations <- sort(unique(mydata$estacao.amostral))
  for(i in 1:length(samplingStations)) {
    assign(paste("station", samplingStations[i], sep=""),
           subset(mydata, estacao.amostral == samplingStations[i]) ) # create subset for each sampling station
    assign(paste("station", samplingStations[i], "$taxon", sep=""), 
           factor(get(paste("station", samplingStations[i], sep=""))[,taxon]) ) # reset factor levels for taxon at each sampling station subset
    
  encounter_rate <- data.frame(matrix(ncol = (2+length(seq(min(get(paste("station", samplingStations[i], sep=""))[,"Ano"]), max(get(paste("station", samplingStations[i], sep=""))[,"Ano"])))),
                                        nrow = length(unique(get(paste("station", samplingStations[i], sep=""))[,taxon]))) )
  
  colnames(encounter_rate <- c("ID", "taxon", sort(unique(seq(min(get(paste("station", samplingStations[i], sep=""))[,"Ano"]), max(get(paste("station", samplingStations[i], sep=""))[,"Ano"]) )))) ) # cria automaticamente os nomes de colunas de anos
  
  names(get(paste("station", samplingStations[i], "_encounter_rate", sep="")) <- c("ID", "taxon", sort(unique(seq(min(get(paste("station", samplingStations[i], sep=""))[,"Ano"]), max(get(paste("station", samplingStations[i], sep=""))[,"Ano"]) )))) )
  
  #colnames(encounter_rate) <- c("ID", "UC", "taxon", sort(unique(seq(min(mydata$Ano), max(mydata$Ano)))))
  encounter_rate$ID <- c(1:nrow(encounter_rate))
  encounter_rate$taxon <- sort(unique(mydata[,taxon]))
  vetor.Ano <- seq(min(mydata$Ano), max(mydata$Ano))
  # passo 2, preencher objeto encounter_rate
  for(i in 1:nrow(encounter_rate))
    for(j in 1:length(vetor.Ano)){
      a <- subset(mydata, mydata[,taxon] == encounter_rate[i,2])
      b <- subset(a, Ano == vetor.Ano[j]) # extrai o ano automaticamente
      cduc <- unique(a$CDUC)
      c <- subset(mydata, CDUC %in% cduc) # effort must be calculated separately per UC
      
      if ( nrow(subset(c, Ano == vetor.Ano[j])) <= 0)  { encounter_rate[i,j+2] <- NA } else {
        if ( nrow(b) == 0)  { encounter_rate[i,j+2] <- 0 
        }
        else {
          encounter_rate[i,j+2] <- round(nrow(b)/(sum(subset(c, Ano == vetor.Ano[j])$esforço, na.rm=TRUE)/10000), 3)
        }
      }}
  #sitename <- deparse(substitute(mydata))
  #assign(paste("encounter_rate", sitename, sep="_"), encounter_rate, .GlobalEnv)
  assign("encounter_rate", encounter_rate, .GlobalEnv)
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

# remove "species" without ID
#removeTheseBinomials <- sort(unique(dadosICMBio$Binomial))[1:31]
#dadosICMBio <- subset(dadosICMBio, !Binomial %in% removeTheseBinomials)


vector.taxon <- vector.taxon[-c(1:7,73:74)] # remover especies "estranhas"
encounter_rate <- subset(encounter_rate, taxon %in% vector.taxon) # keeping only species in species.list
encounter_rate$taxon <- factor(encounter_rate$taxon) # to remove excluded species from factor levels otherwise they will end up as zeros in the paMatrix
encounter_rate <- encounter_rate[,-c(9:10)] # remover coluna "mean"
head(encounter_rate)
