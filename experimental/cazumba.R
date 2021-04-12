# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyr)
library(htmlTable)
library(ggplot2)
library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
source(here("experimental", "population-functions.R"))

# carregar dados
dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2019.rds"))

head(dadosICMBio)

# para os testes usar dados de Cazumba-Iracema
cazumba <- subset(dadosICMBio, nome.UC == "Resex Cazumbá-Iracema")

# ... e somente mamiferos 
cazumba <- subset(cazumba, Classe == "Mammalia")
cazumba[cazumba$Especie == "Dasyprocta cristata", "Especie"] <- "Dasyprocta fuliginosa"
cazumba[cazumba$Especie == "Dasyprocta leporina", "Especie"] <- "Dasyprocta fuliginosa"
cazumba[cazumba$Especie == "Dasyprocta sp.", "Especie"] <- "Dasyprocta fuliginosa"
head(cazumba)
dim(cazumba)

#---------- Parte 1: calcular taxas de encontro ----------


encounter.rate(cazumba, "Especie")
#encounter_rate_cazumba
encounter_rate

# se taxa avistamento for = 0, alterar para 0.001
encounter_rate[encounter_rate == 0] <- 0.001
#encounter_rate

# selecionar somente espécies que atendem a criterios minimos
# (taxa de avistamento >= 0.5) e pelo menos 3 anos de dados
encounter_rate$mean <- rowMeans(encounter_rate[,c(3:8)], na.rm = TRUE)
#encounter_rate$max <- apply(encounter_rate[,c(3:8)], 1, max, na.rm=T)
#encounter_rate$min <- apply(encounter_rate[,c(3:8)], 1, min, na.rm=T)
#use.this <- subset(encounter_rate, mean >= 0.25 & min >= 0.2) # usar somente espécies com taxa de avistamento médio > 0.1
use.this <- subset(encounter_rate, mean >= 0.5) 
encounter_rate <- use.this
encounter_rate$taxon <- factor(encounter_rate$taxon)
#encounter_rate <- encounter_rate[,-c(9:11)]
encounter_rate <- encounter_rate[,-9]
encounter_rate
#head(encounter_rate)
dim(encounter_rate)

encounter_rate$taxon

#---------- Parte 2: modelo bayesiano ----------

# check species names for analysis
species <- encounter_rate$taxon

# criar tabela para receber resultados do loop
table1 <- tibble(
  Especie = character(),
  'r medio' = numeric(),
  IC = character(),
  'Probabilidade de declinio' = character(),
)

#loop para rodar modelo para todas as spp
for(i in 1:length(species)) {
  y  <- as.numeric(encounter_rate[species[i], 3:ncol(encounter_rate)])
  state.space.model(y, n.years)
  # plot trends
  jpeg(here("experimental", paste(species[i], "cazumba.jpg", sep="_")), width=1000, height=600, res=120) # Open jpeg file
  pop.trends()
  dev.off()
  # build table 1
  table1[i,] <- list(species[i],
                round(mean(meanR), 2),
                paste(round(quantile(meanR, probs = c(0.025)), 2), " a ", round(quantile(meanR, probs = c(0.975)), 2), sep=""),
                paste((round(length(which(meanR < 0))/length(meanR), 2)*100), "%", sep= "") )
}
table1
#write.csv(table1, here("experimental", "table1_cazumba.csv"))
class(table1)

#library(pander)
#pandoc.table(table1, keep.line.breaks = TRUE)

#table2 <- as.data.frame(table1)
#table2
#class(table2)

#table1.html <- htmlTable(table2, rnames = FALSE)
#table1.html
#class(table1.html)

# install.packages("xtable")
#library("xtable")
#print(xtable(table1), type="html", file=here("experimental", "table1_cazumba.html"))
#print(xtable(table1, type="html")

