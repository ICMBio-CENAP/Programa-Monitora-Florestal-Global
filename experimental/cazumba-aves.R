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
cazumba <- subset(cazumba, Classe == "Aves")
#cazumba[cazumba$Especie == "Dasyprocta cristata", "Especie"] <- "Dasyprocta fuliginosa"
#cazumba[cazumba$Especie == "Dasyprocta leporina", "Especie"] <- "Dasyprocta fuliginosa"
#cazumba[cazumba$Especie == "Dasyprocta sp.", "Especie"] <- "Dasyprocta fuliginosa"
#head(cazumba)
dim(cazumba)

# tabela 1 com esforco por ano
tabela1 <- tibble(
  'Ano' = sort(unique(cazumba$Ano)),
  'Transectos' = 0,
  'Esforço (km)' = 0,
  'Registros' = 0
)
for(i in 1:nrow(tabela1)) {
  df1 <- subset(cazumba, Ano == unique(cazumba$Ano)[i])
  tabela1[i, "Transectos"] <- length(unique(df1$estacao.amostral))
  tabela1[i, "Esforço (km)"] <- sum(df1$esforço, na.rm=T)/1000
  tabela1[i, "Registros"] <- nrow(df1)
}
tabela1


# tabela com especies registradas
# tabela para o anexo
tabela2 <- cazumba %>% group_by(Especie) %>% summarize(n = n(), ) %>% arrange(Especie)
tabela2
tabela2$'grupos/10km' <- round((tabela2$n/sum(cazumba$esforço, na.rm=T))*10000, 2)
print(tabela2, n=Inf)



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
tabela3 <- tibble(
  Especie = character(),
  'r' = numeric(),
  IC = character(),
  'Prob. declinio' = character(),
  'Prob. aumento' = character(),
)

# criar lista para receber N.est
N.est_all_species <- list()

#loop para rodar modelo para todas as spp
for(i in 1:length(species)) {
  y  <- as.numeric(encounter_rate[species[i], 3:ncol(encounter_rate)])
  state.space.model(y, n.years)
  # plot trends
  jpeg(here("report", paste(gsub(" ", "_", species)[i], "cazumba.jpg", sep="_")), width=1000, height=600, res=120) # Open jpeg file
  pop.trends()
  dev.off()
  # save ssm results to be used for making figures later independent of loop and add them to list
  N.est_all_species[[paste(gsub(" ", "_", species)[i], "ssm", sep="_")]] <- ssm$BUGSoutput$sims.list$N.est
  # build table 1
  tabela3[i,] <- list(species[i],
                round(mean(meanR), 2),
                paste(round(quantile(meanR, probs = c(0.025)), 2), " a ", round(quantile(meanR, probs = c(0.975)), 2), sep=""),
                paste((round(length(which(meanR < 0))/length(meanR), 2)*100), "%", sep= ""),
                paste((round(length(which(meanR > 0))/length(meanR), 2)*100), "%", sep= ""))
}
tabela3
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

