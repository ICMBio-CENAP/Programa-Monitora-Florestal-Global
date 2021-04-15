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

# criar uma coluna especie por uc
dadosICMBio$species_by_UC <- paste(dadosICMBio$Especie, dadosICMBio$nome.UC)

# para os testes usar somente mamiferos
dadosICMBio <- subset(dadosICMBio, Classe == "Mammalia")

dim(dadosICMBio)

tabela1 <- table(dadosICMBio$nome.UC, dadosICMBio$Ano)
tabela1[tabela1 >1 ] <- 1

#---------- Parte 1: calcular taxas de encontro ----------

encounter.rate(dadosICMBio, "species_by_UC")
#encounter_rate_dadosICMBio
encounter_rate

# selecionar somente pecari tajacu
encounter_rate <- encounter_rate[grep("Pecari tajacu", encounter_rate$taxon), ]

# se taxa avistamento for = 0, alterar para 0.001
encounter_rate[encounter_rate == 0] <- 0.001
#encounter_rate


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
  jpeg(here("report", paste(gsub(" ", "_", pops)[i], "dadosICMBio.jpg", sep="_")), width=1000, height=600, res=120) # Open jpeg file
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
#write.csv(tabela3, here("experimental", "tabela3_dadosICMBio.csv"))
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
#print(xtable(tabela3), type="html", file=here("experimental", "tabela3_dadosICMBio.html"))
#print(xtable(tabela3, type="html")

