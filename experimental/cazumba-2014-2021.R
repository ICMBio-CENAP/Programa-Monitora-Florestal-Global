# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyverse)
library(lubridate)
#library(tidyr)
#library(htmlTable)
#library(ggplot2)
#library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
source(here("experimental", "population-functions_cazumba2014-2020.R"))

# para os testes usar dados de Cazumba-Iracema
cazumba <- read_csv(here("data", "cazumba-2014-2021.csv"))

# arrumar os dados
cazumba <- cazumba %>%
  rename(nome_UC = "Local - Nome da Unidade de Conservação",
         estacao_amostral = "Número da Estação Amostral",
         esforco = "Esforço de amostragem tamanho da trilha (m)",
         data = "data da amostragem",
         ano = "Ano",
         classe = "Classe",
         ordem = "Ordem",
         familia = "Família",
         genero = "Gênero",
         especie = "Espécies validadas para análise do ICMBio") %>%
  select(nome_UC, estacao_amostral, esforco, data, ano, classe, ordem,
         familia, genero, especie) %>%
  #mutate(data = mdy(data))
  mutate(data = as.Date(data, "%d/%m/%Y"))
cazumba


# tabela 1: registros por especie
tabela_1 <- cazumba %>%
  filter(classe == "Mammalia") %>%
  group_by(especie) %>%
  count() %>%
  left_join(cazumba %>%
              select(ordem, familia, especie), by="especie") %>%
  distinct(especie, .keep_all = TRUE) %>%
  arrange(ordem, familia, especie) %>%
  select(ordem, familia, especie, n) %>%
  print(n=Inf)


# esforco de amostragem por estacao por dia
cazumba %>%
  filter(! is.na(esforco)) %>%
  distinct(estacao_amostral, data, esforco)

# esforco por ano por estacao
esforco_anual <- cazumba %>%
  group_by(ano, estacao_amostral) %>%
  summarize(esforco = sum(esforco, na.rm = TRUE))
esforco_anual

# numero de registros por ano por estacao
n_registros <- cazumba %>%
  filter(classe == "Mammalia") %>%
  group_by(ano, estacao_amostral, especie) %>%
  count()
n_registros

# taxa de encontro por ano por estacao
taxa_encontro <- left_join(esforco_anual, n_registros, by=c("ano", "estacao_amostral")) %>%
  mutate(taxa_encontro = n/(esforco/1000)) %>%
  select(ano, estacao_amostral, especie, taxa_encontro) %>%
  arrange(especie, ano)
print(taxa_encontro, n=25)

# ajustar tabela acima
#taxa_encontro %>%
#  arrange(desc(especie)) %>%
#  pivot_wider(names_from = ano, values_from = taxa_encontro) #%>%
#  #replace_na()

# grafico da taxa de encontro "naive" separado por estacao amostral
# (teste com dados de Dasyprocta)
dasful <- taxa_encontro %>%
  filter(especie == "Dasyprocta fuliginosa")
dasful

ggplot(data=dasful, aes(x=ano, y=taxa_encontro, group=estacao_amostral)) +
  geom_line(size=0.5, alpha=0.5) +
  geom_point(size=2, alpha=0.5) +
  #geom_smooth(aes(x=ano, y=taxa_encontro), method = "loess") +
  ylim(0, 0.5) +
  #theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  #theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme_classic()

# adicionar linha media
media_anual <- dasful %>% 
  group_by(ano) %>% 
  summarise(media_anual = mean(taxa_encontro))
media_anual

ggplot(dasful, aes(x=ano, y=taxa_encontro, group=estacao_amostral))+ 
  geom_line(size=0.2, alpha=0.5) + 
  #geom_point(size=1, alpha=0.5) +
  geom_line(data = media_anual, aes(ano, media_anual), group = 1, color = "steelblue", size = 3, alpha = 0.5) +
  geom_point(data = media_anual, aes(ano, media_anual), group = 1, color = "steelblue", size = 5) +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  xlab("Ano") + 
  ylab("Taxa de encontro") +
  theme_classic()

# ou smooth...
ggplot(dasful, aes(x=ano, y=taxa_encontro, group=estacao_amostral))+ 
  geom_line(size=0.52, alpha=0.5) + 
  geom_point(size=1, alpha=0.5) +
  geom_smooth(data = media_anual, aes(ano, media_anual), group = 1, color = "steelblue", size = 3, alpha = 0.3) +
  #geom_point(data = media_anual, aes(ano, media_anual), group = 1, color = "steelblue", size = 5) +
  ylim(0, 0.4) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  xlab("Ano") + 
  ylab("Taxa de encontro") +
  theme_classic()

#--------------
# PAREI AQUI!!!

pop.trends.for.rmd <- function(x) { 
  fitted <- lower <- upper <- numeric()
  year <- 2014:2020
  n.years <- length(3:ncol(encounter_rate))
  
  for (i in 1:n.years){
    fitted[i] <- mean(x[,i])
    lower[i] <- quantile(x[,i], 0.025)
    upper[i] <- quantile(x[,i], 0.975)
  }
  m1 <- min(c(fitted, y, lower), na.rm = TRUE)
  m2 <- max(c(fitted, y, upper), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  #plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  #plot(0, 0, ylim = c(m1-0.5, m2+0.5), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  plot(0, 0, ylim = c(0, m2+(mean(fitted)*0.5)), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:n.years, labels = year)
  #polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  #points(y, type = "l", col = "black", lwd = 1, lty = 2)
  #points(fitted, type = "l", col = "blue", lwd = 2)
  points(x = (1:n.years), y = fitted, type = "b", pch = 16, cex = 1.5, lty = 1)
  segments((1:n.years), lower, 1:(n.years), upper, cex=0.5)
}


  summarize(n_registros = sum(esforco, na.rm = TRUE))
n_registros


cazumba <- cazumba %>%
  mutate(ea_dia = paste(estacao_amostral, data, sep="_"))

for(i in 1:length(unique(cazumba$ea_dia))) {
  cazumba[cazumba$ea_dia == unique(cazumba$ea_dia)[i],][1,]$esforco <- 5000
}
cazumba$esforço <- cazumba$esforco
cazumba # check
cazumba$Genero <- gsub( " .*$", "", cazumba$Especie)


# ... e somente mamiferos 
#cazumba <- subset(cazumba, Classe == "Mammalia")
#cazumba[cazumba$Especie == "Dasyprocta cristata", "Especie"] <- "Dasyprocta fuliginosa"
#cazumba[cazumba$Especie == "Dasyprocta leporina", "Especie"] <- "Dasyprocta fuliginosa"
cazumba[cazumba$Especie == "Dasyprocta sp.", "Especie"] <- "Dasyprocta fuliginosa"
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
tabela2 <- cazumba %>% group_by(Genero) %>% summarize(n = n(), ) %>% arrange(Genero)
tabela2
tabela2$'grupos/10km' <- round((tabela2$n/sum(cazumba$esforço, na.rm=T))*10000, 2)
print(tabela2, n=Inf)



#---------- Parte 1: calcular taxas de encontro ----------


encounter.rate(cazumba, "Genero")
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
  Genero = character(),
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

