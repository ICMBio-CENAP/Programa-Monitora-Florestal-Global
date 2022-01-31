# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyr)
library(tidyverse)
library(htmlTable)
library(ggplot2)
library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
source(here("experimental", "population-functions.R"))

# carregar dados
#dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2019.rds"))
#head(dadosICMBio)

# para os testes usar dados de tapArap-Iracema
#tapArap <- subset(dadosICMBio, nome.UC == "Resex Tapajós-Arapiuns")

# alternativa: dados brutos enviados por Rubia e Iasmin jan 2022
tapArap <- read_csv(here("data", "Dados_brutos_Mastoaves_TapajosArapiuns_Validadoate2019_YR_final.csv"))

# alguns ajustes
tapArap <- tapArap %>%
  rename(cduc = "CDUC",
         nome.UC = "Local - Nome da Unidade de Conservação",
         estacao.amostral = "Número da Estação Amostral",
         esforco = "Esforço de amostragem tamanho da trilha (m)",
         ano = "Ano",
         data = "data da amostragem",
         hora.inicio = "horário de início  (h:mm)",
         hora.fim = "horário de término (h:mm)",
         classe = "Classe",
         genero = "Gênero",
         binomial = "Espécies validadas para análise do ICMBio",
         nivel.validacao  = "Clasificação taxonômica validada",
         hora.registro = "horário do avistamento",
         n.individuos = "n° de animais",
         dist.perpendicular = "distância (m)     do animal em relação a trilha") %>%
  unite("populacao", cduc, binomial, sep =" ", remove = FALSE) %>%
  mutate(binomial = as_factor(binomial),
         populacao = as_factor(populacao),
         n.individuos = as.numeric(n.individuos),
         dist.perpendicular = gsub(",", ".", dist.perpendicular),
         dist.perpendicular = as.numeric(dist.perpendicular),
         ano = str_sub(ano, start= -4)) %>%
  mutate(ano = as.numeric(ano)) %>%
  filter(ano != 2020) %>%
  select(cduc, nome.UC, estacao.amostral, esforco, ano, data, hora.inicio, hora.fim,
         classe, genero, binomial, populacao, nivel.validacao, hora.registro, n.individuos, dist.perpendicular)
tapArap  

# verificar e corrigir nomes das especies
sort(unique(tapArap$binomial))
tapArap[tapArap$binomial == "bradypus variegatus", "binomial"] <- "Bradypus variegatus"
tapArap[tapArap$binomial == "Mazama americana", "binomial"] <- "Mazama sp."
tapArap[tapArap$binomial == "Mazama nemorivaga", "binomial"] <- "Mazama sp."
levels(tapArap$binomial) <- c(levels(tapArap$binomial),"Leopardus sp.")
tapArap[tapArap$binomial == "Leopardus pardalis", "binomial"] <- "Leopardus sp."
tapArap[tapArap$binomial == "Leopardus tigrinus", "binomial"] <- "Leopardus sp."
tapArap[tapArap$binomial == "Leopardus wiedii", "binomial"] <- "Leopardus sp."
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Cracidae',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Tinamidae',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'na',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Tayassuidae',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Felidae',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Primates',])
tapArap <- droplevels(tapArap[!tapArap$binomial == 'Sciuridae',])
tapArap$binomial <- factor(tapArap$binomial)
sort(unique(tapArap$binomial))

# somente mamiferos
tapArap <- tapArap %>%
  filter(classe == "Mammalia")


# tabela 1 com esforco por ano
#tabela1 <- tibble(
#  'Ano' = sort(unique(tapArap$ano)),
#  'Transectos' = 0,
#  'Esforço (km)' = 0,
#  'Registros' = 0
#)
#for(i in 1:nrow(tabela1)) {
#  df1 <- subset(tapArap, ano == unique(tapArap$ano)[i])
#  tabela1[i, "Transectos"] <- length(unique(df1$estacao.amostral))
#  tabela1[i, "Esforço (km)"] <- sum(df1$esforco, na.rm=T)/1000
#  tabela1[i, "Registros"] <- nrow(df1)
#}
#tabela1


## construir tabela 1 passo a passo -----

# esforco total
esforco_total <- tapArap %>%
  distinct(estacao.amostral, ano, data, esforco) %>%
  summarize(sum(esforco, na.rm = TRUE)/(1000))
esforco_total

# quantos transectos ativos a cada ano
transectos <- tapArap %>%
  group_by(ano) %>%
  distinct(estacao.amostral) %>%
  count() %>%
  rename(transectos = n)
transectos

# quanto esforco a cada ano em km
esforco <- tapArap %>%
  distinct(estacao.amostral, ano, data, esforco) %>%
  group_by(ano) %>%
  summarize(esforco = sum(esforco, na.rm = TRUE)/1000)
esforco

# quantos registros a cada ano
registros <- tapArap %>%
  group_by(ano) %>%
  summarise(registros = n())
registros

tabela1 <- left_join(transectos, esforco, by="ano") %>%
  left_join(registros, by="ano")
tabela1

# tabela com especies registradas
# tabela para o anexo
tabela2 <- tapArap %>%
  group_by(binomial) %>% 
  summarize(n = n(), ) %>% 
  arrange(desc(n)) %>%
  mutate(grupos_10km = round(n/pull(esforco_total/10), 2))
tabela2

print(tabela2, n=Inf)



#---------- Parte 1: calcular taxas de encontro ----------

# criar objeto para receber taxas de avistamento anuais
mydata <- distinct(tapArap, binomial) %>%
  mutate("2014" = as.numeric(NA),
         "2015" = as.numeric(NA),
         "2016" = as.numeric(NA),
         "2017" = as.numeric(NA),
         "2018" = as.numeric(NA),
         "2019" = as.numeric(NA)) #%>%
#  mutate(id = row_number()) %>%
#  relocate(id, .before = cduc)
mydata


# preencher objeto acima com numero de encontros por ano
vetor_ano <- seq(min(tapArap$ano), max(tapArap$ano))
for(i in 1:nrow(mydata)) {
  for(j in 1:length(vetor_ano)){
    mydata[i, 1+j] <- tapArap %>%
      filter(binomial == pull(mydata[i,"binomial"])) %>%
      filter(ano == vetor_ano[j]) %>%
      count() %>%
      pull(n)
  }}
mydata
#print(mydata, n=100)


# colocar esforco no formato wide
esforco <- esforco %>%
  pivot_wider(names_from = ano, values_from = esforco) %>%
  select("2014", "2015", "2016", "2017", "2018", "2019")
esforco


# corrigir pelo esforco e substituir zeros por NAs quando for o caso
for(i in 1:nrow(mydata)) {
  for(j in 1:length(vetor_ano)){
    mydata[i, j+1] <- mydata[i, j+1]/(esforco[j]/10)
  }}
mydata
#print(mydata, n=100)


# o esforco em 2014 foi muito baixo, excluir esse ano
mydata <- mydata %>%
  select(-'2014')

# se taxa avistamento for = 0, alterar para 0.001
mydata[mydata == 0] <- 0.001

# selecionar somente espécies que atendem a criterios minimos
# (taxa de avistamento >= 0.5) e pelo menos 3 anos de dados
mydata$mean <- rowMeans(mydata[,c(2:6)], na.rm = TRUE)
print(mydata, n=Inf)

mydata <- subset(mydata, mean >= 0.25)

mydata <- mydata %>%
  mutate(binomial = factor(binomial)) %>%
  select(-mean)

mydata

#---------- Parte 2: modelo bayesiano ----------

# check species names for analysis
binomial <- mydata$binomial
binomial

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

encounter_rate <- mydata

for(i in 1:nrow(mydata)) {
  y  <- mydata[i,] %>%
    select(-binomial) %>%
    as.numeric()
  state.space.model(y, n.years)
  # plot trends
  jpeg(here("report", paste(gsub(" ", "_", binomial)[i], "tapArap.jpg", sep="_")), width=1000, height=600, res=120) # Open jpeg file
  pop.trends()
  dev.off()
  # save ssm results to be used for making figures later independent of loop and add them to list
  N.est_all_species[[paste(gsub(" ", "_", binomial)[i], "ssm", sep="_")]] <- ssm$BUGSoutput$sims.list$N.est
  # build table 1
  tabela3[i,] <- list(binomial[i],
                      round(mean(meanR), 2),
                      paste(round(quantile(meanR, probs = c(0.025)), 2), " a ", round(quantile(meanR, probs = c(0.975)), 2), sep=""),
                      paste((round(length(which(meanR < 0))/length(meanR), 2)*100), "%", sep= ""),
                      paste((round(length(which(meanR > 0))/length(meanR), 2)*100), "%", sep= ""))
}
tabela3
#write.csv(tabela3, here("experimental", "tabela3_tapArap.csv"))
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
#print(xtable(tabela3), type="html", file=here("experimental", "tabela3_tapArap.html"))
#print(xtable(tabela3, type="html")

