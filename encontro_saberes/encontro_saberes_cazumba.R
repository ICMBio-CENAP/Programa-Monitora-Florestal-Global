# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyverse)
library(lubridate)
#library(tidyr)
#library(htmlTable)
library(ggplot2)
#library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
#source(here("experimental", "population-functions_cazumba2014-2020.R"))

# para os testes usar dados de Cazumba-Iracema
cazumba <- read_csv(here("data", "cazumba-2014-2021.csv")) 

# arrumar os dados
cazumba <- cazumba %>%
  rename(nome_UC = "Local - Nome da Unidade de Conservação",
         estacao_amostral = "Número da Estação Amostral",
         nome_ea = `Nome da EA`,
         esforco = "Esforço de amostragem tamanho da trilha (m)",
         data = "data da amostragem",
         ano = "Ano",
         classe = "Classe",
         ordem = "Ordem",
         familia = "Família",
         genero = "Gênero",
         especie = "Espécies validadas para análise do ICMBio",
         n_animais = "n° de animais",
         distancia = "distância (m)     do animal em relação a trilha") %>%
  select(nome_UC, estacao_amostral, nome_ea, esforco, data, ano, classe, ordem,
         familia, genero, especie, n_animais, distancia) %>%
  mutate(data = as.Date(data, "%d/%m/%Y"),
         distancia =  as.numeric(str_replace(distancia, ",", ".")) )
cazumba

# agrupar sciurideos
cazumba <- cazumba %>%
  mutate(especie = replace(especie, especie %in% c("Guerlinguetus ignitus", "Urosciurus spadiceus", "Microsciurus flaviventer"),
                           "Sciuridae"))

# tabela 1: registros por especie
tabela_1 <- cazumba %>%
  filter(classe == "Mammalia",
         ordem != "Primates") %>%
  mutate(especie = replace(especie, especie %in% c("Guerlinguetus ignitus", "Urosciurus spadiceus", "Microsciurus flaviventer"),
                           "Sciuridae")) %>%
  group_by(especie) %>%
  summarize(n = n(),
            n_ind = sum(n_animais))%>%
  left_join(cazumba %>%
              select(ordem, familia, especie), by="especie") %>%
  distinct(especie, .keep_all = TRUE) %>%
  arrange(ordem, familia, especie) %>%
  mutate(nome_popular = ifelse(especie == "Mazama americana", "Veado-mateiro",
                        ifelse(especie == "Mazama americana", "Veado-mateiro",
                        ifelse(especie == "Pecari tajacu", "Caititu",
                        ifelse(especie == "Tayassu pecari", "Queixada",
                        ifelse(especie == "Speothos venaticus", "Janauira",
                        ifelse(especie == "Leopardus pardalis", "Maracajá-açu",
                        ifelse(especie == "Leopardus tigrinus", "Gato-do-mato",
                        ifelse(especie == "Leopardus wiedii", "Maracajá-peludo",
                        ifelse(especie == "Panthera onca", "Onça-pintada",
                        ifelse(especie == "Puma concolor", "Onça-vermelha",
                        ifelse(especie == "Puma yagouaroundi", "Gato-mourisco",
                        ifelse(especie == "Eira barbara", "Irara",
                        ifelse(especie == "Galictis vittata", "Furão",
                        ifelse(especie == "Nasua nasua", "Quati",
                        ifelse(especie == "Cabassous unicinctus", "Tatu-de-rabo-mole",
                        ifelse(especie == "Dasypus novencimctus", "Tatu-galinha",
                        ifelse(especie == "Priodontes maximus", "Tatu-canastra",
                        ifelse(especie == "Tapirus terrestris", "Anta",
                        ifelse(especie == "Bradypus variegatus", "Preguiça",
                        ifelse(especie == "Choloepus hoffmanni", "Preguiça-real",
                        ifelse(especie == "Myrmecophaga tridactyla", "Tamanduá-bandeira",
                        ifelse(especie == "Tamandua tetradactyla", "Mambira",
                        ifelse(especie == "Dasyprocta fuliginosa", "Cutia",
                        ifelse(especie == "Myoprocta pratti", "Cutiara",
                        ifelse(especie == "Coendou bicolor", "Quandú",
                        ifelse(especie == "Sciuridae", "Quatipuru",
                        "no" ))))))))))))))))))))))))))) %>%
  select(ordem, familia, especie, nome_popular, n, n_ind) %>%
  rename(Ordem = ordem, Familia = familia, Especie = especie, 
         'Nome popular' = nome_popular, Encontros = n,
         Individuos = n_ind)
tabela_1

total_encontros_por_spp <- tabela_1 %>%
  filter(Encontros > 5) %>%
  ggplot(aes(x = reorder(`Nome popular`, -Encontros), y = Encontros)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Espécie") + ylab("Número de encontros") +
  theme_bw() +
  theme(text = element_text(size=14))
total_encontros_por_spp

# salvar jpeg
ggsave(here("encontro_saberes", "total_encontros_por_spp.jpeg"), 
       width = 10, height = 6, dpi = 150)



# esforco por ano
esforco_anual <- cazumba %>%
  group_by(ano) %>%
  summarize(esforco = sum(esforco, na.rm = TRUE))

# numero de registros por ano
n_registros <- cazumba %>%
  filter(classe == "Mammalia",
         ordem != "Primates") %>%
  group_by(ano, especie) %>%
  count()

# taxa de encontro por ano
taxa_encontro <- left_join(esforco_anual, n_registros, by=c("ano")) %>%
  mutate(taxa_encontro = n/(esforco/1000)) %>%
  select(ano, especie, taxa_encontro) %>%
  arrange(especie, ano)

# colocar em formato wide
taxa_encontro_wide <- taxa_encontro %>%
  pivot_wider(names_from = ano, values_from = taxa_encontro) %>%
  select(especie, `2014`, `2015`, `2016`, `2017`,  `2018`, `2019`, `2020`, `2021`) %>%
  replace(is.na(.), 0.001)

# selecionar as especies mais comuns (taxa media > 0.025)
taxa_encontro_wide <- taxa_encontro_wide %>%
  mutate(means = rowMeans(.[,2:9])) %>%
  arrange(desc(means)) %>%
  filter(means > 0.05) %>%
  #select(-means) %>%
  mutate_at(2:9, round, 2)


# funcao para gerar y por trilha para sp selecionada
calcular.registros.por.trilha <- function(nome_especie) {
  
  # esforco por ano por trilha
  esforco_anual_trilha <- cazumba %>%
    group_by(ano, estacao_amostral) %>%
    summarize(esforco = sum(esforco, na.rm = TRUE)) %>%
    pivot_wider(names_from = ano, values_from = esforco) %>%
    select(estacao_amostral, `2014`, `2015`, `2016`, `2017`,  `2018`, `2019`, `2020`, `2021`)
  esforco_anual_trilha
  
  # numero de registros por ano por trilha
  n_registros_trilha <- cazumba %>%
    filter(especie == nome_especie) %>%
    group_by(ano, estacao_amostral) %>%
    count() %>%
    mutate(n = as.numeric(n)) %>%
    pivot_wider(names_from = ano, values_from = n) %>%
    select(estacao_amostral, `2014`, `2015`, `2016`, `2017`,  `2018`, `2019`, `2020`, `2021`) %>%
    replace(is.na(.), 0.001)
  n_registros_trilha
  
  # combinar registros e esforco para obter taxa de encontro
  esforco_anual_trilha <- esforco_anual_trilha %>%
    select(-estacao_amostral)
  n_registros_trilha <- n_registros_trilha %>%
    ungroup() %>%
    select(-estacao_amostral)
  
  y <- n_registros_trilha/(esforco_anual_trilha/1000)
  
  return(y)
  
}


# funcao para rodar modelo populacional bayesiano
state.space.model <- function(y) {
  # especificar modelo na linguagem BUGS
  sink(here("encontro_saberes", "ssm.jags"))
  cat("
  model {
  # Priors and constraints
  for(ea in 1:3) {
  logN.est[ea, 1] ~ dnorm(0, 0.01)       # Prior for initial population size
  mean.r[ea] ~ dnorm(1, 0.001)             # Prior for mean growth rate
  }
  sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
  sigma2.proc <- pow(sigma.proc, 2)
  tau.proc <- pow(sigma.proc, -2)
  sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
  sigma2.obs <- pow(sigma.obs, 2)
  tau.obs <- pow(sigma.obs, -2)
  
  # Likelihood
  # State process
  for (year in 1:(n_years-1)) {
    for(ea in 1:3) {
      r[ea, year] ~ dnorm(mean.r[ea], tau.proc)
      logN.est[ea, year+1] <- logN.est[ea, year] + r[ea, year]
    }
   }
  # Observation process
    for (year in 1:n_years) {
      for(ea in 1:3) {
        y[ea, year] ~ dnorm(logN.est[ea, year], tau.obs)
      }
   }
   
  # Population sizes on real scale
  for (year in 1:n_years) {
    for(ea in 1:3) {
     N.est[ea, year] <- exp(logN.est[ea, year])/100
    }
  }
  
  # Derived parameter: mean encounter rate
  for(year in 1:n_years) {
      mean_N[year] <- mean(N.est[, year])
  }
    
  }
  ",fill = TRUE)
  sink()
  
  #  definir y
  y <-  y
  
  # definir numero de anos
  n_years <- ncol(y)
  
  # juntar os dados
  jags.data <- list(y = log(y*100), n_years = n_years)
  
  # Initial values
  #inits <- function(){
  #  list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1),
  #       sigma.obs = runif(1, 0, 1),
  #       LogN.est = rbind(c(rnorm(1, -0.5, 0.1), rep(NA, (n.years-1))),
  #                        c(rnorm(1, -0.5, 0.1), rep(NA, (n.years-1))),
  #                        c(rnorm(1, -0.5, 0.1), rep(NA, (n.years-1))) )
  #} 
  #inits()
  
  # Parametros monitorados
  parameters <- c("mean_N", "r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")
  
  # definicoes MCMC
  ni <- 25000
  nt <- 3
  nb <- 10000
  nc <- 3
  
  # chamar o JAGS a partir do R
  ssm <- jags(jags.data, inits=NULL, parameters, here("encontro_saberes", "ssm.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # checar resultados
  print(ssm, digits = 2)
  
  # Probabiliade de N(2021) < N(2014)
  #mean(ssm$BUGSoutput$sims.list$N.est[,8] < ssm$BUGSoutput$mean$N.est[1])
  
  assign("N.est", ssm$BUGSoutput$sims.list$N.est, .GlobalEnv)
  assign("mean_N", ssm$BUGSoutput$sims.list$mean_N, .GlobalEnv) 
  assign("meanR", ssm$BUGSoutput$sims.list$mean.r, .GlobalEnv) 
  assign("ssm", ssm, .GlobalEnv)
  
}


# funcao para rodar modelo e gerar grafico para cada especie 
computar.plotar.taxa_encontro <- function(nome_especie) {
  
  # preparar dados
  y <- calcular.registros.por.trilha(nome_especie)
  n_years <- ncol(y)
  
  # rodar modelo bayesiano
  state.space.model(y)
  
  # computar quantis
  modelo_quantis <- as_tibble(mean_N) %>%
    pivot_longer(1:8, names_to = "year") %>%
    group_by(year) %>%
    summarize(quant025 = quantile(value, probs = 0.025),
              quant25 = quantile(value, probs = 0.25), 
              quant50 = quantile(value, probs = 0.5),
              quant75 = quantile(value, probs = 0.75),
              quant975 = quantile(value, probs = 0.975)) %>%
    mutate(ano = c(2014:2021))
  modelo_quantis
  
  # plotar com ribbons
  plot_tendencia <- ggplot() + 
    geom_ribbon(data=modelo_quantis, aes(x=ano, ymin=quant025, ymax=quant975),fill="coral", alpha=0.2) +
    geom_ribbon(data=modelo_quantis, aes(x=ano, ymin=quant25, ymax=quant75),fill="coral3", alpha=0.3) +
    geom_line(data=modelo_quantis, aes(x=ano, y=quant50), size=0.5, alpha=0.5) +
    geom_point(data=modelo_quantis, aes(x=ano, y=quant50), size=2, alpha=0.5) +
    #stat_smooth(data=modelo_quantis, aes(x=ano, y=quant50), method = "lm", formula = y ~ poly(x, 5), se = FALSE) +
    ylim(0, max(modelo_quantis$quant975)+0.05) +
    xlab("Ano") + 
    ylab("Taxa de encontro") +
    theme_bw() +
    theme(text = element_text(size=14)) #+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  plot_tendencia
  
  # salvar jpeg
  ggsave(here("encontro_saberes", paste(nome_especie, "2014 a 2021.jpeg", sep = " ")) )
  
  
  # plotar por transecto
  dim(N.est)
  N_trilha1 <- apply(N.est[,1,], 2, mean)
  N_trilha2 <- apply(N.est[,2,], 2, mean)
  N_trilha3 <- apply(N.est[,3,], 2, mean)
  N.est_por_trilha <- as_tibble(rbind(N_trilha1,
                                       N_trilha2,
                                       N_trilha3)) %>%
    rename('2014' = V1, '2015' = V2, '2016' = V3, '2017' = V4, '2018' = V5,
           '2019' = V6, '2020' = V7, '2021' = V8) %>%
    mutate(trilha = c("Cazumbá", "Gamas", "Alto Caeté")) %>%
    pivot_longer(1:8, names_to = "ano")
  N.est_por_trilha
  plot_tendencia_trilha <- ggplot() + 
    geom_line(data=N.est_por_trilha, aes(x=ano, y=value, group=trilha, color=trilha), size=2, alpha=0.5) +
    geom_point(data=N.est_por_trilha, aes(x=ano, y=value, group=trilha, color=trilha), size=2, alpha=0.5) +
    xlab("Ano") + 
    ylab("Taxa de encontro") +
    theme_bw() +
    theme(text = element_text(size=14))
  plot_tendencia_trilha
  
  # salvar jpeg
  ggsave(here("encontro_saberes", paste(nome_especie, "2014 a 2021 por trilha.jpeg", sep = " ")) )
  
}

computar.plotar.taxa_encontro("Sciuridae")
computar.plotar.taxa_encontro("Dasyprocta fuliginosa")
computar.plotar.taxa_encontro("Pecari tajacu")
computar.plotar.taxa_encontro("Myoprocta pratti")



#--------------------------------------------------

# Densidade populacional

# tutorial Rdistance

install.packages("Rdistance")
library(Rdistance)

# ver formato dos dados de modelo
data(sparrowDetectionData)
head(sparrowDetectionData)
data(sparrowSiteData)
head(sparrowSiteData)

cutiaDetectionData <- cazumba %>%
  mutate(siteID = estacao_amostral,
         groupsize = n_animais,
         dist = distancia) %>%
  select(siteID, groupsize, dist) %>%
  filter(!is.na(dist))
cutiaDetectionData

cutiaSiteData <- cazumba %>%
  mutate(siteID = estacao_amostral) %>%
  group_by(siteID) %>%
  summarize(length = sum(esforco, na.rm = TRUE))
cutiaSiteData

# explorar distribuicao de distancias
hist(cutiaDetectionData$dist, col="grey", main="",
     xlab="distance (m)")
rug(cutiaDetectionData$dist, quiet = TRUE)

summary(cutiaDetectionData$dist)


# truncar distancias ate 30m
cutiaDetectionData <- cutiaDetectionData %>%
  mutate(dist = ifelse(dist > 25, 25, dist))

# plotar novamente  
hist(cutiaDetectionData$dist, col="grey", main="",
     xlab="distance (m)")
rug(cutiaDetectionData$dist, quiet = TRUE)

# ajustar uma funcao de deteccao
dfuncCutia<- dfuncEstim(formula = dist~1,
                          detectionData = cutiaDetectionData,
                          siteData = cutiaSiteData,
                          likelihood = "halfnorm", w.hi = 50)
plot(dfuncCutia)
dfuncCutia
# The effective strip width (ESW) is the key information from the detection function
# that will be used to next estimate abundance (or density).

# estimar abundancia a partir da funcao de deteccao
# usar area=10000 converte a estimativa para hectares
fit <- abundEstim(dfuncCutia,
                  detectionData=cutiaDetectionData,
                  siteData=cutiaSiteData,
                  area=10000,
                  ci=NULL)
# To save vignette build time, we insert values from
# a separate run with R=500
fit$ci <- c( lowCI=0.6473984, hiCI=1.0167007)
fit

# resultados de interesse, como a estimativa de abundancia e intervalo de confianca
# podem ser extraidos do objeto (aqui chamado fit).
fit$n.hat
fit$ci

# 2 ind/ha equivale a 200 ind/km2 (1km2 = 100ha)

#--------------------------
# colocar tudo numa funcao

estimar.densidade.monitora <- function(x) {
  library(Rdistance)
  
  DetectionData <- cazumba %>%
    filter(especie == x) %>%
    mutate(siteID = estacao_amostral,
           groupsize = n_animais,
           dist = distancia) %>%
    select(siteID, groupsize, dist) %>%
    filter(!is.na(dist))
  DetectionData
  
  SiteData <- cazumba %>%
    mutate(siteID = estacao_amostral) %>%
    group_by(siteID) %>%
    summarize(length = sum(esforco, na.rm = TRUE))
  SiteData
  
  # truncar distancias ate 25m
  DetectionData <- DetectionData %>%
    mutate(dist = ifelse(dist > 25, 25, dist))
  
  # explorar distribuicao de distancias
  hist(DetectionData$dist, col="grey", main="",
       xlab="distance (m)")
  rug(DetectionData$dist, quiet = TRUE)
  
  # ajustar uma funcao de deteccao
  dfunc<- dfuncEstim(formula = dist~1,
                          detectionData = DetectionData,
                          siteData = SiteData,
                          likelihood = "halfnorm", w.hi = 30)
  plot(dfunc)
  dfunc
  # The effective strip width (ESW) is the key information from the detection function
  # that will be used to next estimate abundance (or density).
  
  # estimar abundancia a partir da funcao de deteccao
  # usar area=10000 converte a estimativa para hectares
  # usar area=1000000 converte a estimativa para km2
  fit <- abundEstim(dfunc,
                    detectionData=DetectionData,
                    siteData=SiteData,
                    area=1000000,
                    ci=NULL)
  # To save vignette build time, we insert values from
  # a separate run with R=500
  fit$ci <- c( lowCI=0.6473984, hiCI=1.0167007)
  fit
  
  # resultados de interesse, como a estimativa de abundancia e intervalo de confianca
  # podem ser extraidos do objeto (aqui chamado fit)
  return(fit$n.hat)
  #return(list(fit$n.hat, fit$ci))
  #return(fit$ci)
  
}

estimar.densidade.monitora("Dasyprocta fuliginosa")

estimar.densidade.monitora("Myoprocta pratti")

estimar.densidade.monitora("Pecari tajacu")

# populacao total aproximada: multiplicar pela area de Cazumba ~7000 km2

# cutia
estimar.densidade.monitora("Dasyprocta fuliginosa")*7000

estimar.densidade.monitora("Sapajus macrocephalus")
estimar.densidade.monitora("Alouatta puruensis")
estimar.densidade.monitora("Cebus unicolor")


