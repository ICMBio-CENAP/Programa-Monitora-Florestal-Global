
# Script para cálculo do LPI
# Baseado nas instruções disponíveis em https://github.com/Zoological-Society-of-London/rlpi
# Adaptado por Elildo Carvalho Jr, ICMBio/CENAP
# Função prepara dados do icmbio e chama funções do pacote rlpi
# Em fase de TESTE

# carregar pacotes
library(here)
library(tidyverse)
library(dplyr)
library(stringr)

# abrir planilha de dados
dados <- read_csv(here("data", "Planilha consolidada mastoaves até 2019 - FINAL.csv"))
dados
names(dados)

# alguns ajustes
dados2 <- dados %>%
  rename(cduc = "CDUC",
         nome.UC = "Local - Nome da Unidade de Conservação",
         estacao.amostral = "Número da Estação Amostral",
         esforço = "Esforço de amostragem tamanho da trilha (m)",
         ano = "Ano",
         hora.inicio = "horário de início  (h:mm)",
         hora.fim = "horário de término (h:mm)",
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
  filter(ano != "" | ! is.na(ano) | ano != 2020) %>%
  select(cduc, nome.UC, estacao.amostral, esforço, ano, hora.inicio, hora.fim,
         genero, binomial, populacao, nivel.validacao, hora.registro, n.individuos, dist.perpendicular)
dados2  


# Selecionr registros somente das espécies de interesse para o LPI
# pattern matching
# completar a lista abaixo com todas as espécies a serem incluídas
#toMatch <- c("Cebus", "Dasypr", "Allou", "Crax", "Psophia", "Penelo") # a list with at least part of the names of species for inclusion in LPI analysis
#toMatch <- c("Sapaj", "Dasypr", "Penelo") # versão com apenas três gêneros

#row.list <- grep(paste(toMatch,collapse="|"), dados2$Binomial , fixed=F) # que linhas incluir?
#dados3 <- dados2[c(row.list),] # subset a partir do passo anterior
# salvar como rds para proxima etapa
saveRDS(dados2, here("data", "dadosICMBio_2014a2019.rds"))
