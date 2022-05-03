
# Script para cálculo do LPI
# Baseado nas instruções disponíveis em https://github.com/Zoological-Society-of-London/rlpi
# Adaptado por Elildo Carvalho Jr, ICMBio/CENAP
# Função prepara dados do icmbio e chama funções do pacote rlpi
# Em fase de TESTE

# carregar pacotes
library(here)
library(tidyverse)
library(lubridate)
library(hms)
#library(dplyr)
#library(stringr)

#--------------------------------
# ler dados, versao nao disponibilizada, mas mais atualizada
dados <- read_csv(here("data", "Planilha consolidada mastoaves até 2019 - FINAL.csv"))
dados
names(dados)

# ajustar dados: versao nao disponibilizada, mas mais atualizada
dados <- dados %>%
  rename(cnuc = "CDUC",
         nome_UC = "Local - Nome da Unidade de Conservação",
         estacao_amostral = "Número da Estação Amostral",
         nome_ea = "Nome da EA",
         esforco = "Esforço de amostragem tamanho da trilha (m)",
         data = "data da amostragem",
         hora_inicio = "horário de início  (h:mm)",
         hora_fim = "horário de término (h:mm)",
         #ano = "Ano",
         classe = "Classe",
         ordem = "Ordem",
         familia = "Família",
         genero = "Gênero",
         binomial = "Espécies validadas para análise do ICMBio",
         n_animais = "n° de animais",
         distancia = "distância (m)     do animal em relação a trilha") %>%
  mutate(data = as.Date(data, "%d/%m/%Y"),
         ano = year(data),
         distancia =  as.numeric(str_replace(distancia, ",", ".")),
         populacao = paste(binomial, cnuc, sep ="_"),
         populacao = str_replace(populacao, " ", "_") ) %>%
  select(cnuc, nome_UC, estacao_amostral, nome_ea, esforco, ano, data,
         hora_inicio, hora_fim, classe, ordem, familia, genero, binomial,
         n_animais, distancia, populacao)
dados

# salvar como rds para proxima etapa
saveRDS(dados, here("data", "dadosICMBio_2014a2019.rds"))


#--------------------------------
# ler dados versao oficial disponibilizada
dados <- read_csv(here("data", "Dados_Florestal_14a18_disponibilizacao.csv"))

# ajustar dados: versao oficial disponibilizada
dados <- dados %>%
  rename(cnuc = "Cadastro Nacional de Unidades de Conservação (CNUC)",
         nome_UC = "Unidade de Conservação (UC)",
         estacao_amostral = "Número da Estação Amostral",
         nome_ea = "Nome da Estação Amostral",
         esforco = "Esforço de amostragem (metros percorridos por dia)",
         data = "data da amostragem (dd/mm/aaaa)",
         hora_inicio = "Horário de início  (hh:mm)",
         hora_fim = "Horário de término (hh:mm)",
         ano = "Ano",
         classe = "Classe",
         ordem = "Ordem",
         familia = "Família",
         genero = "Gênero",
         binomial = "Espécies validadas pelo ICMBio",
         n_animais = "N° de animais",
         distancia = "Distância perpendicular  (m) do animal em relação a trilha") %>%
  mutate(data = as.Date(data, "%d/%m/%Y"),
         distancia =  as.numeric(str_replace(distancia, ",", ".")),
         populacao = paste(binomial, cnuc, sep ="_"),
         populacao = str_replace(populacao, " ", "_") ) %>%
  select(cnuc, nome_UC, estacao_amostral, nome_ea, esforco, ano, data,
         hora_inicio, hora_fim, classe, ordem, familia, genero, binomial,
         n_animais, distancia, populacao)
dados


# salvar como rds para proxima etapa
saveRDS(dados, here("data", "dadosICMBio_2014a2018.rds"))
