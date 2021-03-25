
# Script para cálculo do LPI
# Baseado nas instruções disponíveis em https://github.com/Zoological-Society-of-London/rlpi
# Adaptado por Elildo Carvalho Jr, ICMBio/CENAP
# Função prepara dados do icmbio e chama funções do pacote rlpi
# Em fase de TESTE

# pegando diretório do script executado (Necessário usar ctrl/command-shift-s)
#setwd(dirname(parent.frame(2)$ofile))
#setwd("/media/elildojr/Dados/r/monitora")
library(here)

library(dplyr)
library(stringr)

# abrir planilha de dados mastoaves
dados <- read.csv("Planilha consolidada mastoaves até 2018  julho 2019 - FINAL.csv", sep=",") # mais atual, mas com problemas

dados2 <- dados[dados$Espécies.validadas.para.análise.do.ICMBio != "",] # remover linhas com nome de espécie em branco
dados2$Binomial <- dados2$Espécies.validadas.para.análise.do.ICMBio
dados2$Binomial <- gsub(" ", "_", dados2$Binomial) # necessário para rodar LPIMain
dados2$Binomial <- paste(dados2$Binomial, dados2$CDUC, sep="_") # NOVO em 20190405, para separar populacoes diferentes da mesma spp
dados2$Binomial <- as.factor(dados2$Binomial) # pode ser desnecessário, verificar
dados2 <- dados2[dados2$Ano != "",]
dados2$Ano <- str_sub(dados2$Ano, start= -4) # porque o ano está como factor/data
dados2 <- dados2[complete.cases(dados2$Ano), ] # remover linhas com NA na coluna Anos
dados2$Ano <- as.numeric(dados2$Ano)
#dados2 <- dados[dados$Ano != is.na,] # remover linhas com ano em branco
#dados2$Ano <- as.numeric(substr(dados2$Ano, 7, 10)) extrai ano quatro digitos
colnames(dados2)[5] <- "esforço"


# Selecionando registros somente das espécies de interesse para o LPI
# pattern matching
# completar a lista abaixo com todas as espécies a serem incluídas
#toMatch <- c("Cebus", "Dasypr", "Allou", "Crax", "Psophia", "Penelo") # a list with at least part of the names of species for inclusion in LPI analysis
#toMatch <- c("Sapaj", "Dasypr", "Penelo") # versão com apenas três gêneros

#row.list <- grep(paste(toMatch,collapse="|"), dados2$Binomial , fixed=F) # que linhas incluir?
#dados3 <- dados2[c(row.list),] # subset a partir do passo anterior
#View(dados3)

## carregar função "lpi_icmbio.R"
# a função agora separa as populações da mesma espécie de acordo com o código CNUC, ou seja, elas entram separadamente no cálculo
source("lpi_icmbio4.R")

# Testando a função

lpi_icmbio(dados2) # cálculo do LPI para todo o conjunto de dados do ICMBio
mydata2

  # selecionando populações para cálculo (taxa de avistamento >= 0.1)
  mydata2$mean <- rowMeans(mydata2[,c(3:7)], na.rm = TRUE)
  use.this <- subset(mydata2, mean >= 0.1) # usar somente espécies com taxa de avistamento médio > 0.1
  use.this <- use.this[rowSums(is.na(use.this)) <= 2, ]
  vector.Binomial <- as.character(factor(use.this$Binomial)) # criar vetor com "Binomials" (populações) a serem incluídas na análise
  dados4 <- filter(dados2, Binomial %in% vector.Binomial) # keeping only species in species.list
  dados4$Binomial <- factor(dados4$Binomial) # to remove excluded species from factor levels otherwise they will end up as zeros in the paMatrix

# excluir especies noturnas
# excluir especies com baixa detectabilidade/raras com base na literatura previa
# e.g felideos, preguiça etc

  
lpi_icmbio(dados4)
mydata2

# quantas spp de mamiferos e aves (cuttoff 0.1)
aves <- subset(dados4, Classe == "Aves")
  sort(unique(aves$Espécies.validadas.para.análise.do.ICMBio))
mamiferos <- subset(dados4, Classe == "Mamíferos")
  sort(unique(mamiferos$Espécies.validadas.para.análise.do.ICMBio))

lpi_icmbio(dados4,z="Mamíferos") # cálculo somente para mamíferos
lpi_icmbio(dados4,z="Aves") # somente para aves

# quantas UCs
sort(unique(dados4$Local...Nome.da.Unidade.de.Conservação))


lpi_icmbio(dados2, "Resex Cazumbá-Iracema", "Aves")
mydata2
lpi_icmbio(dados2, "Esec Terra do Meio")
mydata2
mydata
dim(mydata)
View(mydata)
class(mydata$Ano)

lpi_icmbio(dados2, "Resex Tapajós-Arapiuns", "Mamíferos")
lpi_icmbio(dados3, "Resex Tapajós-Arapiuns", "Aves")


lpi_icmbio(dados3, "Resex Tapajós-Arapiuns")


lpi_icmbio(dados3, "Esec Terra do Meio", "Mamíferos")
lpi_icmbio(dados3, "Esec Terra do Meio", "Aves")
lpi_icmbio(dados3, "Esec Terra do Meio")

lpi_icmbio(dados2, "Flona do Jamari", "Mamíferos")
lpi_icmbio(dados3, "Flona do Jamari", "Aves")
lpi_icmbio(dados3, "Flona do Jamari")

lpi_icmbio(dados2, "Parna do Juruena", "Mamíferos")
lpi_icmbio(dados2, "Parna do Juruena", "Aves")
lpi_icmbio(dados2, "Parna do Juruena")

lpi_icmbio(dados2, "Rebio Guaribas", "Mamíferos")
lpi_icmbio(dados2, "Rebio Guaribas", "Aves")
lpi_icmbio(dados2, "Rebio Guaribas")