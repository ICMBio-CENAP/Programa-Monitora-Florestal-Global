
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
#dados <- read.csv(here("data", "Planilha consolidada mastoaves até 2018  julho 2019 - FINAL.csv"), sep=",") # mais atual, mas com problemas
dados <- read.csv(here("data", "Planilha consolidada mastoaves até 2019 - FINAL.csv"))

# alguns ajustes
dados2 <- dados
dados2$Binomial <- gsub( " .*$", "", dados2$Espécies.validadas.para.análise.do.ICMBio) # trabalhar somente com genero para minimizar efeito de erros
dados2$Binomial <- paste(dados2$Binomial, dados2$CDUC, sep="_") # separa populacoes diferentes da mesma spp (neste caso genero)
dados2$Binomial <- as.factor(dados2$Binomial) # pode ser desnecessário, verificar
dados2 <- dados2[dados2$Ano != "",]
dados2$Ano <- str_sub(dados2$Ano, start= -4) # porque o ano está como factor/data
dados2 <- dados2[complete.cases(dados2$Ano), ] # remover linhas com NA na coluna Anos
dados2$Ano <- as.numeric(dados2$Ano)
#dados2 <- dados[dados$Ano != is.na,] # remover linhas com ano em branco
dados2 <- subset(dados2, Ano != 2020)

# Selecionando registros somente das espécies de interesse para o LPI
# pattern matching
# completar a lista abaixo com todas as espécies a serem incluídas
#toMatch <- c("Cebus", "Dasypr", "Allou", "Crax", "Psophia", "Penelo") # a list with at least part of the names of species for inclusion in LPI analysis
#toMatch <- c("Sapaj", "Dasypr", "Penelo") # versão com apenas três gêneros

#row.list <- grep(paste(toMatch,collapse="|"), dados2$Binomial , fixed=F) # que linhas incluir?
#dados3 <- dados2[c(row.list),] # subset a partir do passo anterior
#View(dados3)
colnames(dados2)[2] <- "nome.UC"
colnames(dados2)[3] <- "estacao.amostral"
colnames(dados2)[5] <- "esforço"
colnames(dados2)[9] <- "hora.inicio"
colnames(dados2)[10] <- "hora.fim"
colnames(dados2)[17] <- "Genero"
colnames(dados2)[18] <- "Especie"
colnames(dados2)[23] <- "hora.registro"
colnames(dados2)[24] <- "n.individuos"
colnames(dados2)[26] <- "dist.perpendicular"

# simplificar dados icmbio
to.remove <- c("Nome.da.EA", "Estação.do.ano", "condição.climática..aberto..nublado.e.chuvoso", 
               "nome.dos.observadores", "Número.do.animal.no.guia", "Clasificação.taxonômica..espécie..gênero..família.ou.ordem.",
                         "Espécies.validadas.para.análise.do.ICMBio", "Número.do.animal.no.guia.validado",
                         "Clasificação.taxonômica.validada", "teve.problema.durante.a.amostragem.",
                         "Observações", "X", "X.1", "X.2", "X.3")
`%ni%` <- Negate(`%in%`)
dados2 <- subset(dados2, select = names(dados2) %ni% to.remove)

# salvar como rds para proxima etapa
saveRDS(dados2, here("data", "dadosICMBio_2014a2019.rds"))
