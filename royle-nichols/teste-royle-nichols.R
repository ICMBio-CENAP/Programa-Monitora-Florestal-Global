# Modelos populacionais Bayesianos
# Baseado no livro Bayesian Population Analysis Using WinBUGS (Kery e Schaub)
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP


# Carregar pacotes
library(here)
library(tidyr)
library(stringr)
library(htmlTable)
library(ggplot2)
library(dplyr)
library(R2jags)

# carregar funcoes
#source(here("bin", "lpi_icmbio.R"))
source(here("royle-nichols", "population-functions-royle-nichols.R"))

# para os testes usar dados de Cazumba-Iracema
cazumba <- read.csv(here("experimental", "cazumba_2014-2020.csv"))
head(cazumba)

# adicionar ano e corrigir formato da data
cazumba$ano <- as.numeric(str_sub(cazumba$data.da.amostragem,-4,-1))
head(cazumba)

cazumba$data.da.amostragem <- as.Date(cazumba$data.da.amostragem, format = "%d/%m/%Y")

# adicionar esforco de 5km por sitio-dia
cazumba$site.day <- paste(cazumba$estacao.amostral, cazumba$data.da.amostragem, sep="_")
for(i in 1:length(unique(cazumba$site.day))) {
  cazumba[cazumba$site.day == unique(cazumba$site.day)[i],][1,]$esforco <- 5000
}
cazumba$esforço <- cazumba$esforco
cazumba # check

# como esses dados ainda nao foram validados, ha spp com nomes errados etc
# portanto usar somente genero nesse teste
cazumba$Genero <- gsub( " .*$", "", cazumba$Especie)
head(cazumba)

# criar uma coluna especie-trilha
cazumba$spp_trilha <- paste(cazumba$Genero, cazumba$estacao.amostral)

# criar coluna trilha_ano
cazumba$trilha_ano <- paste(cazumba$estacao.amostral, cazumba$ano)

#---------- Parte 1: Historico de deteccao ----------

# a funcao detection.history gera historico de deteccao para
# cada trilha e para cada dia de amostragem
# pode ser util se modelo tiver covariaveis de ocasiao
# no entanto por enquanto o modelo nao tera ocasiao, entao vamos desligar:
#teste <- detection.history(cazumba)
#teste["Dasyprocta"] # checar uma especie


# a funcao make.binomial.table cria df com dados para modelo binomial
# para cada visita (dia) a especie pode ser detectada (sucesso) ou nao (fracasso)
# pode ser util, mas por enquanto o formato nao e ideal para indexar no modelo do jags
# desligar por enquanto
#make.binomial.table(cazumba)
#binomial.table


# a funcao prepare.data.RNmodel gera dataframe y com numero de sucessos por estacao por ano
prepare.data.RNmodel(cazumba, "Dasyprocta")
y
prepare.data.RNmodel(cazumba, "Mazama") # lembrando que R-N nao e bom para especies gregarias
y

# a funcao acima contabiliza o numero de sucessos, mas nao o de tentativas
# para Cazumba isso nao e problema ja que as tentativas (esforco em numero de visitas) foram sempre 10
# mas para outros sitios onde o esforco foi variavel sera preciso calcular tentativas
# alem disso o modelo do jags abaixo vai precisar de index para "v": v[i,j]
# idealmente, adaptar funcao prepare.data.RNmodel para que gere um dataframe adicional com tentativas por sitio por ano
# neste teste vamos manter o numero de tentativas fixo = 10


## Modelo Role-Nichols
# Escrever modelo no formato do jags
sink(here("royle-nichols", "model.txt"))
cat("
	model{
	
	# Priors
	for(i in 1:nsite){
	  for(j in 1:nyear){
	    lambda[i,j] ~ dlnorm(0, 0.001)
	    r[i,j] ~ dunif(0,1)
	  } # j
	} # i

	# Likelihood
	for(i in 1:nsite){
	  for(j in 1:nyear){
	    N[i,j] ~ dpois(lambda[i,j])
	    p[i,j] <- 1-pow(1-r[i,j], N[i,j])
	    y[i,j] ~ dbin(p[i,j], v) # por enquanto fixando em 10, depois mudar para esforco variavel
		} #j
	} #i
	
	# Quantidades derivadas 
	for (j in 1:nyear){
   totalN[j] <- sum(N[,j])   # Estimate total pop. size across all sites
   mean.abundance[j] <- mean(lambda[,j])
   mean.N[j] <- mean(N[,j])
   mean.detection[j] <- mean(p[,j])
   } # j
}
",fill=TRUE)
sink()

y <- data.frame(y[,-1]) # remover primeira coluna com nome das estacoes
nsite <- nrow(y)
nyear <- ncol(y)
data <- list(y=y, nsite=nrow(y), nyear=ncol(y), v=10)
inits <- function()list(N=as.matrix(y),
                        lambda=matrix(rlnorm(nsite*nyear,0.5), nrow = nsite, ncol = nyear),
                        #lambda=as.matrix(y),
                        r=matrix(runif(nsite*nyear,0.2), nrow = nsite, ncol = nyear) )
parameters <- c("lambda","r", "totalN", "mean.abundance", "mean.N", "mean.detection")

# rodar o modelo
# obs: parece que o modelo e bem sensivel ao valor inicial de lambda
# as vezes precisa ajustar caso a caso dependendo da spp, verificar isso depois
out <- jags(data, inits, parameters, here("royle-nichols", "model.txt"),
            n.chain=3, n.burnin=10000, n.iter=20000)
out

# Summarize posteriors
print(out, dig = 2)

# usar funcao pop.trends para criar grafico
pop.trends(out, y)

