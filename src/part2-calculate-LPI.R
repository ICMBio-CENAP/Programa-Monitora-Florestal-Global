# Script para cálculo do LPI
# Baseado nas instruções disponíveis em https://github.com/Zoological-Society-of-London/rlpi
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP
# Função prepara dados do icmbio e chama funções do pacote rlpi
# Em fase de TESTE
# NB! Funcao lpi_icmbio esta salvando arquivos na pasta principal, editar para que destino seja subpasta "results"


# Carregar pacotes
library(here)
library(tidyr)
library(ggplot2)
library(dplyr)

# carregar funcoes
source(here("bin", "lpi_icmbio.R"))

# carregar dados
dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2019.rds"))

# rodar funcao rlpi
lpi_icmbio(dadosICMBio) # cálculo do LPI para todo o conjunto de dados do ICMBio
#mydata2

# salvar o grafico
ggsave(file = here("results", "lpi-global2016-2019.jpg"), plot = ggplot_lpi_modif(mydata_lpi, col="darkcyan"))


  # selecionando populações para cálculo (taxa de avistamento >= 0.1) e pelo menos 
  mydata2$mean <- rowMeans(mydata2[,c(3:8)], na.rm = TRUE)
  use.this <- subset(mydata2, mean >= 0.1) # usar somente espécies com taxa de avistamento médio > 0.1
  use.this <- use.this[rowSums(is.na(use.this)) <= 3, ]
  vector.Binomial <- as.character(factor(use.this$Binomial)) # criar vetor com "Binomials" (populações) a serem incluídas na análise
  dados4 <- subset(dadosICMBio, Binomial %in% vector.Binomial) # keeping only species in species.list
  dados4$Binomial <- factor(dados4$Binomial) # to remove excluded species from factor levels otherwise they will end up as zeros in the paMatrix

# excluir especies noturnas
# excluir especies com baixa detectabilidade/raras com base na literatura previa
# e.g felideos, preguiça etc

lpi_icmbio(dados4)
ggsave(file = here("results", "lpi-spp-selecionadas-2016-2019.jpg"), plot = ggplot_lpi_modif(mydata_lpi, col="darkcyan"))


# tendencias de especies individuais (grafico espaguete)
mydata2 <- mydata2[3:nrow(mydata2),] # as duas primeiras linhas estavam com "falsas" especies
speciesTrends <- mydata2
speciesTrends <- speciesTrends %>% gather(Ano, value, X2014, X2015, X2016, X2017, X2018, X2019, na.rm = FALSE, convert = FALSE)
speciesTrends$Ano <- gsub("X", "", speciesTrends$Ano)
speciesTrends$value[speciesTrends$value == 0] <- NA # remover zeros apenas para deixar grafico bonito

# Plot
plotSpeciesTrends <- speciesTrends %>%
  ggplot( aes(x=Ano, y=log(value), group=Binomial, color=Binomial)) +
  geom_line(size=0.6, alpha=0.5) + theme_classic()
#plotSpeciesTrends
plotSpeciesTrends  + ggplot2::ylab("Taxa de avistamento (escala log)") + theme(legend.position = "none")


# quantas spp de mamiferos e aves (cuttoff 0.1)
allPops <- as.factor(mydata2$Binomial)
levels(allPops)
allSpp <- as.factor(gsub( "_.*$", "", mydata2$Binomial))
levels(allSpp)
aves <- subset(dados4, Classe == "Aves")
  sort(unique(aves$Binomial))
  avesSpp <- as.factor(gsub( "_.*$", "", aves$Binomial))
  levels(avesSpp)
mamiferos <- subset(dados4, Classe == "Mammalia")
  sort(unique(mamiferos$Binomial))
  mammaliaSpp <- as.factor(gsub( "_.*$", "", mamiferos$Binomial))
  levels(mammaliaSpp)
ucs <- as.factor(sort(unique(dados4$Local...Nome.da.Unidade.de.Conservação)))
  levels(ucs)
  

lpi_icmbio(dados4,z="Mammalia") # cálculo somente para mamíferos
ggsave(file = here("results", "lpi-mamiferos-selecionados-2016-2019.jpg"), plot = ggplot_lpi_modif(mydata_lpi, col="cornflowerblue"))

lpi_icmbio(dados4,z="Aves") # somente para aves
ggsave(file = here("results", "lpi-aves-selecionadas.jpg-2016-2019"), plot = ggplot_lpi_modif(mydata_lpi, col="cornflowerblue"))


# quantas UCs
sort(unique(dados4$Local...Nome.da.Unidade.de.Conservação))
aves$Local...Nome.da.Unidade.de.Conservação <- as.factor(aves$Local...Nome.da.Unidade.de.Conservação)
levels(aves$Local...Nome.da.Unidade.de.Conservação)
sort(unique(aves$Local...Nome.da.Unidade.de.Conservação))


lpi_icmbio(dados4, "Resex Cazumbá-Iracema", "Aves")
lpi_icmbio(dados4, "Resex Cazumbá-Iracema", "Mammalia")

lpi_icmbio(dados4, "Resex Cazumbá-Iracema")

lpi_icmbio(dados4, "Flona do Jamari", "Aves")
lpi_icmbio(dados4, "Flona do Jamari", "Mammalia")

lpi_icmbio(dados4, "Rebio do Uatumã", "Mammalia")

lpi_icmbio(dados4, "Resex Tapajós-Arapiuns", "Mammalia")
lpi_icmbio(dados4, "Parna Montanhas do Tumucumaque", "Mammalia")

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