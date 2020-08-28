
# Carregar pacotes
library(here)

# carregar funcoes
source(here("bin", "lpi_icmbio4.R"))

# carregar dados
dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2018.rds"))

# rodar funcao rlpi
lpi_icmbio(dadosICMBio) # cálculo do LPI para todo o conjunto de dados do ICMBio
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