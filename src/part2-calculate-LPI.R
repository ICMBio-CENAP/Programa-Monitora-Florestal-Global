# Script para cálculo do LPI
# Baseado nas instruções disponíveis em https://github.com/Zoological-Society-of-London/rlpi
# Adaptado por Elildo Carvalho Jr @ ICMBio/CENAP
# Função prepara dados do icmbio e chama funções do pacote rlpi
# Em fase de TESTE
# NB! Funcao lpi_icmbio esta salvando arquivos na pasta principal, editar para que destino seja subpasta "results"


# Carregar pacotes
library(here)
library(tidyverse)
library(rlpi)
#library(tidyr)
#library(ggplot2)
#library(dplyr)

# carregar funcoes
source(here("bin", "lpi_icmbio.R"))

# carregar dados
# versao oficial disponibilizada
#dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2018.rds"))
# versao nao disponibilizada, mas mais atualizada
dadosICMBio <- readRDS(here("data", "dadosICMBio_2014a2019.rds"))
mydata <- dadosICMBio

# quantas UCs
mydata %>%
  group_by(nome_UC) %>%
  count()

# quantas especies
mydata %>%
  group_by(binomial) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(! binomial  %in% c(NA, "E", "Tinamidae", "Cracidae", "Sciuridae", "Callitrichidae",
                          "Pitheciidae", "Atelidae", "Cebidae", "Primates", "Felidae",
                          "Bradypodidae", "Cervidae")) %>%
  filter(! grepl('sp.', binomial)) %>%
  print(n=Inf)

# quantas populacoes
mydata %>%
  filter(! binomial  %in% c(NA, "E", "Tinamidae", "Cracidae", "Sciuridae", "Callitrichidae",
                            "Pitheciidae", "Atelidae", "Cebidae", "Primates", "Felidae",
                            "Bradypodidae", "Cervidae")) %>%
  group_by(populacao) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(! grepl('sp.', populacao)) %>%
  filter() %>%
  print(n=Inf)

# quantos transectos
mydata %>%
  distinct(nome_UC, estacao_amostral) %>%
  count()

# esforco total
mydata %>%
  summarize(sum(esforco, na.rm = TRUE)/1000)


# calcular esforco anual por UC em km
esforco <- mydata %>%
  group_by(cnuc, ano) %>%
  summarize(esforco = sum(esforco, na.rm = TRUE)/1000) %>%
  filter(esforco > 149) 
esforco


# selecionar UCs com > 1 ano de amostragem e > 149 km de esforco anual
ucs_selecionadas <- esforco %>%
  group_by(cnuc) %>% 
  count() %>%
  filter(n > 1) %>%
  pull(cnuc)
ucs_selecionadas


# filtrar mydata mantendo somente UCs selecionadas
mydata <- mydata %>%
  filter(cnuc %in% ucs_selecionadas)
mydata


# remover especies com "sp." pois podem representar > 1 especie
mydata <- mydata %>%
  filter(! grepl('sp.', binomial))
mydata


# esforco por ano em km
esforco_anual <- mydata %>%
  group_by(cnuc, ano) %>%
  summarize(esforco = sum(esforco, na.rm = TRUE)) %>%
  mutate(esforco = esforco/1000)
esforco_anual

# formato wide
esforco_anual_wide <- esforco_anual %>%
  pivot_wider(names_from = ano, values_from = esforco) %>%
  select(cnuc, `2014`, `2015`, `2016`, `2017`,  `2018`, `2019`) %>%
  ungroup()
esforco_anual_wide

# numero de registros por ano
n_registros <- mydata %>%
  group_by(cnuc, ano, populacao) %>%
  count() %>%
  mutate(n = as.numeric(n))
n_registros

# formato wide
n_registros_wide <- n_registros %>%
  pivot_wider(names_from = ano, values_from = n) %>%
  select(cnuc, populacao, `2014`, `2015`, `2016`, `2017`,  `2018`, `2019`) %>%
  replace(is.na(.), 0) %>%
  ungroup()
n_registros_wide

# taxa de encontro anual
# dividir numero de registros por esforco, NAs devem ser mantidos!
taxas_anuais <- data.frame(n_registros_wide)
esforcos_anuais <- data.frame(esforco_anual_wide)
for(i in 1:nrow(taxas_anuais)) {
  cnuc_temp <- taxas_anuais[i, "cnuc"]
  esforco_temp <- esforcos_anuais %>% filter(cnuc == cnuc_temp) %>% select(-cnuc)
  taxas_anuais[i, 3:8] <- round(taxas_anuais[i, 3:8]/esforco_temp, 3)
}
taxas_anuais <- as_tibble(taxas_anuais)
taxas_anuais %>%
  print(n=Inf)


# selecionar populacoes com taxa de avistamento media acima de 0.1 (1 ind/10km)
taxas_anuais <- taxas_anuais %>%
  mutate(media = rowMeans(.[,3:8])) %>%
  filter(media > 0.1) %>%
  select(-media)
taxas_anuais
#print(taxas_anuais, n=Inf)


#---------------------------------------------
# Problemas para rodar LPI!!!
# linhas abaixo sao tentativa incompleta para resolver...

# alguns ajustes para rodar lpi
mydata_lpi <- taxas_anuais %>%
  mutate(ID = as.numeric(1:nrow(taxas_anuais))) %>%
  rename(Binomial = populacao) %>%
  select(cnuc, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
mydata_lpi

# criar objetos que necessarios para lpi
index_vector = rep(TRUE, nrow(mydata_lpi))
year_vector <- 2014:2019

# criar infile
mydata_infile <- create_infile(mydata_lpi, index_vector=index_vector, 
                               #name=here("lpi", "mydata_data"),
                               name=here("mydata_data"),
                               start_col_name = colnames(mydata_lpi)[4], end_col_name = tail(colnames(mydata_lpi), n=1),
                               CUT_OFF_YEAR = year_vector[1])

# Calcular LPI com 100 bootstraps
#source(here("bin", "LPIMain.R")) # adicionado por não ter pacote rlpi
mydata_lpi <- LPIMain(#infile = here("lpi", "mydata_data_infile.txt"),
                      infile = here("mydata_data_infile.txt"),
                      #basedir = here("lpi"),
                      REF_YEAR = year_vector[1],
                      PLOT_MAX = tail(year_vector, n=1),
                      BOOT_STRAP_SIZE = 100, VERBOSE=FALSE)

# ... parei aqui, nao resolvido ainda nas linhas acima!
#---------------------------------------------

# Remover NAs (anos seguidos sem dados)
mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]

# carregar função ggplot_lpi_modif
source(here("bin", "ggplot_lpi_modif.R"))

# Gerar gráfico mais bonito usando função ggplot_lpi_modif
ggplot_lpi_modif(mydata_lpi, col="darkcyan")

# salvar mydata_lpi na area de trabalho
assign("mydata_lpi", mydata_lpi, .GlobalEnv)


#---------------------------------------------

#---------------------------------------------



# verificar quais especies ficaram e selecionar nao-sociais
sort(unique(mydata4$binomial))
nao_sociais <- c("Dasyprocta croconota", "Guerlinguetus aestuans", "Myoprocta acouchy",
                 "Dasyprocta fuliginosa", "Myoprocta pratti", "Dasyprocta leporina")
  
mydata4 <- mydata4 %>%
  filter(binomial %in% nao_sociais)
mydata4

mydata4 <- mydata4 %>%
  left_join(distinct(mydata[, c("cnuc", "nome_UC")]), by="cnuc") %>%
  relocate(nome.UC, .before = binomial)
mydata4

write.csv(mydata4, here("data", "ICMBio_Monitora_dados_selecionados.csv"), row.names = FALSE)

#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------

# Função lpi_icmbio original:

lpi_icmbio <- function(x,y,z) { # x = dados, y = UC, z = Classe 
  library(rlpi) # carrega pacote rlpi
  source(here("bin", "CalcLPI.R"))
  source(here("bin", "create_infile.R"))
  source(here("bin", "ProcessFile.R"))
  source(here("bin", "debug_print.R"))
  source(here("bin", "calculate_index.R"))
  source(here("bin", "bootstrap_lpi.R"))
  source(here("bin", "plot_lpi.R"))
  
  mydata <- x
  
  if(missing(y)) {    mydata <- mydata }
  else {  mydata <- subset(x, nome.UC == y) } # seleciona UC
  
  for.effort <- mydata # sera usado mais a frente
  
  if(missing(z)) {    mydata <- mydata }
  else {  mydata <- subset(mydata, Classe == z) } # seleciona taxon
  
  # colocar dados no formato exigido pela função LPIMain
  # passo 1, criar objeto a ser preenchido
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata$Binomial))))
  colnames(mydata2) <- c("ID", "Binomial", sort(unique(seq(min(mydata$Ano), max(mydata$Ano))))) # cria automaticamente os nomes de colunas de anos
  mydata2$ID <- c(1:nrow(mydata2))
  mydata2$Binomial <- sort(unique(mydata$Binomial))
  vetor_ano <- seq(min(mydata$Ano), max(mydata$Ano))
  
  # passo 2, preencher objeto criado acima
  for(i in 1:nrow(mydata2))
    for(j in 1:length(vetor_ano)){
      a <- subset(mydata, Binomial == mydata2[i,2])
      b <- subset(a, Ano == vetor_ano[j]) # extrai o ano automaticamente
      cduc <- unique(a$CDUC)
      c <- subset(for.effort, CDUC %in% cduc) # effort must be calculated separately per UC
      
      if ( nrow(subset(c, Ano == vetor_ano[j])) <= 0)  { mydata2[i,j+2] <- NA } else {
        if ( nrow(b) == 0)  { mydata2[i,j+2] <- 0 } else {
          mydata2[i,j+2] <- nrow(b)/(sum(subset(c, Ano == vetor_ano[j])$esforço, na.rm=TRUE)/10000)}
      }}
  
  
  #mydata2[, colSums(is.na(mydata2)) != nrow(mydata2)] # remover coluna de NAs se houver
  require(stringr) # necessario para passo abaixo
  colnames(mydata2)[3:ncol(mydata2)] <- str_c( "X", colnames(mydata2)[3:ncol(mydata2)]) # adiciona um "X" no nome das colunas de anos porque o create_infile exige
  assign("mydata", mydata, envir=globalenv())
  assign("mydata2", mydata2, envir=globalenv())
  
  # criar vetor indice com TRUE para todas as linhas, pois todas as espécies serão incluídas
  index_vector = rep(TRUE, nrow(mydata2))
  
  # criar infile
  mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = colnames(mydata2)[3], end_col_name = tail(colnames(mydata2), n=1), CUT_OFF_YEAR = vetor_ano[1])
  #mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = "X2014", end_col_name = "X2016", CUT_OFF_YEAR = 2014)
  #mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = "2014", end_col_name = "2017", CUT_OFF_YEAR = "2014")
  
  
  # Calcular LPI com 100 bootstraps
  source(here("bin", "LPIMain.R")) # adicionado por não ter pacote rlpi
  mydata_lpi <- LPIMain(mydata_infile, REF_YEAR = vetor_ano[1], PLOT_MAX = tail(vetor_ano, n=1), BOOT_STRAP_SIZE = 100, VERBOSE=FALSE)
  
  # Remover NAs (anos seguidos sem dados)
  mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
  
  # carregar função ggplot_lpi_modif
  source(here("bin", "ggplot_lpi_modif.R"))
  
  # Gerar gráfico mais bonito usando função ggplot_lpi_modif
  ggplot_lpi_modif(mydata_lpi, col="darkcyan")
  
  # salvar mydata_lpi na area de trabalho
  assign("mydata_lpi", mydata_lpi, .GlobalEnv)
  
} # Fim da função


#----------------------------------------
#----------------------------------------
#----------------------------------------
#----------------------------------------


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