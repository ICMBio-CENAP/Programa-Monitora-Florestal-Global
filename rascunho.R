  taxas_anuaisLonger <- taxas_anuais[,-1]
  taxas_anuais <-
    taxas_anuaisLonger %>%
    pivot_longer(!populacao, names_to = "ano", values_to = "index")
  
  taxas_anuais$cnuc <- substring(taxas_anuais$populacao,nchar(taxas_anuais$populacao)-2,nchar(taxas_anuais$populacao))
  taxas_anuais <- merge(taxas_anuais,ucs, by =  "cnuc")
  taxas_anuais$especie <- substr(taxas_anuais$populacao,1,nchar(taxas_anuais$populacao)-4)
  taxas_anuais$especie <- gsub("_"," ", taxas_anuais$especie)
  taxas_anuais$ano <- gsub("X","", taxas_anuais$ano)
  #taxas_anuaisTeste <- taxas_anuais[,2:8]
  library(here)
    library(dbplyr)
  library(tidyverse)
  
  mydata <- readRDS(here("03_dadosDeSaida/dados", "dadosICMBio_2014a2019.rds"))
  #Contagem de UCs
  ucs <- mydata %>%
    group_by(cnuc,nome_UC) %>%
    count()
  
  taxas_anuais <- read.csv(here("03_dadosDeSaida/dados", "taxas_anuais.csv"))
  taxas_anuais <- merge(taxas_anuais,ucs, by =  "cnuc")  
  
  mydata <- readRDS("../../03_dadosDeSaida/dados/dadosICMBio_2014a2019.rds")
  ucs <- mydata %>%
    group_by(cnuc,nome_UC) %>%
    count()
  taxas_anuais <- merge(taxas_anuais,ucs, by =  "cnuc")
  source("02_script/funcoes/LPIMain.R")
  source("02_script/funcoes/create_infile.R")
  source("02_script/funcoes/CalcLPI.R")
  source("02_script/funcoes/ProcessFile.R")
  source("02_script/funcoes/debug_print.R")
  source("02_script/funcoes/calculate_index.R")
  source("02_script/funcoes/bootstrap_lpi.R")
  source("02_script/funcoes/plot_lpi.R")
  source("02_script/funcoes/ggplot_lpi_modif.R")
  
  getwd()
  taxas_anuais <- read.csv("03_dadosDeSaida/dados/taxas_anuais.csv")
  
  mydata_lpi <- taxas_anuais %>%
    mutate(ID = as.numeric(1:nrow(taxas_anuais))) %>%
    rename(Binomial = populacao) %>%
    select(cnuc, nome_UC, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
  
  lpi <- data.frame(uc = character(),  ano = integer(), LPI_final = double(),CI_low = double(), CI_high = double())    
  
#Fazer o cálculo por UC  
    for(UC in unique(taxas_anuais$cnuc)){
      #cria um objeto temporário para receber o cálculo de LPI, extrai os anos a partir do nome das linhas
      lpiTemp <- tibble::rownames_to_column(mydata_lpi) %>%
        #renomea a coluna para ano
        rename(ano = rowname)
      #cria a coluna de UC e adiciona o nome armazenado
      lpiTemp$uc <- temp
      #adiciona os dados no objeto
      lpi <- rbind(lpitemp,lpi)


  
  lpi <- data.frame(ano = integer(), LPI_final = double(),CI_low = double(), CI_high = double())    
  
  tibble::rownames_to_column(mydata_lpi)
    

  ggplot(lpi, aes(x=factor(ano), y= LPI_final, group=1), height = 15)+
    geom_point()+
    geom_line()+
    geom_ribbon(aes(ymin = CI_low, ymax=CI_high), alpha=0.5, fill = "blue" ,color = "dark grey")+
    
    geom_smooth(method = "loess", se = FALSE,
                data =lpi,
                aes(x= factor(ano), y = LPI_final, y=Temp), 
                alpha=1,  size=0.5, color = "#00798c")
    
    theme(strip.text.y = element_text(size=7, angle = 360, face="bold"))+
    facet_grid(vars(uc))+#uc ~ .)+
    #facet_wrap(facets = vars(uc))+
  #  options(facet_size_manual = list(width = c(1,5), height = 5))+
    xlab("UC")+
    
    ylab("LPI")+

        theme(legend.position="top")+
    
    theme(axis.text.x = element_text(
      angle=45,
      vjust=1,
      hjust = 1,
      colour="black",
      size=rel(1))
    )
  
  
  theme(strip.text.x = element_text(size=8, angle=75),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))
  
  
  }
  

geom_smooth(method="loess", se=T) 


#PEGUE TAXAS ANUAIS PARA FAZER UMA TABELA
taxas_anuaisTable <- taxas_anuais
#RETIRE CNUC
#taxas_anuaisLonger <- taxas_anuaisTable[,-1]
#TRANSPÕE
#taxas_anuaisTable <-
#   taxas_anuaisLonger %>%
#  tidyr::pivot_longer(
#   !populacao,
#  names_to = "ano",
# values_to = "index"
#)
#PEGUE O CÓDIGO DA UC
#taxas_anuaisTable$cnuc <- substring(
# taxas_anuaisTable$populacao,
#nchar(taxas_anuaisTable$populacao)-2,
#nchar(taxas_anuaisTable$populacao)
#)
#RETIRE _
#taxas_anuaisTable$cnuc <- as.numeric(gsub("_","", taxas_anuaisTable$cnuc))
#PEGUE O NOME DAS UCs
taxas_anuaisTable <- merge(taxas_anuaisTable,ucs, by =  "cnuc")
#PEGUE O NOME DAS ESPÉCIES
taxas_anuaisTable$especie <- substr(
  taxas_anuaisTable$populacao,
  1,
  nchar(taxas_anuaisTable$populacao)-4
)
#RETIRE _
taxas_anuaisTable$especie <- gsub("_"," ", taxas_anuaisTable$especie)
#RETIRE X DE ANO
#taxas_anuaisTable$ano <- as.factor(gsub("X","", taxas_anuaisTable$ano))

taxas_anuaisTable <- taxas_anuaisTable %>%
  dplyr::select(
    nome_UC,
    ano,
    especie,
    index
  )#%>%
dplyr::rename(
  "Unidade de Conservação" = nome_UC,
  Ano = ano,
  "Espécie" = especie
)

