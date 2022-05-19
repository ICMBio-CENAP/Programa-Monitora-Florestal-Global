
taxas_anuais1 <- pivot_longer(taxas_anuais,
                              cols = starts_with("X"),
                              names_to = "populacao") 
taxas_anuaisLonger <- taxas_anuais[, -c(1,2)]
taxas_anuais <-
taxas_anuaisLonger %>%
    pivot_longer(!populacao, names_to = "ano", values_to = "index")

taxas_anuais$cnuc <- substring(taxas_anuais$populacao,nchar(taxas_anuais$populacao)-2,nchar(taxas_anuais$populacao))
taxas_anuais <- merge(taxas_anuais,ucs, by =  "cnuc")
taxas_anuais$especie <- substr(taxas_anuais$populacao,1,nchar(taxas_anuais$populacao)-4)
taxas_anuais$especie <- gsub("_"," ", taxas_anuais$especie)
taxas_anuais$ano <- gsub("X","", taxas_anuais$ano)
#taxas_anuaisTeste <- taxas_anuais[,2:8]

install.packages("ggsubplot")
library(ggsubplot)
ggplot(taxas_anuais)+
  geom_subpl


  mydata_lpiUC <- taxas_anuais[taxas_anuais$cnuc==118,] %>%
    mutate(ID = 1:nrow(taxas_anuais[taxas_anuais$cnuc==118,]))%>%
    rename(Binomial = populacao) %>%
    select(cnuc, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
  
  mydata_lpi <- mydata_lpiUC[mydata_lpiUC$Binomial=="Ateles_chamek_118",]
  #print(mydata_lpi)
  temp <- mydata_lpi$Binomial
  index_vector = rep(TRUE, nrow(mydata_lpi))
  year_vector <- 2014:2018
  
  
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
  
  source(here("02_script/funcoes", "LPIMain.R"))
  source(here("02_script/funcoes", "create_infile.R"))
  source(here("02_script/funcoes", "CalcLPI.R"))
  source(here("02_script/funcoes", "ProcessFile.R"))
  source(here("02_script/funcoes", "debug_print.R"))
  source(here("02_script/funcoes", "calculate_index.R"))
  source(here("02_script/funcoes", "bootstrap_lpi.R"))
  source(here("02_script/funcoes", "plot_lpi.R"))
  source(here("02_script/funcoes", "ggplot_lpi_modif.R"))
  
  
  library(ggplot2)


  
  lpi <- data.frame(ano = integer(), LPI_final = double(),CI_low = double(), CI_high = double())    
  
  tibble::rownames_to_column(mydata_lpi)
    
    for(UC in unique(taxas_anuais$cnuc)){
  
lpiTemp <- rbind(lpi,mydata_lpi)

      
    mydata_lpi <- taxas_anuais[taxas_anuais$cnuc==UC,] %>%
      mutate(ID = 1:nrow(taxas_anuais[taxas_anuais$cnuc==UC,]))%>%
      rename(Binomial = populacao) %>%
      select(cnuc,nome_UC, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
    #print(mydata_lpi)
    #temp <- mydata_lpi$Binomial
    temp <- unique(mydata_lpi$nome_UC)
    index_vector = rep(TRUE, nrow(mydata_lpi))
    year_vector <- 2014:2018
    
    mydata_infile <- create_infile(mydata_lpi, index_vector=index_vector, 
                                   #name=here("lpi", "mydata_data"),
                                   name = "Infile",
                                   start_col_name = colnames(mydata_lpi)[4], end_col_name = tail(colnames(mydata_lpi), n=1),
                                   CUT_OFF_YEAR = year_vector[1])
    
    mydata_lpi <- LPIMain(#infile = here("lpi", "mydata_data_infile.txt"),
      infile = "Infile_infile.txt",
      #basedir = here("lpi"),
      REF_YEAR = year_vector[1],
    #  PLOT_MAX = tail(year_vector, n=1),
      BOOT_STRAP_SIZE = 1000) #, VERBOSE=FALSE)
    
    
      #mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
      png(paste('"../../03_dadosDeSaida/plots',temp,'.png', sep = ""))#, width = 10, height = 7, units = 'cm', res = 72)
      g = ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = temp)
      
      #g = ggplot(invTemp) + geom_histogram(aes(dap2012), binwidth = 10) + 
       # xlab('Centro de classe (cm)') + ylab('# árvores') + xlim(0, 200) +
        #theme_light() 
      print(g)
      dev.off()
    }
    #mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
    #ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = UC)
  }
  