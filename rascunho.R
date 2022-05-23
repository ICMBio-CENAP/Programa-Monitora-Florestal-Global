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
  
  
  library(ggplot2)
  
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

    mydata_lpi <- taxas_anuais[taxas_anuais$cnuc==UC,] %>%
      mutate(ID = 1:nrow(taxas_anuais[taxas_anuais$cnuc==UC,]))%>%
      rename(Binomial = populacao) %>%
      select(cnuc,nome_UC, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
    #armazena o nome da uc
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
     ## png(paste('"../../03_dadosDeSaida/plots',temp,'.png', sep = ""))#, width = 10, height = 7, units = 'cm', res = 72)
      ##g = ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = temp)
      
      #g = ggplot(invTemp) + geom_histogram(aes(dap2012), binwidth = 10) + 
       # xlab('Centro de classe (cm)') + ylab('# árvores') + xlim(0, 200) +
        #theme_light() 
      ##print(g)
      ##dev.off()
    }
    #mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
    #ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = UC)

  
  
  
  
  
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
  
  