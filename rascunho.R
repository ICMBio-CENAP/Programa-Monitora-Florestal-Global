

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

for(UC in unique(taxas_anuais$nome_UC)){
  for(SP in unique(taxas_anuais$especie)){
  dev.new()
print(
ggplot(taxas_anuais[taxas_anuais$nome_UC==UC,], aes(x=ano, y=indeX))+
  geom_bar(
    position = position_dodge2(preserve ="single"),
    stat="identity",
    width = 0.5,
    size =0.3,
    fill = "blue"
  )+
  
  #theme_bw()+
  
  #facet_wrap(facets = vars(especie))+
  
  xlab("CASTANHAL")+
  
  ylab("MÉDIA DE SEMENTES")+
  
  theme(legend.position="top")+
  
  theme(axis.text.x = element_text(
    angle=45,
    vjust=1,
    hjust = 1,
    colour="black",
    size=rel(1))
  )
  }
  
  ls("package:tidyverse")
  search()

  
  
  for(UC in unique(taxas_anuais$cnuc)){
    mydata_lpi <- taxas_anuais[taxas_anuais$cnuc==UC,] %>%
      mutate(ID = 1:nrow(taxas_anuais[taxas_anuais$cnuc==UC,]))%>%
      rename(Binomial = populacao) %>%
      select(cnuc, ID, Binomial, X2014, X2015, X2016, X2017, X2018, X2019)
    #print(mydata_lpi)
    
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
      BOOT_STRAP_SIZE = 100) #, VERBOSE=FALSE)
    
    for(pop in 1:nrow(mydata_lpi)){
      mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
      png(paste('"../../03_dadosDeSaida/plots',pop,'.png', sep = ""), width = 10, height = 7, units = 'cm', res = 72)
      g = ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = pop)
      
      #g = ggplot(invTemp) + geom_histogram(aes(dap2012), binwidth = 10) + 
       # xlab('Centro de classe (cm)') + ylab('# árvores') + xlim(0, 200) +
        #theme_light() 
      print(g)
      dev.off()
    }
    #mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
    #ggplot_lpi_modif(mydata_lpi, col="darkcyan", title = UC)
  }
  