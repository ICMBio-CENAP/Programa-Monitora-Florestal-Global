# Função lpi_icmbio:
lpi_icmbio <- function(x,y,z) # x = dados, y = UC, z = Classe
{
  #library(rlpi) # carrega pacote rlpi
  source("/media/elildojr/Dados/r/lpi/rlpi-master/R/CalcLPI.R") # carrega função direto da pasta
  source("create_infile.R") # idem
  source("ProcessFile.R")
  source("debug_print.R")
  source("calculate_index.R")
  source("bootstrap_lpi.R")
  source("plot_lpi.R")

  for.effort <- x # sera usado mais a frente
  mydata <- x
  
  if(missing(y)) {    mydata <- mydata }
  else {  mydata <- subset(x, Local...Nome.da.Unidade.de.Conservação == y) } # seleciona UC
  
  if(missing(z)) {    mydata <- mydata }
  else {  mydata <- subset(mydata, Classe == z) } # seleciona taxon

  # colocar dados no formato exigido pela função LPIMain
  # passo 1, criar objeto a ser preenchido
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata$Binomial))))
  colnames(mydata2) <- c("ID", "Binomial", sort(unique(seq(min(mydata$Ano), max(mydata$Ano))))) # cria automaticamente os nomes de colunas de anos
  mydata2$ID <- c(1:nrow(mydata2))
  mydata2$Binomial <- sort(unique(mydata$Binomial))
  vetor.Ano <- seq(min(mydata$Ano), max(mydata$Ano))

  # passo 2, preencher objeto criado acima
  for(i in 1:nrow(mydata2))
    for(j in 1:length(vetor.Ano)){
      a <- subset(mydata, Binomial == mydata2[i,2])
      b <- subset(a, Ano == vetor.Ano[j]) # extrai o ano automaticamente
        
            # esforço anual por UC 
            c <- filter(for.effort, CDUC %in% b$CDUC)
            annual.effort <- rep(NA, length(vetor.Ano))
                for(k in 1:length(annual.effort)){ annual.effort[k] <- sum(subset(c, Ano == vetor.Ano[k])$esforço, na.rm=TRUE)/10000 }
            
      if ( annual.effort[j] <= 0) { mydata2[i,j+2] <- NA }
      if ( nrow(b) == 0) { mydata2[i,j+2] <- 0 }
      else { mydata2[i,j+2] <- nrow(b)/annual.effort[j]}
    }
  
  #mydata2[, colSums(is.na(mydata2)) != nrow(mydata2)] # remover coluna de NAs se houver
  require(stringr) # necessario para passo abaixo
  colnames(mydata2)[3:ncol(mydata2)] <- str_c( "X", colnames(mydata2)[3:ncol(mydata2)]) # adiciona um "X" no nome das colunas de anos porque o create_infile exige
  assign("mydata", mydata, envir=globalenv())
  assign("mydata2", mydata2, envir=globalenv())
  
  # criar vetor indice com TRUE para todas as linhas, pois todas as espécies serão incluídas
  index_vector = rep(TRUE, nrow(mydata2))
  
  # criar infile
  mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = colnames(mydata2)[3], end_col_name = tail(colnames(mydata2), n=1), CUT_OFF_YEAR = vetor.Ano[1])
  #mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = "X2014", end_col_name = "X2016", CUT_OFF_YEAR = 2014)
  #mydata_infile <- create_infile(mydata2, index_vector=index_vector, name="mydata_data", start_col_name = "2014", end_col_name = "2017", CUT_OFF_YEAR = "2014")
  
  
  # Calcular LPI com 100 bootstraps
  source("LPIMain.R") # adicionado por não ter pacote rlpi
  mydata_lpi <- LPIMain(mydata_infile, REF_YEAR = vetor.Ano[1], PLOT_MAX = tail(vetor.Ano, n=1), BOOT_STRAP_SIZE = 100, VERBOSE=FALSE)
  
  # Remover NAs (anos seguidos sem dados)
  mydata_lpi <- mydata_lpi[complete.cases(mydata_lpi), ]
  
  # carregar função ggplot_lpi_modif
  source("ggplot_lpi_modif.R")
  
  # Gerar gráfico mais bonito usando função ggplot_lpi_modif
  ggplot_lpi_modif(mydata_lpi, col="cornflowerblue")
  
} # Fim da função