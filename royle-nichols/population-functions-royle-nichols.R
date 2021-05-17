# Funcoes para analise de tendencias populacionais

encounter.rate <- function(mydata, taxon) {
  mydata[,taxon] <- factor(mydata[,taxon])
  # passo 1, criar objeto mydata2
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$ano), max(mydata$ano)))), nrow = length(unique(mydata[,taxon]))) )
  colnames(mydata2) <- c("ID", "taxon", sort(unique(seq(min(mydata$ano), max(mydata$ano))))) # cria automaticamente os nomes de colunas de anos
  mydata2$ID <- c(1:nrow(mydata2))
  mydata2$taxon <- sort(unique(mydata[,taxon]))
  vetor.ano <- seq(min(mydata$ano), max(mydata$ano))
  # passo 2, preencher objeto mydata2
  for(i in 1:nrow(mydata2))
    for(j in 1:length(vetor.ano)){
      a <- subset(mydata, mydata[,taxon] == mydata2[i,2])
      b <- subset(a, ano == vetor.ano[j]) # extrai o ano automaticamente
      cduc <- unique(a$CDUC)
      c <- subset(mydata, CDUC %in% cduc) # effort must be calculated separately per UC
      
      if ( nrow(subset(c, ano == vetor.ano[j])) <= 0)  { mydata2[i,j+2] <- NA } else {
        if ( nrow(b) == 0)  { mydata2[i,j+2] <- 0 
        }
        else {
          mydata2[i,j+2] <- round(nrow(b)/(sum(subset(c, ano == vetor.ano[j])$esforço, na.rm=TRUE)/10000), 3)
        }
      }}
  #sitename <- deparse(substitute(mydata))
  #assign(paste("encounter_rate", sitename, sep="_"), mydata2, .GlobalEnv)
  assign("encounter_rate", mydata2, .GlobalEnv)
}


make.binomial.table <- function(x) {
  binomial.table <- tibble( station.year=unique(x$trilha_ano),
                            station= substr(station.year, start = 1, stop = 4),
                            year=as.numeric(substr(station.year, start = 6, stop = nchar(station.year))), 
                            trials=as.numeric(NA), 
                            success=as.numeric(0) )
  for(i in 1:nrow(binomial.table)) {
    df1 <- x %>% filter(trilha_ano == binomial.table$station.year[i])
    binomial.table$trials[i] <- length(unique(df1$data.da.amostragem))
    for(j in 1:length(unique(df1$data.da.amostragem))) {
      if("Dasyprocta" %in% df1$Genero[j]) {
        binomial.table$success[i] <- binomial.table$success[i]+1
      }
    }
  }
  return(binomial.table)
  assign("binomial.table", binomial.table, .GlobalEnv)
  #print(binomial.table, n=Inf)
  #y <- binomial.table %>% pivot_wider(names_from = year, values_from = success)
}


prepare.data.RNmodel <- function(x, taxon) {
  y <- tibble(station=unique(x$estacao.amostral))
  vec <- unique(x$ano)
  df <- bind_rows(setNames(rep(0, length(vec)), vec))
  y <- bind_cols(y, df)
  y
  for(i in 1:nrow(y)) {
    for(j in 2:ncol(y)) {
      df1 <- x %>% filter(estacao.amostral == y[[i,1]])
      df2 <- df1 %>% filter(ano == colnames(y)[j])
      days <- unique(df2$data.da.amostragem)
      for(k in 1:length(unique(days))) {
        df3 <- df2 %>% filter(data.da.amostragem == days[k])
        if(taxon %in% df3$Genero) { # ATENCAO: por enquanto escolha de taxon manual
          y[[i,j]] <- y[[i,j]]+1
        }
      }
    }}
  #return(y) # checar
  assign("y", y, .GlobalEnv)
}




state.space.model <- function(y, n.years) {
  
  # Specify model in BUGS language
  sink(here("experimental", "ssm.jags"))
  cat("
model {
# Priors and constraints
logN.est[1] ~ dnorm(0, 0.01)       # Prior for initial population size
mean.r ~ dnorm(1, 0.001)             # Prior for mean growth rate
sigma.proc ~ dunif(0, 1)             # Prior for sd of state process
sigma2.proc <- pow(sigma.proc, 2)
tau.proc <- pow(sigma.proc, -2)
sigma.obs ~ dunif(0, 1)              # Prior for sd of observation process
sigma2.obs <- pow(sigma.obs, 2)
tau.obs <- pow(sigma.obs, -2)

# Likelihood
# State process
for (t in 1:(T-1)){
   r[t] ~ dnorm(mean.r, tau.proc)
   logN.est[t+1] <- logN.est[t] + r[t]
   }
# Observation process
for (t in 1:T) {
   y[t] ~ dnorm(logN.est[t], tau.obs)
   }

# Population sizes on real scale
for (t in 1:T) {
   N.est[t] <- exp(logN.est[t])/100
   }
}
",fill = TRUE)
  sink()
  
  # definir numero de anos
  n.years <- length(3:ncol(encounter_rate))
  
  # Bundle data
  jags.data <- list(y = log(y*100), T = n.years)
  
  # Initial values
  inits <- function(){list(sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), 
                           sigma.obs = runif(1, 0, 1),
                           LogN.est = c(rnorm(1, -0.5, 0.1), rep(NA, (n.years-1))))} 
  #LogN.est = c(rnorm(n.years, 5, 0.1)) )} 
  
  # Parameters monitored
  parameters <- c("r", "mean.r", "sigma2.obs", "sigma2.proc", "N.est")
  
  # MCMC settings
  ni <- 25000
  nt <- 3
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R (BRT <1 min)
  ssm <- jags(jags.data, inits, parameters, here("experimental", "ssm.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # ccheck results
  print(ssm, digits = 2)
  
  # Probability of N(2019) < N(2014)
  mean(ssm$BUGSoutput$sims.list$N.est[,6] < ssm$BUGSoutput$mean$N.est[1])
  
  assign("meanR", ssm$BUGSoutput$sims.list$mean.r, .GlobalEnv) 
  assign("ssm", ssm, .GlobalEnv)

}

# Draw figure
pop.trends <- function() { 
  fitted <- lower <- upper <- numeric()
  year <- 2014:2019
  n.years <- length(3:ncol(encounter_rate))
  
  for (i in 1:n.years){
    fitted[i] <- mean(ssm$BUGSoutput$sims.list$N.est[,i])
    lower[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.025)
    upper[i] <- quantile(ssm$BUGSoutput$sims.list$N.est[,i], 0.975)
  }
  m1 <- min(c(fitted, y, lower), na.rm = TRUE)
  m2 <- max(c(fitted, y, upper), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  #plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  #plot(0, 0, ylim = c(m1-0.5, m2+0.5), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  plot(0, 0, ylim = c(0, m2+(mean(fitted)*0.5)), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:n.years, labels = year)
  #polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  #points(y, type = "l", col = "black", lwd = 1, lty = 2)
  #points(fitted, type = "l", col = "blue", lwd = 2)
  points(x = (1:n.years), y = fitted, type = "b", pch = 16, cex = 1.5, lty = 1)
  segments((1:n.years), lower, 1:(n.years), upper, cex=0.5)
}


# Draw figure
pop.trends.for.rmd <- function(x) { 
  fitted <- lower <- upper <- numeric()
  year <- 2014:2019
  n.years <- length(3:ncol(encounter_rate))
  
  for (i in 1:n.years){
    fitted[i] <- mean(x[,i])
    lower[i] <- quantile(x[,i], 0.025)
    upper[i] <- quantile(x[,i], 0.975)
  }
  m1 <- min(c(fitted, y, lower), na.rm = TRUE)
  m2 <- max(c(fitted, y, upper), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  #plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  #plot(0, 0, ylim = c(m1-0.5, m2+0.5), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  plot(0, 0, ylim = c(0, m2+(mean(fitted)*0.5)), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:n.years, labels = year)
  #polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  #points(y, type = "l", col = "black", lwd = 1, lty = 2)
  #points(fitted, type = "l", col = "blue", lwd = 2)
  points(x = (1:n.years), y = fitted, type = "b", pch = 16, cex = 1.5, lty = 1)
  segments((1:n.years), lower, 1:(n.years), upper, cex=0.5)
}


#----
# Funcao para criar historico de deteccao

detection.history <- function(data){ # occ_length is length of sampling occasions
  require(dplyr)
  require(lubridate)
  
  #results object
  res<-list()
  
  #get the dimensions of the matrix
  
  #list of sampling units
  taxa <- sort(unique(data$Genero))
  station <- as.factor(unique(data$estacao.amostral))
  station <- sort(station)
  rows <- length(station)
  cols <- length(unique(data$data.da.amostragem))
  
  #sampling period
  date.header <- sort(unique(data$data.da.amostragem))
  mat<-matrix(NA,rows,cols,dimnames=list(station,as.character(date.header)))
  
  mat.template<-mat
  
  dates <- data$data.da.amostragem[indx]
  sites <- data$estacao.amostral[indx]
  dates.sites<-data.frame(dates, sites)
  
  #outline the days in which which site was sampled
  df1 <- distinct(data, estacao.amostral, data.da.amostragem)
  for(i in 1:length(station)) {
    for(j in 1:ncol(mat.template)){
    indx <- which(df1$estacao.amostral==station[i])
    dates <- df1$data.da.amostragem[indx]
    station <- df1$estacao.amostral[indx]
    dates.station<-data.frame(dates, station)
    dates.station<-unique(dates.station)
    col<-which(date.header==dates.station[j,1])
    mat.template[i,col]<-0
      }
    }
    
  mat <- mat.template
  
  #insert detections for each taxa i
  for(i in 1:length(taxa)){
    indx <- which(data$Genero==taxa[i])
    #dates and sites when/where the individuals were recorded
    dates <- data$data.da.amostragem[indx]
    sites <- data$estacao.amostral[indx]
    dates.sites<-data.frame(dates, sites)
    #unique combination of dates and sites 
    dates.sites<-unique(dates.sites)
    #fill in the matrix
    for(j in 1:length(dates.sites[,1])){
      col<-which(date.header==dates.sites[j,1])
      row<-which(trap==dates.sites[j,2])
      mat[row,col]<-1
    }
    mat.nas<-is.na(mat)
    sum.nas<-apply(mat.nas,2,sum)
    indx.nas<-which(sum.nas==rows)
    if(length(indx.nas)>0){
      mat<-mat[,-indx.nas]
    }
    # reduce the size of the matrix
    #mat <- f.shrink.matrix(mat, occ_length, data)
    res<-c(res,list(mat))
    #return the matrix to its original form
    mat<-mat.template
  }
  
  names(res) <- taxa
  #res<-lapply(res,f.dum)
  res
  
}


# Draw figure
pop.trends <- function(x,y) { 
  fitted <- lower <- upper <- numeric()
  year <- as.numeric(gsub("X", "", colnames(y)))
  n.years <- ncol(y)
  
  for (i in 1:n.years){
    fitted[i] <- x$BUGSoutput$mean$mean.abundance[i]
    lower[i] <- quantile(x$BUGSoutput$sims.list$mean.abundance[,i], 0.025)
    upper[i] <- quantile(x$BUGSoutput$sims.list$mean.abundance[,i], 0.975)
  }
  m1 <- min(c(fitted, lower), na.rm = TRUE)
  m2 <- max(c(fitted, upper), na.rm = TRUE)
  par(mar = c(4.5, 4, 1, 1))
  #plot(0, 0, ylim = c(m1, m2), xlim = c(1, n.years), ylab = "Taxa de encontro", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  #plot(0, 0, ylim = c(m1-0.5, m2+0.5), xlim = c(1, n.years), ylab = "Taxa de encontro (Ind/10km)", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  plot(0, 0, ylim = c(0, m2+(mean(fitted)*0.5)), xlim = c(1, n.years), ylab = "Abundância média", xlab = "Ano", col = "black", type = "l", lwd = 2, axes = FALSE, frame = FALSE)
  axis(2, las = 1)
  axis(1, at = 1:n.years, labels = year)
  #polygon(x = c(1:n.years, n.years:1), y = c(lower, upper[n.years:1]), col = "gray90", border = "gray90")
  #points(y, type = "l", col = "black", lwd = 1, lty = 2)
  #points(fitted, type = "l", col = "blue", lwd = 2)
  points(x = (1:n.years), y = fitted, type = "b", pch = 16, cex = 1.5, lty = 1)
  segments((1:n.years), lower, 1:(n.years), upper, cex=0.5)
}

