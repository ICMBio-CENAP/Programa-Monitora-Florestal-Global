# Funcoes para analise de tendencias populacionais

encounter.rate <- function(mydata, taxon) {
  mydata[,taxon] <- factor(mydata[,taxon])
  # passo 1, criar objeto mydata2
  mydata2 <- data.frame(matrix(ncol = (2+length(seq(min(mydata$Ano), max(mydata$Ano)))), nrow = length(unique(mydata[,taxon]))) )
  colnames(mydata2) <- c("ID", "taxon", sort(unique(seq(min(mydata$Ano), max(mydata$Ano))))) # cria automaticamente os nomes de colunas de anos
  mydata2$ID <- c(1:nrow(mydata2))
  mydata2$taxon <- sort(unique(mydata[,taxon]))
  vetor.Ano <- seq(min(mydata$Ano), max(mydata$Ano))
  # passo 2, preencher objeto mydata2
  for(i in 1:nrow(mydata2))
    for(j in 1:length(vetor.Ano)){
      a <- subset(mydata, mydata[,taxon] == mydata2[i,2])
      b <- subset(a, Ano == vetor.Ano[j]) # extrai o ano automaticamente
      cduc <- unique(a$CDUC)
      c <- subset(mydata, CDUC %in% cduc) # effort must be calculated separately per UC
      
      if ( nrow(subset(c, Ano == vetor.Ano[j])) <= 0)  { mydata2[i,j+2] <- NA } else {
        if ( nrow(b) == 0)  { mydata2[i,j+2] <- 0 
        }
        else {
          mydata2[i,j+2] <- round(nrow(b)/(sum(subset(c, Ano == vetor.Ano[j])$esforço, na.rm=TRUE)/10000), 3)
        }
      }}
  #sitename <- deparse(substitute(mydata))
  #assign(paste("encounter_rate", sitename, sep="_"), mydata2, .GlobalEnv)
  assign("encounter_rate", mydata2, .GlobalEnv)
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
