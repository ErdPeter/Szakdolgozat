library(survival)
library(gaston)
library(dplyr)
library(Rlab)
# library(SPACox)
# library(gwasurvivr)

setwd("c:/Tanulas/Szakdolgozat/Ero2")

### Lambda szimulálása
set.seed(1234)

N = 10000 # mintaszám
gamma <- 0 # genotipic effect
lambdak <- runif(1000, 0, 20)  # event rate

ER <- data.frame(lambda = lambdak,
                 ER = 0)

for(l in 1:length(lambdak)){
  
  lambda <- lambdak[l]
  res <- data.frame(C = rweibull(n = N, shape = 1, scale = 0.15),
                    Tx = 0)
  U <- runif(N, 0, 1)
  X1 <- rnorm(N)
  X2 <- rbern(N, 0.5)
  
  for (i in 1:N){
    nu <- (0.5 * X1[i]) + (0.5 * X2[i])
    res[i, "Tx"] <- lambda * (sqrt(-log(U[i]) / exp(nu)))
  }
  
  res <- res %>% 
    mutate(event = ifelse(C <= Tx, 0, 1))
  ER$ER[l] <- sum(res$event)/N
}

ER <- ER[order(ER$lambda),]
plot(ER$ER, ER$lambda)

# 1000 minta, 10000 lambda
save(ER, file = paste0("lambda_becsles.RData"))
#load(file = paste0("lambda_becsles.RData"))

# 10000 minta, 1000 lambda
save(ER, file = paste0("lambda_becsles2.RData"))
#load(file = paste0("lambda_becsles2.RData"))

#becsles 2-ből
ER1 <- ER[ER$ER >= 0.005 & ER$ER <= 0.015, ]
ER5 <- ER[ER$ER >= 0.001 & ER$ER <= 0.003, ]
#becsles 1-ből
ER2 <- ER[ER$ER >= 0.08 & ER$ER <= 0.12, ]
ER3 <- ER[ER$ER >= 0.1 & ER$ER <= 0.3, ]
ER4 <- ER[ER$ER >= 0.3 & ER$ER <= 0.7, ]


with (ER1, plot(lambda ~ ER, pch = 20))
expreg <- lm(log(lambda) ~ ER, data = ER1)
with (ER1, plot(log(lambda) ~ ER, pch = 20))
lines(ER1$ER, predict(expreg), col = "red", lwd = 2)
new_data <- data.frame(ER = 0.01)
lambda1 <- predict(expreg, newdata = new_data)

with (ER2, plot(lambda ~ ER, pch = 20))
expreg <- lm(log(lambda) ~ ER, data = ER2)
with (ER2, plot(log(lambda) ~ ER, pch = 20))
lines(ER2$ER, predict(expreg), col = "red", lwd = 2)
new_data <- data.frame(ER = 0.1)
lambda2 <- predict(expreg, newdata = new_data)

with (ER3, plot(lambda ~ ER, pch = 20))
expreg <- lm(log(lambda) ~ ER, data = ER3)
with (ER3, plot(log(lambda) ~ ER, pch = 20))
lines(ER3$ER, predict(expreg), col = "red", lwd = 2)
new_data <- data.frame(ER = 0.2)
lambda3 <- predict(expreg, newdata = new_data)

with (ER4, plot(lambda ~ ER, pch = 20))
expreg <- lm(log(lambda) ~ ER, data = ER4)
with (ER4, plot(log(lambda) ~ ER, pch = 20))
lines(ER4$ER, predict(expreg), col = "red", lwd = 2)
new_data <- data.frame(ER = 0.5)
lambda4 <- predict(expreg, newdata = new_data)

with (ER5, plot(lambda ~ ER, pch = 20))
expreg <- lm(log(lambda) ~ ER, data = ER5)
with (ER5, plot(log(lambda) ~ ER, pch = 20))
lines(ER5$ER, predict(expreg), col = "red", lwd = 2)
new_data <- data.frame(ER = 0.002)
lambda5 <- predict(expreg, newdata = new_data)



lambdak <- data.frame(lambda = exp(c(lambda5, lambda1, lambda2, lambda3, lambda4)),
                      ER = c(0.002, 0.01, 0.1, 0.2, 0.5))

save(lambdak, file = paste0("lambdak.RData"))
#load(file = paste0("lambdak.RData"))


###############################################################################
###############################################################################
###############################################################################
### Szimulált adatok
set.seed(1234)

load(file = paste0("lambdak.RData"))

# MAF: 0.01, 0.05, 0.3
db_adatszett <- 50
MAF = 0.01 # MAF
N = 100000 # mintaszám
gammak_l <- list(seq(from = 0, to = 3, length.out = 11),
                 seq(from = 0, to = 1.5, length.out = 11),
                 seq(from = 0, to = 0.6, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.25, length.out = 11))# genotipic effect
geno_alt <- 20


for (l in 1:5){
  rm(list=setdiff(ls(), c("MAF", "N", "gammak_l", "geno_alt", "lambdak", "l",
                          "db_adatszett")))
  gc()
  
  lambda <- lambdak[l, 1] 
  data <- data.frame()
  gammak <- gammak_l[[l]]
  
  start <- Sys.time()
  
  # Genotípusok szimulálása
  geno <- matrix(rbinom(N * geno_alt * db_adatszett, size = 2, prob = MAF),
                 nrow = N * db_adatszett,
                 ncol = geno_alt)
  geno <- as.data.frame(geno)
  colnames(geno) <- paste0("G", 1:ncol(geno))
  
  for(g in 1:length(gammak)){
    gc()
    gamma <- gammak[g]
    
    # Idő szimulálása
    res <- data.frame(C = rweibull(n = N * db_adatszett, shape = 1, scale = 0.15),
                      X1 = rnorm(N * db_adatszett),
                      X2 = rbern(N * db_adatszett, 0.5),
                      gamma = gamma,
                      lambda = lambda,
                      U = runif(N * db_adatszett, 0, 1))
    res$SGX <- rowSums(geno)
    res$nu <- (0.5 * res$X1) + (0.5 * res$X2) + (res$gamma * res$SGX)
    res$Tx <- res$lambda * sqrt(-log(res$U) / exp(res$nu))
    res <- res %>% 
      mutate(time = ifelse(C <= Tx, C, Tx),
             event = ifelse(C <= Tx, 0, 1))
    res[, c("SGX", "nu", "U", "lambda", "C", "Tx")] <- NULL
    
    data <- rbind(data, res)
  }
  
  filename <- paste0("MAF", MAF, "_Lambda", lambdak[l, 2])
  filename <- sub("\\.", "_", filename)
  filename <- sub("\\.", "_", filename)
  saveRDS(data, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, ".rds"),
          compress = "xz")
  saveRDS(geno, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, "_geno.rds"),
          compress = "xz")
  
  stop <- Sys.time()
  print(stop - start)
  
}



###############################################################################
### Szimulált adatok
set.seed(1234)

load(file = paste0("lambdak.RData"))

# MAF: 0.01, 0.05, 0.3
db_adatszett <- 50
MAF = 0.05 # MAF
N = 100000 # mintaszám
gammak_l <- list(seq(from = 0, to = 2, length.out = 11),
                 seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.3, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11))# genotipic effect
geno_alt <- 20

for (l in 1:5){
  rm(list=setdiff(ls(), c("MAF", "N", "gammak_l", "geno_alt", "lambdak", "l",
                          "db_adatszett")))
  gc()
  
  lambda <- lambdak[l, 1]
  data <- data.frame()
  gammak <- gammak_l[[l]]
  
  start <- Sys.time()
  
  # Genotípusok szimulálása
  geno <- matrix(rbinom(N * geno_alt * db_adatszett, size = 2, prob = MAF),
                 nrow = N * db_adatszett,
                 ncol = geno_alt)
  geno <- as.data.frame(geno)
  colnames(geno) <- paste0("G", 1:ncol(geno))
  
  for(g in 1:length(gammak)){
    gc()
    gamma <- gammak[g]
    
    # Idő szimulálása
    res <- data.frame(C = rweibull(n = N * db_adatszett, shape = 1, scale = 0.15),
                      X1 = rnorm(N * db_adatszett),
                      X2 = rbern(N * db_adatszett, 0.5),
                      gamma = gamma,
                      lambda = lambda,
                      U = runif(N * db_adatszett, 0, 1))
    res$SGX <- rowSums(geno)
    res$nu <- (0.5 * res$X1) + (0.5 * res$X2) + (res$gamma * res$SGX)
    res$Tx <- res$lambda * sqrt(-log(res$U) / exp(res$nu))
    res <- res %>%
      mutate(time = ifelse(C <= Tx, C, Tx),
             event = ifelse(C <= Tx, 0, 1))
    res[, c("SGX", "nu", "U", "lambda", "C", "Tx")] <- NULL
    
    data <- rbind(data, res)
  }
  
  filename <- paste0("MAF", MAF, "_Lambda", lambdak[l, 2])
  filename <- sub("\\.", "_", filename)
  filename <- sub("\\.", "_", filename)
  saveRDS(data, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, ".rds"),
          compress = "xz")
  saveRDS(geno, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, "_geno.rds"),
          compress = "xz")
  
  stop <- Sys.time()
  print(stop - start)
  
}


###############################################################################
### Szimulált adatok
set.seed(1234)

load(file = paste0("lambdak.RData"))

# MAF: 0.01, 0.05, 0.3
db_adatszett <- 50
MAF = 0.3 # MAF
N = 100000 # mintaszám
gammak_l <- list(seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.5, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11),
                 seq(from = 0, to = 0.1, length.out = 11),
                 seq(from = 0, to = 0.06, length.out = 11))# genotipic effect
geno_alt <- 20


for (l in 1:5){
  rm(list=setdiff(ls(), c("MAF", "N", "gammak_l", "geno_alt", "lambdak", "l",
                          "db_adatszett")))
  gc()
  
  lambda <- lambdak[l, 1]
  data <- data.frame()
  gammak <- gammak_l[[l]]
  
  start <- Sys.time()
  
  # Genotípusok szimulálása
  geno <- matrix(rbinom(N * geno_alt * db_adatszett, size = 2, prob = MAF),
                 nrow = N * db_adatszett,
                 ncol = geno_alt)
  geno <- as.data.frame(geno)
  colnames(geno) <- paste0("G", 1:ncol(geno))
  
  for(g in 1:length(gammak)){
    gc()
    gamma <- gammak[g]
    
    # Idő szimulálása
    res <- data.frame(C = rweibull(n = N * db_adatszett, shape = 1, scale = 0.15),
                      X1 = rnorm(N * db_adatszett),
                      X2 = rbern(N * db_adatszett, 0.5),
                      gamma = gamma,
                      lambda = lambda,
                      U = runif(N * db_adatszett, 0, 1))
    res$SGX <- rowSums(geno)
    res$nu <- (0.5 * res$X1) + (0.5 * res$X2) + (res$gamma * res$SGX)
    res$Tx <- res$lambda * sqrt(-log(res$U) / exp(res$nu))
    res <- res %>%
      mutate(time = ifelse(C <= Tx, C, Tx),
             event = ifelse(C <= Tx, 0, 1))
    res[, c("SGX", "nu", "U", "lambda", "C", "Tx")] <- NULL
    
    data <- rbind(data, res)
  }
  
  filename <- paste0("MAF", MAF, "_Lambda", lambdak[l, 2])
  filename <- sub("\\.", "_", filename)
  filename <- sub("\\.", "_", filename)
  saveRDS(data, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, ".rds"),
          compress = "xz")
  saveRDS(geno, file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/",
                              filename, "_geno.rds"),
          compress = "xz")
  
  stop <- Sys.time()
  print(stop - start)
  
}