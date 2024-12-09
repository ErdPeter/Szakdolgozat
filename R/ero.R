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

#dat <- readRDS("MAF0_05_Lambda0_01.rds")


###############################################################################
###############################################################################
###############################################################################
## Geno info
set.seed(1234)
mut <- function(x){
  list <- setdiff(c("A", "T", "G", "C"), x)
  return(sample(list, 1))
  
}

geno_info <- data.frame(id = paste0("G", 1:20),
                        chr = sample(as.character(1:22), 20, replace = T),
                        dist = 0,
                        pos = round(runif(20, 1, 99999999), 0),
                        A1 = sample(c("A", "T", "G", "C"), 20, replace = T))
geno_info$A2 <- unlist(lapply(geno_info$A1, mut))

saveRDS(geno_info,
        file = paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/geno_info.rds"),
        compress = "xz")


## Bed filek
set.seed(1234)
setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")

rm(list = ls())
gc()

geno_info <- readRDS("geno_info.rds")

fam <- data.frame(famid = NA,
                  id = paste0("S", 1:100000),
                  father = NA,
                  mother = NA,
                  sex = NA,
                  pheno = NA)

bim <- data.frame(chr = as.character(geno_info$chr),
                  id = geno_info$id,
                  dist = 0,
                  pos = geno_info$pos,
                  A1 = geno_info$A1,
                  A2 = geno_info$A2)

files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "_geno")

maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)

for(g in 1:15){
data2 <- readRDS(files[g])
print(paste0(g, ". beolvasas kész!"))

maf_act <- maf[g]
lambda_act <- lambda[g]

  for(i in 1:50){
    if(i%%10 == 0){print(i)}
    start_index <- (i - 1) * 100000 + 1
    end_index <- min(i * 100000, nrow(data2))
    
    geno <- as.matrix(data2[start_index:end_index,])
    rownames(geno) <- paste0("S", 1:100000)
  
    bed <- gaston::as.bed.matrix(geno, fam, bim)
    write.bed.matrix(bed,
                     basename = paste("c:/Tanulas/Szakdolgozat/Ero2/Data/BED/bed_maf",
                                      maf_act,
                                      "lambda",
                                      lambda_act,
                                      i,
                                      sep = "_"))
  }
}


## GWASSURVIVR
library(gwasurvivr)

setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")

rm(list = ls())
gc()

maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)

bedek <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data/BED",
                    full.names = F,
                    pattern = ".bed")
bedek <- sort(bedek)
bedek <- paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/BED/", bedek)

files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

files <- files[!grepl("geno", files)]


for(h in 13:15){
filesorszam <- h
data <- readRDS(paste0("c:/Tanulas/Szakdolgozat/Ero2/Data/", files[filesorszam]))

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

for(k in 1:11){
start <- Sys.time()

gc()
run <- data[rows1[k]:rows[k + 1], ]

result <- data.frame(p = numeric(0),
                     p5 = numeric(0),
                     p8 = numeric(0),
                     r = numeric(0))

bedek_tmp <- bedek[((h - 1) * 50 + 1):(min(h * 50, nrow(bedek)))]
bedek_tmp <- bedek_tmp[c(1, 12, 23, 34, 45, 47:50, 2:11, 13:22, 24:33, 35:44, 46)]

for(i in 1:50){
start_index <- (i - 1) * 100000 + 1
end_index <- min(i * 100000, nrow(run))
run_seq <- run[start_index:end_index, ]
run_seq$ID <- paste0("S", 1:100000)
bedfile <- bedek_tmp[i]

gwasurvivr::plinkCoxSurv(bed.file = bedfile,
                         covariate.file = run_seq[, c(1:2, 4:6)],
                         id.column = "ID",
                         sample.ids = run_seq$ID,
                         time.to.event = "time",
                         event = "event",
                         covariates = c("X1", "X2"),
                         print.covs = "only",
                         out.file = paste("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2/GWAS_result", maf[h], lambda[h], k, i, sep = "_"),
                         chunk.size = 20,
                         maf.filter = NULL)

gwas <- data.table::fread(paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2/GWAS_result_", maf[h], "_", lambda[h], "_", k, "_", i, ".coxph"))
result[i, "p"] <- length(gwas$PVALUE[gwas$PVALUE > 0.00005])
result[i, "p5"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00005])
result[i, "p8"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00000005])
result[i, "r"] <- i

print(paste0("diagram: ", h, "., pont: ", k, "., ciklus szám: ", i))

}

result$point <- k
plot_result <- rbind(plot_result, result)

save(result, file = "c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2_res/TEMP.RData")

stop <- Sys.time()
print(stop - start)

}

rm(filename)
filename <- files[filesorszam]
filename <- gsub(".rds", "", filename)
filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr2_res/", filename, "_gwasurvivr.rds")
saveRDS(plot_result, file = filename, compress = "xz")
}


###############################################################################
###############################################################################
###############################################################################
## SPACox
library(SPACox)
library(survival)

rm(list = ls())
gc()

setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")
files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

geno_files <- files[grepl("geno", files)]
files <- files[!grepl("geno", files)]

for(h in 1:15){
data <- readRDS(files[h])
data2 <- readRDS(geno_files[h])

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

for(k in 1:11){
  start <- Sys.time()
  
  gc()
  run <- data[rows1[k]:rows[k + 1], ]
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  
  for(i in 1:50){
    start_index <- (i - 1) * 100000 + 1
    end_index <- min(i * 100000, nrow(run))
    run_seq <- run[start_index:end_index, ]
    run_seq$ID <- paste0("S", 1:100000)
    
    geno <- data2[start_index:end_index, ]
    geno <- as.matrix(geno)
    rownames(geno) <- run_seq$ID

    run_seq$X2 <- as.factor(run_seq$X2)
    obj.null <- SPACox_Null_Model(Surv(time, event) ~ X1 + X2,
                                  data = run_seq,
                                  pIDs = run_seq$ID,
                                  gIDs = rownames(geno),
                                  length.out = 10000)
    SPACox.res <- SPACox(obj.null = obj.null, Geno.mtx = geno, G.model = "Add")
    SPACox.res <- as.data.frame(SPACox.res)
    result[i, "p"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa > 0.00005])
    result[i, "p5"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00005])
    result[i, "p8"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00000005])
    result[i, "r"] <- i
    
  }
  
  result$point <- k
  plot_result <- rbind(plot_result, result)
  
  stop <- Sys.time()
  print(stop - start)
  
}

filename <- files[h]
filename <- gsub(".rds", "", filename)
filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/spacox/", filename, "_spacox.rds")
saveRDS(plot_result, file = filename, compress = "xz")
}

###############################################################################
###############################################################################
###############################################################################
## LR teszt
library(survival)

rm(list = ls())
gc()

setwd("c:/Tanulas/Szakdolgozat/Ero2/Data")
files <- list.files("c:/Tanulas/Szakdolgozat/Ero2/Data",
                    full.names = F,
                    pattern = "MAF")

geno_files <- files[grepl("geno", files)]
files <- files[!grepl("geno", files)]

for(h in 1:15){
  data <- readRDS(files[h])
  data2 <- readRDS(geno_files[h])
  
  rows <- seq(0, nrow(data), by = 100000*50)
  rows1 <- rows + 1
  
  plot_result <- data.frame(p = numeric(0),
                            p5 = numeric(0),
                            p8 = numeric(0),
                            r = numeric(0),
                            point = numeric(0))
  
  for(k in 1:11){
    start <- Sys.time()
    
    gc()
    run <- data[rows1[k]:rows[k + 1], ]
    
    result <- data.frame(p = numeric(0),
                         p5 = numeric(0),
                         p8 = numeric(0),
                         r = numeric(0))
    
    for(i in 1:50){
      start_index <- (i - 1) * 100000 + 1
      end_index <- min(i * 100000, nrow(run))
      run_seq <- run[start_index:end_index, ]
      run_seq$ID <- paste0("S", 1:100000)
      
      geno <- data2[start_index:end_index, ]
      rownames(geno) <- run_seq$ID
      
      run_seq$X2 <- as.factor(run_seq$X2)
      
      fit0 <- coxph(Surv(time, event) ~ X1 + X2,  
                    data = run_seq)
      pek <- c()
      for(g in 1:ncol(geno)){
        gcol <- colnames(geno)[g]
        fmla <- as.formula(paste0("Surv(time, event) ~ X1 + X2 + ", gcol))
        run_seq_tmp <- cbind(run_seq, geno[gcol])
        fit <- coxph(fmla,  
                      data = run_seq_tmp)
        pek[g] <- anova(fit0, fit, test = "Chisq")$`Pr(>|Chi|)`[2]
      }
      
      result[i, "p"] <- sum(pek > 0.00005)
      result[i, "p5"] <- sum(pek <= 0.00005)
      result[i, "p8"] <- sum(pek <= 0.00000005)
      result[i, "r"] <- i
      
    }
    
    result$point <- k
    plot_result <- rbind(plot_result, result)
    
    stop <- Sys.time()
    print(stop - start)
    
  }
  
  filename <- files[h]
  filename <- gsub(".rds", "", filename)
  filename <- paste0("c:/Tanulas/Szakdolgozat/Ero2/ResData/Score/", filename, "_score.rds")
  saveRDS(plot_result, file = filename, compress = "xz")
}



###############################################################################
###############################################################################
###############################################################################
## Plot gwasurvivr
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("c:/Tanulas/Szakdolgozat/Ero2/ResData/gwasurvivr")

files <- list.files(pattern = ".rds")

gammak_l <- list(seq(from = 0, to = 3, length.out = 11),
                 seq(from = 0, to = 1.5, length.out = 11),
                 seq(from = 0, to = 0.6, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.25, length.out = 11),
                 seq(from = 0, to = 2, length.out = 11),
                 seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.3, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11),
                 seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.5, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11),
                 seq(from = 0, to = 0.1, length.out = 11),
                 seq(from = 0, to = 0.06, length.out = 11))


maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)

for(num in 1:15){

dat <- readRDS(files[num])
dat$r <- NULL
dat2 <- dat %>% 
  group_by(point) %>% 
  summarise(across(everything(), mean))  %>%
  as.data.frame()
dat2$p52 <- dat2$p5 / 20
dat2$p82 <- dat2$p8 / 20
dat2$gamma <- gammak_l[[num]]

p <- ggplot(dat2, aes(gamma, p82)) +
  geom_point() +
  geom_line() +
  labs (x = "Genetikai hatás",
        y = "Empirikus erő",
        title = paste0("p <= 5*10e-8, MAF = ", maf[num], " ER = ", lambda[num], "%"))
assign(paste0("p", num), p)

pp <- ggplot(dat2, aes(gamma, p52)) +
  geom_point() +
  geom_line() +
  labs (x = "Genetikai hatás",
        y = "Empirikus erő",
        title = paste0("p <= 5*10e-5, MAF = ", maf[num], " ER = ", lambda[num], "%"))
assign(paste0("pp", num), pp)
}


jpeg(file="GWASSURVIVR_p8.jpeg")
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15
dev.off()

jpeg(file="GWASSURVIVR_p5.jpeg")
pp1+pp2+pp3+pp4+pp5+pp6+pp7+pp8+pp9+pp10+pp11+pp12+pp13+pp14+pp15
dev.off()

## Plot SPACox
library(dplyr)
library(ggplot2)
library(patchwork)

setwd("c:/Tanulas/Szakdolgozat/Ero2/ResData/spacox")

files <- list.files(pattern = ".rds")

gammak_l <- list(seq(from = 0, to = 3, length.out = 11),
                 seq(from = 0, to = 1.5, length.out = 11),
                 seq(from = 0, to = 0.6, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.25, length.out = 11),
                 seq(from = 0, to = 2, length.out = 11),
                 seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.4, length.out = 11),
                 seq(from = 0, to = 0.3, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11),
                 seq(from = 0, to = 1, length.out = 11),
                 seq(from = 0, to = 0.5, length.out = 11),
                 seq(from = 0, to = 0.15, length.out = 11),
                 seq(from = 0, to = 0.1, length.out = 11),
                 seq(from = 0, to = 0.06, length.out = 11))


maf <- c(rep(0.01, 5), rep(0.05, 5), rep(0.3, 5))
lambda <- rep(c(0.002, 0.01, 0.1, 0.2, 0.5), 3)



for(num in 1:15){
  
  dat <- readRDS(files[num])
  dat$r <- NULL
  dat2 <- dat %>% 
    group_by(point) %>% 
    summarise(across(everything(), mean))  %>%
    as.data.frame()
  dat2$p52 <- dat2$p5 / 20
  dat2$p82 <- dat2$p8 / 20
  dat2$gamma <- gammak_l[[num]]
  
  p <- ggplot(dat2, aes(gamma, p82)) +
    geom_point() +
    geom_line() +
    labs (x = "Genetikai hatás",
          y = "Empirikus erő",
          title = paste0("p <= 5*10e-8, MAF = ", maf[num], " ER = ", lambda[num], "%"))
  assign(paste0("p", num), p)
  
  pp <- ggplot(dat2, aes(gamma, p52)) +
    geom_point() +
    geom_line() +
    labs (x = "Genetikai hatás",
          y = "Empirikus erő",
          title = paste0("p <= 5*10e-5, MAF = ", maf[num], " ER = ", lambda[num], "%"))
  assign(paste0("pp", num), pp)
}

jpeg(file="SPACox_p8.jpeg")
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15
dev.off()

jpeg(file="SPACox_p5.jpeg")
pp1+pp2+pp3+pp4+pp5+pp6+pp7+pp8+pp9+pp10+pp11+pp12+pp13+pp14+pp15
dev.off()









