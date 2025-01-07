## LR
library(survival)
library(parallel)

rm(list = ls())
gc()

setwd("/srv/__PERSONAL__/Szakdoga/Szakdoga/mintaszam")


data2 <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5_geno.rds")
data <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5.rds")

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1


k = 11

run <- data[rows1[k]:rows[k + 1], ]

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

pek2 <- c()

mintaszam_lista <- c(1000, 5000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 7000080000, 90000)


num_cores <- 8


start <- Sys.time()

for(mintaszam in mintaszam_lista){
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  
  for(i in 1:50){
    sample <- sample(1:nrow(run), mintaszam, replace = F)
    
    run_seq <- run[sample, ]
    run_seq$ID <- paste0("S", 1:mintaszam)
    
    geno <- data2[sample, ]
    rownames(geno) <- run_seq$ID
    
    run_seq$X2 <- as.factor(run_seq$X2)
    
    fit0 <- coxph(Surv(time, event) ~ X1 + X2,  
                  data = run_seq)
    
    pek <- mclapply(colnames(geno), function(g) {
      fit <- coxph(as.formula(paste0("Surv(time, event) ~ X1 + X2 + ", g)),  
                   data = cbind(run_seq, geno[g]))
      p_values <- anova(fit0, fit, test = "Chisq")$`Pr(>|Chi|)`[2]
      return(p_values)
    }, mc.cores = num_cores)  
    
    
    result[i, "p"] <- sum(pek > 0.00005)
    result[i, "p5"] <- sum(pek <= 0.00005)
    result[i, "p8"] <- sum(pek <= 0.00000005)
    result[i, "r"] <- i
    pek2 <- c(pek2, pek)
    
    
  }
  
  result$mintaszam <- mintaszam
  plot_result <- rbind(plot_result, result)
  
  save(plot_result, pek2, file = "LR_mintaszam3.RData")
}

stop <- Sys.time()
print(stop - start)


plot_result2 <- plot_result


##################################################
## gwasurvivr
library(gwasurvivr)
library(parallel)
library(gaston)

rm(list = ls()[!ls() %in% c("data", "data2")])

gc()

setwd("/srv/__PERSONAL__/Szakdoga/Szakdoga/mintaszam")


data2 <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5_geno.rds")
data <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5.rds")

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1


k = 11

run <- data[rows1[k]:rows[k + 1], ]

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

pek2 <- c()

mintaszam_lista <- c(1000, 5000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 7000080000, 90000)


geno_info <- readRDS("geno_info.rds")
bim <- data.frame(chr = as.character(geno_info$chr),
                  id = geno_info$id,
                  dist = 0,
                  pos = geno_info$pos,
                  A1 = geno_info$A1,
                  A2 = geno_info$A2)

start <- Sys.time()

for(mintaszam in mintaszam_lista){
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  
  pek <- c()
  
  for(i in 1:50){
    sample <- sample(1:nrow(run), mintaszam, replace = F)
    
    run_seq <- run[sample, ]
    run_seq$ID <- paste0("S", 1:mintaszam)
    
    geno <- data2[sample, ]
    geno <- as.matrix(geno)
    rownames(geno) <- run_seq$ID
    
    
    fam <- data.frame(famid = NA,
                      id = run_seq$ID,
                      father = NA,
                      mother = NA,
                      sex = NA,
                      pheno = NA)
    
    bed <- gaston::as.bed.matrix(geno, fam, bim)
    write.bed.matrix(bed,
                     basename = "mintaszam")
    
    cl <- makeForkCluster(nnodes = getOption("mc.cores", 8))
    
    gwasurvivr::plinkCoxSurv(bed.file = "mintaszam.bed",
                             covariate.file = run_seq,
                             id.column = "ID",
                             sample.ids = run_seq$ID,
                             time.to.event = "time",
                             event = "event",
                             covariates = c("X1", "X2"),
                             print.covs = "only",
                             out.file = "gwasurvivr",
                             chunk.size = 20,
                             maf.filter = NULL,
                             clusterObj = cl)
    
    gwas <- data.table::fread(paste0("/srv/__PERSONAL__/Szakdoga/Szakdoga/mintaszam/gwasurvivr.coxph"))
    
    result[i, "p"] <- length(gwas$PVALUE[gwas$PVALUE > 0.00005])
    result[i, "p5"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00005])
    result[i, "p8"] <- length(gwas$PVALUE[gwas$PVALUE <= 0.00000005])
    result[i, "r"] <- i
    pek <- c(pek, gwas$PVALUE)
    
    
  }
  
  result$mintaszam <- mintaszam
  plot_result <- rbind(plot_result, result)
  
  pek2 <- c(pek2, pek)
  
  save(plot_result, pek2, file = "gwasurvivr_mintaszam3.RData")
}

stop <- Sys.time()
print(stop - start)

########################################################
## adult
library(LTFHPlus)
library(dplyr)
library(MASS)
library(future.apply)
library(tidyr)
library(survival)
library(bigstatsr)

rm(list = ls()[!ls() %in% c("data", "data2")])
gc()

setwd("/srv/__PERSONAL__/Szakdoga/Szakdoga/mintaszam")


data2 <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5_geno.rds")
data <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5.rds")

h2 <- 0.5

nthreads <- 16
plan(tweak(multisession, workers = nthreads))

data$time <- round(data$time, 7)
data$aoo <- ifelse(data$event == 1, data$time, NA)
data$role <- "g"
data$fid <- as.character(1:nrow(data))
data$pid <- as.character(1:nrow(data))

covariates <- as.matrix(data[, c("X1", "X2")])

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1


k = 11

run <- data[rows1[k]:rows[k + 1], ]
cov <- covariates[rows1[k]:rows[k + 1], ]

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

pek2 <- c()

mintaszam_lista <- c(1000, 5000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 7000080000, 90000)

start <- Sys.time()

for(mintaszam in mintaszam_lista){
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  
  pek <- c()
  
  for(i in 1:50){
    sample <- sample(1:nrow(run), mintaszam, replace = F)
    
    run_seq <- run[sample, ]
    cov_seq <- cov[sample, ]
    
    # CIP meghatározás
    run_seq$seged <- 1
    
    cip <- run_seq %>%
      arrange(time) %>% 
      group_by(time) %>%
      summarise(event = sum(event),
                seged = sum(seged))
    
    cip$cip <- cip$event / nrow(run_seq)
    cip$cip <- cumsum(cip$cip)
    cip$cip <- round(cip$cip, 7)
    
    cip$cip <- ifelse(cip$cip == 0,
                      0,
                      ifelse(cip$cip < 2e-04, 2e-04, cip$cip))
    
    run_seq$seged <- NULL
    cip[, c("event", "seged")] <- NULL
    
    sajat_feno <- prepare_LTFHPlus_input(.tbl = run_seq,
                                         CIP = cip, 
                                         age_col = "time",
                                         aoo_col = "aoo",
                                         CIP_merge_columns = c("time"),
                                         CIP_cip_col = "cip",
                                         status_col = "event",
                                         use_fixed_case_thr = F,
                                         fam_id_col = "fid",
                                         personal_id_col = "pid",
                                         interpolation = NA,
                                         min_CIP_value = 1e-5)
    
    sajat_feno <- unique(sajat_feno)
    
    
    sajat_liab <- estimate_liability(.tbl = sajat_feno,
                                     h2 = h2,
                                     pid = "pid",
                                     fam_id = "fid",
                                     role = "role",
                                     tol = 0.1,
                                     out = c(1))
    
    geno <- data2[sample, ]
    geno <- as_FBM(geno)
    
    ########
    gwas <- big_univLinReg(X = geno,
                           y.train = sajat_liab$genetic_est,
                           ncores = 1,
                           covar.train = cov_seq)
    
    gwas$p.value <- predict(gwas, log10 = FALSE)
    
    
    
    
    
    result[i, "p"] <- sum(gwas$p.value > 0.00005)
    result[i, "p5"] <- sum(gwas$p.value <= 0.00005)
    result[i, "p8"] <- sum(gwas$p.value <= 0.00000005)
    result[i, "r"] <- i
    pek <- c(pek, gwas$p.value)
    
    
  }
  
  result$mintaszam <- mintaszam
  plot_result <- rbind(plot_result, result)
  
  pek2 <- c(pek2, pek)
  
  save(plot_result, pek2, file = "adult_mintaszam.RData")
}

stop <- Sys.time()
print(stop - start)

###########################################################
## SPACox
library(SPACox)
library(survival)

rm(list = ls()[!ls() %in% c("data", "data2")])

gc()

setwd("/srv/__PERSONAL__/Szakdoga/Szakdoga/mintaszam")


data2 <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5_geno.rds")
data <- readRDS("~/__PERSONAL__/Szakdoga/Szakdoga/Forras/MAF0_3_Lambda0_5.rds")

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1

k = 11

rows <- seq(0, nrow(data), by = 100000*50)
rows1 <- rows + 1


run <- data[rows1[k]:rows[k + 1], ]

plot_result <- data.frame(p = numeric(0),
                          p5 = numeric(0),
                          p8 = numeric(0),
                          r = numeric(0),
                          point = numeric(0))

pek <- c()

mintaszam_lista <- c(1000, 5000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 7000080000, 90000)

start <- Sys.time()

for(mintaszam in mintaszam_lista){
  
  result <- data.frame(p = numeric(0),
                       p5 = numeric(0),
                       p8 = numeric(0),
                       r = numeric(0))
  
  for(i in 1:50){
    sample <- sample(1:nrow(run), mintaszam, replace = F)
    
    run_seq <- run[sample, ]
    run_seq$ID <- paste0("S", 1:mintaszam)
    
    geno <- data2[sample, ]
    geno <- as.matrix(geno)
    rownames(geno) <- run_seq$ID
    
    run_seq$X2 <- as.factor(run_seq$X2)
    obj.null <- SPACox_Null_Model(Surv(time, event) ~ X1 + X2,
                                  data = run_seq,
                                  pIDs = run_seq$ID,
                                  gIDs = rownames(geno),
                                  length.out = mintaszam)
    SPACox.res <- SPACox(obj.null = obj.null, Geno.mtx = geno, G.model = "Add")
    SPACox.res <- as.data.frame(SPACox.res)
    pek <- c(pek, SPACox.res$p.value.spa)
    result[i, "p"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa > 0.00005])
    result[i, "p5"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00005])
    result[i, "p8"] <- length(SPACox.res$p.value.spa[SPACox.res$p.value.spa <= 0.00000005])
    result[i, "r"] <- i
    
  }
  
  result$mintaszam <- mintaszam
  plot_result <- rbind(plot_result, result)
  
}

stop <- Sys.time()
print(stop - start)

save(plot_result, pek, file = "spcox_mintaszam3.RData")


