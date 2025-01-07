library(gwasurvivr)
library(LTFHPlus)
library(dplyr)
library(Rlab)
library(SPACox)
library(survival)
library(ggplot2)
library(gaston)
library(MASS)
library(future.apply)
library(tidyr)
library(bigstatsr)
library(parallel)



setwd("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/")

### Tesztadat szimulálása
set.seed(1234)

N = 10000 # mintaszám
MAF = 0.01 # MAF
geno_alt = 100000
ER = 0.1

geno <- matrix(rbinom(N * geno_alt * 2, size = 2, prob = MAF),
               nrow = N * 2,
               ncol = geno_alt)
#geno <- as.data.frame(geno)
colnames(geno) <- paste0("G", 1:ncol(geno))

saveRDS(geno, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/teljesitmeny_geno.RDS")


feno <- data.frame(Time = runif(N * 2),
                   Event = rbinom(N * 2, size = 1, prob = ER),
                   X1 = rnorm(N * 2),
                   X2 = rbern(N * 2, 0.5),
                   ID = paste0("S", 1:(N * 2)))

saveRDS(feno, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/teljesitmeny_feno.RDS")


geno <- readRDS(file = "teljesitmeny_geno.RDS")
feno <- readRDS(file = "teljesitmeny_feno.RDS")
rownames(geno) <- feno$ID
feno$X2 <- as.factor(feno$X2)

###############################################################################
## SPACox
set.seed(1234)
#geno <- as.matrix(geno)


n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 

SPAcox_time_teszt <- data.frame(n = numeric(0),
                                ido = numeric(0))



for (i in 1:50){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- as.matrix(geno[samplek, ])
  
  start <- Sys.time()
  obj.null <- SPACox_Null_Model(Surv(Time, Event) ~ X1 + X2,
                                data = feno_tmp,
                                pIDs = samplek,
                                gIDs = samplek,
                                length.out = 10000)
  SPACox.res <- SPACox(obj.null = obj.null,
                       Geno.mtx = geno_tmp,
                       Cutoff = 2,
                       impute.method = "fixed",
                       missing.cutoff = 0.15,
                       min.maf = 1e-04,
                       CovAdj.cutoff = 5e-05)
  stop <- Sys.time()
  SPAcox_time[i, "n"] <- db
  SPAcox_time[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  rm(stop, start, SPACox.res, obj.null, geno_tmp, feno_tmp, db, samplek)
  print(i)
}

save(SPAcox_time, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/SPAcox_time_origi.RData")

SPAcox_time[SPAcox_time$n == 10000, "ido"] <- SPAcox_time[SPAcox_time$n == 10000, "ido"] * 60
SPAcox_time[SPAcox_time$n == 10000, "ido"] <- SPAcox_time[SPAcox_time$n == 10000, "ido"] * 60

SPAcox_time_p <- SPAcox_time

summary_SPAcox_time <- SPAcox_time_p %>%
  group_by(n) %>%
  summarise(
    mean_ido = mean(ido),
    sd_ido = sd(ido)
  )

g_SPAcox_time <- ggplot(summary_SPAcox_time, aes(x = n, y = mean_ido)) +
  geom_line() +  # Görbék rajzolása
  geom_point() +  # Pontok rajzolása
  geom_errorbar(aes(ymin = mean_ido - sd_ido, ymax = mean_ido + sd_ido), width = 0.2) +  
  labs(
    x = "Minta szám", 
    y = "Idő" 
  ) +
  theme_minimal() +  # Egyszerű téma
  ggtitle("SPAcox idő")

save(SPAcox_time_p, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/SPAcox_time_p.RData")

#load(file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny/SPAcox_time_mp.RData")


###############################################################################
## GWASSURVIVR

## Geno info
set.seed(1234)
mut <- function(x){
  list <- setdiff(c("A", "T", "G", "C"), x)
  return(sample(list, 1))
  
}

geno_info <- data.frame(id = paste0("G", 1:geno_alt),
                        chr = sample(as.character(1:22), 20, replace = T),
                        dist = 0,
                        pos = round(runif(geno_alt, 1, 99999999), 0),
                        A1 = sample(c("A", "T", "G", "C"), geno_alt, replace = T))
geno_info$A2 <- unlist(lapply(geno_info$A1, mut))

saveRDS(geno_info,
        file = paste0("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/geno_info.rds"),
        compress = "xz")

#geno_info <- readRDS(file = paste0("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny/geno_info.rds"))


## Bed filek
set.seed(1234)

geno_info <- readRDS("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/geno_info.rds")

n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 

bim <- data.frame(chr = as.character(geno_info$chr),
                  id = geno_info$id,
                  dist = 0,
                  pos = geno_info$pos,
                  A1 = geno_info$A1,
                  A2 = geno_info$A2)


for(i in 1:50){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- as.matrix(geno[samplek, ])
  
  fam <- data.frame(famid = NA,
                    id = feno_tmp$ID,
                    father = NA,
                    mother = NA,
                    sex = NA,
                    pheno = NA)
  
  bed <- gaston::as.bed.matrix(geno_tmp, fam, bim)
  write.bed.matrix(bed,
                   basename = paste("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/bed_n",
                                    db,
                                    "i",
                                    i,
                                    sep = "_"))
}


# Idő mérés

bedek <- list.files("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2",
                    full.names = T,
                    pattern = ".bed")
bedek <- sort(bedek)
bedek <- bedek[c(1, 3:10, 2, 31:40, 11:20, 41:50, 21:30)]

famok <- list.files("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2",
                    full.names = T,
                    pattern = "fam")
famok <- sort(famok)
famok <- famok[c(1, 3:10, 2, 31:40, 11:20, 41:50, 21:30)]

bimek <- list.files("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2",
                    full.names = T,
                    pattern = "bim")
bimek <- sort(bimek)
bimek <- bimek[c(1, 3:10, 2, 31:40, 11:20, 41:50, 21:30)]

feno$X2 <- as.numeric(feno$X2)

GWASSURVIVR_time <- data.frame(n = numeric(0),
                               ido = numeric(0))


for(i in 1:50){
  db <- n[i]
  bedfile <- bedek[i]
  famfile <- famok[i]
  bimfile <- bimek[i]
  basename <- gsub("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/", "", bedfile)
  basename <- gsub(".bed", "", basename)
  
  
  bed <- read.bed.matrix(basename = basename,
                         bed = bedfile, 
                         fam = famfile,
                         bim = bimfile
  )
  
  samplek <- bed@ped$id
  feno_tmp <- feno[feno$ID %in% samplek, ]
  
  cl <- makeForkCluster(nnodes = getOption("mc.cores", 1))
  
  start <- Sys.time()
  gwasurvivr::plinkCoxSurv(bed.file = bedfile,
                           covariate.file = feno_tmp,
                           id.column = "ID",
                           sample.ids = feno_tmp$ID,
                           time.to.event = "Time",
                           event = "Event",
                           covariates = c("X1", "X2"),
                           print.covs = "only",
                           out.file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/GWAS_result",
                           chunk.size = 1000,
                           maf.filter = NULL,
                           clusterObj = NULL)
  
  stop <- Sys.time()
  GWASSURVIVR_time[i, "n"] <- db
  GWASSURVIVR_time[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  save(GWASSURVIVR_time, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/GWASSURVIVR_time_origi.RData")
  print(i)
  rm(db, bedfile, famfile, bimfile, basename, bed, samplek, feno_tmp, start, stop)
  
}


GWASSURVIVR_time[GWASSURVIVR_time$n == 10000, "ido"] <- GWASSURVIVR_time[GWASSURVIVR_time$n == 10000, "ido"] * 60

summary_GWASSURVIVR_time <- GWASSURVIVR_time %>%
  group_by(n) %>%
  summarise(
    mean_ido = mean(ido),
    sd_ido = sd(ido)
  )

g_GWASSURVIVR_time <- ggplot(summary_GWASSURVIVR_time, aes(x = n, y = mean_ido)) +
  geom_line() +  # Görbék rajzolása
  geom_point() +  # Pontok rajzolása
  geom_errorbar(aes(ymin = mean_ido - sd_ido, ymax = mean_ido + sd_ido), width = 0.2) +  
  labs(
    x = "Minta szám", 
    y = "Idő" 
  ) +
  theme_minimal() +  # Egyszerű téma
  ggtitle("SPAcox idő")

###############################################################################
## gwasurvivr 4

GWASSURVIVR_time4 <- data.frame(n = numeric(0),
                                ido = numeric(0))



for(i in 1:50){
  db <- n[i]
  bedfile <- bedek[i]
  famfile <- famok[i]
  bimfile <- bimek[i]
  basename <- gsub("/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/", "", bedfile)
  basename <- gsub(".bed", "", basename)
  
  
  bed <- read.bed.matrix(basename = basename,
                         bed = bedfile, 
                         fam = famfile,
                         bim = bimfile
  )
  
  samplek <- bed@ped$id
  feno_tmp <- feno[feno$ID %in% samplek, ]
  
  
  cl <- makeForkCluster(nnodes = getOption("mc.cores", 4))
  
  start <- Sys.time()
  gwasurvivr::plinkCoxSurv(bed.file = bedfile,
                           covariate.file = feno_tmp,
                           id.column = "ID",
                           sample.ids = feno_tmp$ID,
                           time.to.event = "Time",
                           event = "Event",
                           covariates = c("X1", "X2"),
                           print.covs = "only",
                           out.file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/GWAS_result",
                           chunk.size = 1000,
                           maf.filter = NULL,
                           clusterObj = cl)
  
  stop <- Sys.time()
  GWASSURVIVR_time4[i, "n"] <- db
  GWASSURVIVR_time4[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  save(GWASSURVIVR_time4, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/GWASSURVIVR_time_origi4_b.RData")
  print(i)
  rm(db, bedfile, famfile, bimfile, basename, bed, samplek, feno_tmp, start, stop)
  
}



summary_GWASSURVIVR_time <- GWASSURVIVR_time %>%
  group_by(n) %>%
  summarise(
    mean_ido = mean(ido),
    sd_ido = sd(ido)
  )

g_GWASSURVIVR_time <- ggplot(summary_GWASSURVIVR_time, aes(x = n, y = mean_ido)) +
  geom_line() +  # Görbék rajzolása
  geom_point() +  # Pontok rajzolása
  geom_errorbar(aes(ymin = mean_ido - sd_ido, ymax = mean_ido + sd_ido), width = 0.2) +  
  labs(
    x = "Minta szám", 
    y = "Idő" 
  ) +
  theme_minimal() +  # Egyszerű téma
  ggtitle("GWASurvivr idő")


###############################################################################
## Adult 1

h2 <- 0.5

nthreads <- 1
plan(tweak(multisession, workers = nthreads))

#data$time <- round(data$time, 7)
feno$aoo <- ifelse(feno$Event == 1, feno$Time, NA)
feno$role <- "g"
feno$fid <- as.character(1:nrow(feno))
feno$pid <- as.character(1:nrow(feno))

feno$X2 <- as.numeric(feno$X2)
covariates <- feno[, c("X1", "X2", "ID")]

Adult_time <- data.frame(n = numeric(0),
                         ido = numeric(0))

n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 

for (i in 1:length(n)){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- geno[samplek, ]
  cov_tmp <- covariates[covariates$ID %in% samplek, ]
  cov_tmp$ID <- NULL
  cov_tmp <- as.matrix(cov_tmp)
  
  
  feno_tmp$seged <- 1
  
  cip <- feno_tmp %>%
    arrange(Time) %>% 
    group_by(Time) %>%
    summarise(Event = sum(Event),
              seged = sum(seged))
  
  cip$cip <- cip$Event / nrow(feno_tmp)
  cip$cip <- cumsum(cip$cip)
  cip$cip <- round(cip$cip, 7)
  
  cip$cip <- ifelse(cip$cip == 0,
                    0,
                    ifelse(cip$cip < 2e-04, 2e-04, cip$cip))
  
  feno_tmp$seged <- NULL
  cip[, c("Event", "seged")] <- NULL
  
  sajat_feno <- prepare_LTFHPlus_input(.tbl = feno_tmp,
                                       CIP = cip, 
                                       age_col = "Time",
                                       aoo_col = "aoo",
                                       CIP_merge_columns = c("Time"),
                                       CIP_cip_col = "cip",
                                       status_col = "Event",
                                       use_fixed_case_thr = F,
                                       fam_id_col = "fid",
                                       personal_id_col = "pid",
                                       interpolation = NA,
                                       min_CIP_value = 1e-5)
  
  sajat_feno <- unique(sajat_feno)
  
  start <- Sys.time()
  sajat_liab <- estimate_liability(.tbl = sajat_feno,
                                   h2 = h2,
                                   pid = "pid",
                                   fam_id = "fid",
                                   role = "role",
                                   tol = 0.1,
                                   out = c(1))
  
  
  
  geno_tmp <- as_FBM(geno_tmp)
  gwas <- big_univLinReg(X = geno_tmp,
                         y.train = sajat_liab$genetic_est,
                         ncores = 1,
                         covar.train = cov_tmp)
  
  stop <- Sys.time()
  Adult_time[i, "n"] <- db
  Adult_time[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  save(Adult_time, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/Adult_time_origi.RData")
  print(i)
  rm(stop, start, cip, obj.null, geno_tmp, feno_tmp, db, samplek, cov_tmp, sajat_feno, sajat_liab, gwas)
  
  
  
}

#load(file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/Adult_time_origi.RData")


summary_Adult_time <- Adult_time %>%
  group_by(n) %>%
  summarise(
    mean_ido = mean(ido),
    sd_ido = sd(ido)
  )

g_Adult_time <- ggplot(summary_Adult_time, aes(x = n, y = mean_ido)) +
  geom_line() +  # Görbék rajzolása
  geom_point() +  # Pontok rajzolása
  geom_errorbar(aes(ymin = mean_ido - sd_ido, ymax = mean_ido + sd_ido), width = 0.2) +  
  labs(
    x = "Minta szám", 
    y = "Idő" 
  ) +
  theme_minimal() +  # Egyszerű téma
  ggtitle("SPAcox idő")



###############################################################################
## Adult 4

h2 <- 0.5

nthreads <- 4
plan(tweak(multisession, workers = nthreads))

#data$time <- round(data$time, 7)
feno$aoo <- ifelse(feno$Event == 1, feno$Time, NA)
feno$role <- "g"
feno$fid <- as.character(1:nrow(feno))
feno$pid <- as.character(1:nrow(feno))

feno$X2 <- as.numeric(feno$X2)
covariates <- feno[, c("X1", "X2", "ID")]

Adult_time4 <- data.frame(n = numeric(0),
                          ido = numeric(0))

n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 

for (i in 1:length(n)){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- geno[samplek, ]
  cov_tmp <- covariates[covariates$ID %in% samplek, ]
  cov_tmp$ID <- NULL
  cov_tmp <- as.matrix(cov_tmp)
  
  
  feno_tmp$seged <- 1
  
  cip <- feno_tmp %>%
    arrange(Time) %>% 
    group_by(Time) %>%
    summarise(Event = sum(Event),
              seged = sum(seged))
  
  cip$cip <- cip$Event / nrow(feno_tmp)
  cip$cip <- cumsum(cip$cip)
  cip$cip <- round(cip$cip, 7)
  
  cip$cip <- ifelse(cip$cip == 0,
                    0,
                    ifelse(cip$cip < 2e-04, 2e-04, cip$cip))
  
  feno_tmp$seged <- NULL
  cip[, c("Event", "seged")] <- NULL
  
  sajat_feno <- prepare_LTFHPlus_input(.tbl = feno_tmp,
                                       CIP = cip, 
                                       age_col = "Time",
                                       aoo_col = "aoo",
                                       CIP_merge_columns = c("Time"),
                                       CIP_cip_col = "cip",
                                       status_col = "Event",
                                       use_fixed_case_thr = F,
                                       fam_id_col = "fid",
                                       personal_id_col = "pid",
                                       interpolation = NA,
                                       min_CIP_value = 1e-5)
  
  sajat_feno <- unique(sajat_feno)
  
  start <- Sys.time()
  sajat_liab <- estimate_liability(.tbl = sajat_feno,
                                   h2 = h2,
                                   pid = "pid",
                                   fam_id = "fid",
                                   role = "role",
                                   tol = 0.1,
                                   out = c(1))
  
  
  
  geno_tmp <- as_FBM(geno_tmp)
  gwas <- big_univLinReg(X = geno_tmp,
                         y.train = sajat_liab$genetic_est,
                         ncores = 1,
                         covar.train = cov_tmp)
  
  stop <- Sys.time()
  Adult_time4[i, "n"] <- db
  Adult_time4[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  save(Adult_time4, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/Adult_time4_origi.RData")
  print(i)
  rm(stop, start, cip, obj.null, geno_tmp, feno_tmp, db, samplek, cov_tmp, sajat_feno, sajat_liab, gwas)
  
  
  
}

Adult_time4[Adult_time4$n > 1000, "ido"] <- Adult_time4[Adult_time4$n > 1000, "ido"] * 60
Adult_time4_mp <- Adult_time4
save(Adult_time4_mp, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/Adult_time4_mp.RData")







############################################################################
### LR teszt

set.seed(1234)
n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 
Score_time <- data.frame(n = numeric(0),
                         ido = numeric(0))

for (i in 1:50){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- as.data.frame(geno[samplek, ])
  
  start <- Sys.time()
  fit0 <- coxph(Surv(Time, Event) ~ X1 + X2,  
                data = feno_tmp)
  
  for(g in colnames(geno_tmp)){
    fit <- coxph(as.formula(paste0("Surv(Time, Event) ~ X1 + X2 + ", g)),  
                 data = cbind(feno_tmp, geno_tmp[g]))
    pek <- anova(fit0, fit, test = "Chisq")$`Pr(>|Chi|)`[2]
  }
  stop <- Sys.time()
  Score_time[i, "n"] <- db
  Score_time[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  rm(stop, start, pek, fit0, geno_tmp, feno_tmp, db, samplek, fit)
  save(Score_time, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/Score_time_origi_mind.RData")
  print(i)
  
}


summary_Score_time <- Score_time %>%
  group_by(n) %>%
  summarise(
    mean_ido = mean(ido),
    sd_ido = sd(ido)
  )

g_Score_time <- ggplot(summary_Score_time, aes(x = n, y = mean_ido)) +
  geom_line() +  # Görbék rajzolása
  geom_point() +  # Pontok rajzolása
  geom_errorbar(aes(ymin = mean_ido - sd_ido, ymax = mean_ido + sd_ido), width = 0.2) +  
  labs(
    x = "Minta szám", 
    y = "Idő" 
  ) +
  theme_minimal() +  # Egyszerű téma
  ggtitle("SPAcox idő")


##########################
### LR teszt 4

set.seed(1234)
n <- rep(c(100, 500, 1000, 5000, 10000), each = 10) 
Score_time4 <- data.frame(n = numeric(0),
                          ido = numeric(0))

num_cores <- 4

for (i in 48:50){
  db <- n[i]
  samplek <- sample(feno$ID, db, replace = F)
  feno_tmp <- feno[feno$ID %in% samplek, ]
  geno_tmp <- as.data.frame(geno[samplek, ])
  
  start <- Sys.time()
  fit0 <- coxph(Surv(Time, Event) ~ X1 + X2,  
                data = feno_tmp)
  
  
  
  p_values <- mclapply(colnames(geno_tmp), function(g) {
    fit <- coxph(as.formula(paste0("Surv(Time, Event) ~ X1 + X2 + ", g)),  
                 data = cbind(feno_tmp, geno_tmp[g]))
    pek <- anova(fit0, fit, test = "Chisq")$`Pr(>|Chi|)`[2]
    return(pek)
  }, mc.cores = num_cores)  
  
  
  stop <- Sys.time()
  Score_time4[i, "n"] <- db
  Score_time4[i, "ido"] <- as.numeric(difftime(stop, start, units = "mins"))
  rm(stop, start, pek, fit0, geno_tmp, feno_tmp, db, samplek, fit)
  save(Score_time4, file = "/srv/__PERSONAL__/Szakdoga/Szakdoga/Teljesitmeny2/LR_time_origi4_mind.RData")
  print(i)
  
}
