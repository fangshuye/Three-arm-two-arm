# setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_III_Three_vs_Two/Simulation/')
# setwd('C:/Users/fye//OneDrive\ -\ Iowa\ State\ University/Research/Project_III_Three_vs_Two/Simulation/')
source(file = "./Project3_fun.R")
library(purrr)
library(logistf)
option_list <- list(
  optparse::make_option("--r", default = 1,
                        help = "which row")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
r=opt$r # r 1:9

# setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_III_Three_vs_Two/')
BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_old <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_pair,sm="OR",comb.fixed = T,comb.random = F)
SE_matrix <- nma_old$seTE.fixed

# get p_nac
p_nac <- BRD_long %>% filter(t == "No active control") %>% summarise(p = sum(r)/sum(n)) %>% pull

# new trt: C
# get the risk table
lor_nac_2_all <- nma_old$TE.fixed[8, ]
risk_all <- lor2prob(p_nac,lor_nac_2_all)
risk_all <- as.data.frame(risk_all)
risk_all$trt <- rownames(risk_all)
rownames(risk_all) <- rep(1:nrow(risk_all))
colnames(risk_all)[1] <- "p"

# 9 rows in total
dat_para <- expand.grid(p_new = seq(0.35,0.43,by = 0.01),
                        s = 100,
                        trt1 = 'Tulathromycin')
trt2 <- 'Trimethoprim'

dat_para <- dat_para[r,]
p_new <- dat_para$p_new[1]
s <- dat_para$s[1]
trt1 <- as.character(dat_para$trt1[1]) 
p_old_1 <- risk_all$p[risk_all$trt == trt1]

sigma2_prev <- SE_matrix[trt1,trt2]^2
p_old_2 <- risk_all$p[risk_all$trt == trt2]

# the optimal allocation in three-arm trial (Z, B, A)
pig_alloc_three <- power_diff(p_new = p_new,
                             p_old_1 = p_old_1,
                             p_old_2 = p_old_2,
                             s = s,
                             sigma2_prev = sigma2_prev,
                             min_B = 10)[2:4]

set.seed(429, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)
nrep <- 50000

table_with_prev <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  # conduct new three-arm study (Z, B, A)
  event_num <- rbinom(3,size = pig_alloc_three,
                      prob = c(p_new, p_old_1, p_old_2))
  elrmdata <- data.frame(count=event_num, x=c('New trt',
                                              trt1, trt2), n=pig_alloc_three)
  binary_dat = pmap_dfr(elrmdata, 
                        function(count, x, n) {
                          data.frame(x = x,
                                     dead = c( rep(1, count),
                                               rep(0, n - count) ) )
                        }
  )
  
  fitmod <- logistf(data=binary_dat, formula=dead~x, pl=FALSE)
  
  est_ZB <- as.numeric(coef(fitmod)[2])  # this is estimate of log(p_B/(1-p_B)) - log(p_Z/(1-p_Z))
  est_ZA <- as.numeric(coef(fitmod)[3]) # this is estimate of log(p_A/(1-p_A)) - log(p_Z/(1-p_Z))
  se_ZB <- sqrt(fitmod$var[2,2])
  se_ZA <- sqrt(fitmod$var[3,3])
  
  binary_dat$x <- factor(binary_dat$x)
  binary_dat <- within(binary_dat, 
                       x <- relevel(x, ref = 'Trimethoprim'))
  fitmod <- logistf(data=binary_dat, formula=dead~x, pl=FALSE)
  est_BA <- as.numeric(coef(fitmod)[3]) # this is estimate of log(p_A/(1-p_A)) - log(p_B/(1-p_B))
  se_BA <- sqrt(fitmod$var[3,3])
  # in BRD_pair, TE is log(p_t1/(1-p_t1)) - log(p_t2/(1-p_t2))
  BRD_new_pair <- BRD_pair
  BRD_new_pair[115,1] <- 100
  BRD_new_pair[115,2:3] <- c(trt1,'New trt')
  BRD_new_pair[115,4:5] <- c(est_ZA, se_ZA)
  
  BRD_new_pair[116,1] <- 100
  BRD_new_pair[116,2:3] <- c(trt2,'New trt')
  BRD_new_pair[116,4:5] <- c(est_ZB, se_ZB)
  
  BRD_new_pair[117,1] <- 100
  BRD_new_pair[117,2:3] <- c(trt1,trt2)
  BRD_new_pair[117,4:5] <- c(est_BA, se_BA)
  
  
  nma_updated <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_pair,sm="OR",comb.fixed = T,comb.random = F)

  p_value_three <- nma_updated$pval.fixed["New trt",trt1]
  res_three_arm <- ifelse(p_value_three > 0.05,0,1)
  
  res_three_arm
  
}

BE <- table_with_prev


dat_para[1,4] <- trt2
dat_para[1,5:7] <- pig_alloc_three
dat_para[1,8] <- power_three_arm_given_allo(p_new,p_old_1,p_old_2,sigma2_prev,pig_alloc_three)
dat_para[1,9] <- apply(BE, 2, mean)*100
colnames(dat_para)[4:9] <- c('trt B','n_Z','n_A','n_B',
                             'formula_power_three','sim_power_three')

write.csv(dat_para, 
          file = paste0("./Res/power_three_arm",'_V19_',r,".csv"), 
          row.names = F)
