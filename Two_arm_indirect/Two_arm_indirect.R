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
r=opt$r # r 1:40

# setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_III_Three_vs_Two/')
BRD <- read.csv("./data/updated_dataset.csv", stringsAsFactors = F)
BRD_r <- BRD 
BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] <- BRD_r$Number.of.Event.in.arm.2[BRD_r$Study.number == 2] - 0.5
BRD_long <- wide2long(BRD_r)
BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, studlab = id, data = BRD_long, allstudies = T, sm = "OR")
BRD_pair <- BRD_pair[,c('studlab','treat1','treat2','TE','seTE')]
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


# the optimal allocation in two-arm trial (Z, A)
pig_alloc_two_AZ <- JDEoptim(lower = c(1,1), upper = c(s,s), 
                             fn = var_fn_two_arm, 
                             constr = con_fn_two_arm, 
                             meq = 1, triter = 50, 
                             p_new = p_new,
                             p_old = p_old_2,
                             s = s)
# (n_new, n_old_2)
pig_alloc_two_AZ <- round(pig_alloc_two_AZ$par,0)



set.seed(400, kind = "L'Ecuyer-CMRG")
registerDoParallel(cores = 16)
nrep <- 50000

table_with_prev <- foreach (rep=1:nrep, .combine = rbind)%dopar%{
  
  event_num <- rbinom(2,size = pig_alloc_two_AZ,
                      prob = c(p_new,p_old_2))
  elrmdata <- data.frame(count=event_num, x=c(0,1), n=pig_alloc_two_AZ)
  binary_dat = pmap_dfr(elrmdata, 
                        function(count, x, n) {
                          data.frame(x = x,
                                     dead = c( rep(1, count),
                                               rep(0, n - count) ) )
                        }
  )
  fitmod <- logistf(data=binary_dat, formula=dead~x, pl=FALSE)
  est_ZB <- coef(fitmod)[2] # this is estimate of log(p_B/(1-p_B)) - log(p_Z/(1-p_Z))
  se_ZB <- sqrt(fitmod$var[2,2])
  # in BRD_pair, TE is log(p_t1/(1-p_t1)) - log(p_t2/(1-p_t2))
  BRD_new_pair <- BRD_pair
  BRD_new_pair[115,1] <- 100
  BRD_new_pair[115,2:3] <- c(trt2,'New trt')
  BRD_new_pair[115,4:5] <- c(est_ZB, se_ZB)
  
  nma_updated <- netmeta(TE,seTE,treat1,treat2,studlab,data=BRD_new_pair,sm="OR",comb.fixed = T,comb.random = F)
  p_value_BZ <- nma_updated$pval.fixed["New trt",trt1]
  res_two_arm_ind <- ifelse(p_value_BZ > 0.05,0,1)

  res_two_arm_ind
  
}

BE <- table_with_prev

dat_para[1,4] <- trt2
dat_para[1,5:6] <- pig_alloc_two_AZ
dat_para[1,7] <- power_two_arm_AZ(p_new = p_new,
                                  p_old_1 = p_old_1,
                                  p_old_2 = p_old_2,
                                  pig_alloc_c = pig_alloc_two_AZ,
                                  sigma2_prev = sigma2_prev)
dat_para[1,8] <- apply(BE, 2, mean)*100
colnames(dat_para)[4:8] <- c('trt A','n_Z','n_B',
                             'formula_power_two_BZ','sim_power_two_BZ')

write.csv(dat_para, 
          file = paste0("./Res/power_two_arm_indirect",'_V19_',r,".csv"), 
          row.names = F)
