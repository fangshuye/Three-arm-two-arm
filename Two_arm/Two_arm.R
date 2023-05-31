packages <- c("dplyr", "parallel", "foreach", "doParallel", 
              "doRNG","DEoptimR","purrr","logistf")
lapply(packages, library, character.only = TRUE)
con_fn_two_arm <- function(n, p_new, p_old,s){
  n1 <- n[1]
  n2 <- n[2]
  n1 + n2 - s
}

var_fn_two_arm <- function(n, p_new, p_old, s){
  mu_1 <- p_new*(1-p_new)
  mu_2 <- p_old*(1-p_old)
  1/(mu_1 * n[1]) + 1/(mu_2 * n[2])
}

power_two_arm_given_allo <- function(p_new,p_old_1,n){
  
  var_BZ <- var_fn_two_arm(n, p_new, p_old_1, s = 0)
  se <- sqrt(var_BZ)
  size <- log(p_new/(1-p_new)) - log(p_old_1/(1-p_old_1))
  z <- size/se
  power <- pnorm(z-qnorm(0.975))+pnorm(-z-qnorm(0.975))
  return(power*100)
  
}

###
option_list <- list(
  optparse::make_option("--r", default = 1,
                        help = "which row")
)
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
r=opt$r

# setwd('/Users/fangshu/OneDrive\ -\ Iowa\ State\ University/Research/Project_III_Three_vs_Two/')
# setwd('C:\\Users\\fye\\OneDrive - Iowa State University\\Research\\Project_III_Three_vs_Two')
risk_all <-  read.csv("./data/risk_nma.csv")

dat_para <- expand.grid(p_new = seq(0.35,0.43,by = 0.01),
                        s = 100,
                        trt1 = 'Tulathromycin')

dat_para <- dat_para[r,]
p_new <- dat_para$p_new[1]
s <- dat_para$s[1]
trt1 <- as.character(dat_para$trt1[1]) 
p_old_1 <- risk_all$p[risk_all$trt == trt1]


# the optimal allocation in two-arm trial (Z, B)
pig_alloc_two <- JDEoptim(lower = c(1,1), upper = c(s,s), 
                          fn = var_fn_two_arm, 
                          constr = con_fn_two_arm, 
                          meq = 1, triter = 50, 
                          p_new = p_new,
                          p_old = p_old_1,
                          s = s)
# (n_new, n_old_1)
pig_alloc_two <- round(pig_alloc_two$par,0)



set.seed(111)
#registerDoParallel(cores = 16)
nrep <- 100000

table_with_prev <- c()
for(i in 1:nrep){
  event_num <- rbinom(2,size = pig_alloc_two,prob = c(p_new,p_old_1))
  elrmdata <- data.frame(count=event_num, x=c(0,1), n=pig_alloc_two)
  binary_dat = pmap_dfr(elrmdata, 
                        function(count, x, n) {
                          data.frame(x = x,
                                     dead = c( rep(1, count),
                                               rep(0, n - count) ) )
                        }
  )
  fitmod <- logistf(data=binary_dat, formula=dead~x, pl=FALSE)
  est_ZA <- coef(fitmod)[2] # this is estimate of log(p_B/(1-p_B)) - log(p_Z/(1-p_Z))
  se_ZA <- sqrt(fitmod$var[2,2])
  p_value <- (1-pnorm(abs(est_ZA)/se_ZA))*2
  table_with_prev <- rbind(table_with_prev,ifelse(p_value>0.05,0,1))
}

dat_para[1,4] <- power_two_arm_given_allo(p_new,p_old_1,pig_alloc_two)
dat_para[1,5] <- apply(table_with_prev,2,mean)*100

dat_para[1,6:7] <- pig_alloc_two

colnames(dat_para)[4:7] <- c('formula_power_two','sim_power_two',
                             'n_Z','n_A')
write.csv(dat_para, 
          file = paste0("./Res/power_two_arm_direct",'_V19_',r,".csv"), 
          row.names = F)
