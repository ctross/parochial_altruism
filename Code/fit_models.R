################################################################# Fit Models                                                             
Iter <- 1000
Warmup <- 500   
Thin <- 1
Cores <- 1
Chains <- 1
Refresh <- 10
Seed <- 123
MTD <- 13
AD <- 0.96

model_dat <- vector("list",4)

model_dat_Coast$Q <- c(1,1)
model_dat_Inland$Q <- c(1,1)
model_dat[[1]] <- model_dat_Coast 
model_dat[[2]] <- model_dat_Inland

model_dat_Coast$Q <- c(0,1)
model_dat_Inland$Q <- c(0,1)
model_dat[[3]] <- model_dat_Coast
model_dat[[4]] <- model_dat_Inland

fit <- mclapply(1:4,function(z){
         stan(file="Code/SRM.stan",data=model_dat[[z]],thin=Thin,iter=Iter,warmup=Warmup,cores=Cores,chains=Chains,refresh=Refresh,
         	  control=list(max_treedepth=MTD,adapt_delta=AD))
             },
 mc.cores = 1)    
 
 save.image("results.RData")
 