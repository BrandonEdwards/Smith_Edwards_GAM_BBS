###################################################
# Adam C. Smith & Brandon P.M. Edwards
# Script to re-run the models and cross-validation for 4 of the species 
# The original model runs for these species included at least one model for which
# > 1% of the n-values failed to converge 
# these species include the model GAM for Carolina Wren, Chestnut-collared Longspur, and Pine Siskin
# as well s the DIFFERENCE model for Chestnut-collare Longspur and Cooper's Hawk
# Because the BBS data has had additional 2018 observations added since the original models were run
# All models for each of these species need to be re-run to keep the data identical across all models for a given species
# 
#
# This script checks for convergence using the posterior package, if any of the 
# n-values have failed to converge, it re-starts the MCMC sampler using the final 
# values from the earlier run as initial values, and increase the thinning rate and number of iterations
#
# This script also uses a different version of the DIFFERENCE model that more closely follows
# the centered version outlined in Link and Sauer 2020.
# This alternate version of the DIFFERENCE model centers the time-series in the mid-year
# and estimates each year-effect as a function of the previous year in the second half of the time-series
# and a function of the next-year in the early half of the time-series
###################################################
remove(list = ls())
K <- 15 #k = 15 fold X-valid
n_iter = 20000
n_thin = 20
n_burnin = 20000
n_chains = 3
n_adapt = NULL



dir.create("output", showWarnings = F)

# Install v1.1.2 from Github 

library(bbsBayes)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)
library(posterior)
library(tidyverse)

# Only need to run this if you don't have BBS data saved
# in directory on your computer
# yes immediately following to agree to terms and conditions
####################################
# fetch_bbs_data()
# yes

stratified_data <- stratify(by = "bbs_usgs")


species_to_run <- c("Chestnut-collared Longspur",
                    "Cooper's Hawk",
                    "Carolina Wren",
                    "Pine Siskin") # a selection of species for which > 1% of n values failed to converge for at least one of the models


models <- c("firstdiff", "gam", "gamye","slope")


###################################################
# Analysis by Species X Model Combination
###################################################

for (species in species_to_run)
{
  sp.dir = paste0("output/", species)
  dir.create(sp.dir)
  
  
  #### identifying the K folds for cross-validation
    ## selecting stratified samples that remove 10% of data within each stratum
  jags_data <- prepare_jags_data(strat_data = stratified_data,
                                 species_to_run = species,
                                 min_max_route_years = 3,
                                 model = "slope",
                                 max_year = 2018)
  
  sp.k = paste0(sp.dir, "/sp_k.RData")
  
  if(file.exists(sp.k) == FALSE){
  
    kk = vector(mode = "integer",length = jags_data$ncounts)
    for(i in 1:jags_data$nstrata){
      set.seed(2019)
      wstrat = which(jags_data$strat == i)
      
      kk[wstrat] = as.integer(ceiling(runif(length(wstrat),0,K))) 
    }
    
    ### saving the k-fold identifiers so that they're the same across all models
    save(kk,file = sp.k)
  }
    
    # Set up parallel stuff
  # requires only 4 cores to run each model in parallel for species
    n_cores <- length(models)
    cluster <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cluster)
    
    message(paste0("Beginning full model analysis for species ", species))
    
foreach(m = 1:4,
        .packages = c('bbsBayes','tidyverse','posterior'),
        .inorder = FALSE,
        .errorhandling = "pass") %dopar%
    {
      
      model = models[m]
      model_dir <- paste0(sp.dir,
                        "/",
                        model)
    dir.create(model_dir)
    
    jags_data <- prepare_jags_data(strat_data = stratified_data,
                                   species_to_run = species,
                                   min_max_route_years = 3,
                                   model = model,
                                   max_year = 2018)
    
    
    if(model == "firstdiff"){
      jags_data[["fixedyear"]] <- jags_data$ymin+round((jags_data$ymax-jags_data$ymin)/2)
    
      }
    load(sp.k)
    jags_data$ki <- kk
    
    save(jags_data,file = paste0(model_dir,"/jags_data.RData"))
    ##################### FULL MODEL RUN ##########################
  

        #inits = function()
    model_file <- paste0("loo-models/", model, "-t.jags") #using hte heavy-tailed error distribution
    if(model == "firstdiff" & species == "Cooper's Hawk"){
      model_file <- paste0("loo-models/", model, "alt-t.jags") #using the alternate first difference structure centered on the midyear    
      }
    jags_mod_full <- run_model(jags_data = jags_data,
                               model_file_path = model_file,
                               n_iter = n_iter,
                               #n_adapt = n_adapt,
                               n_burnin = n_burnin,
                               n_chains = n_chains,
                               n_thin = n_thin,
                               parallel = FALSE,
                               modules = NULL)
 
     n_iter2 <- n_iter+10000
     n_thin2 <- n_thin+10
     attempts <- 0
     
     tmp = posterior::as_draws(jags_mod_full$samples)
     sumt = posterior::summarise_draws(tmp) %>% 
       as.data.frame
     
     while(any(sumt$rhat > 1.1) & attempts < 3){
       attempts <- attempts+1
       n_iter2 <- n_iter2+10000
       n_thin2 <- n_thin2+10
       
       new.inits <- get_final_values(jags_mod_full)
       jags_mod_full <- run_model(jags_data = jags_data,
                                  model_file_path = model_file,
                                  n_iter = n_iter2,
                                  inits = new.inits,
                                  n_burnin = 0,
                                  n_chains = n_chains,
                                  n_thin = n_thin2,
                                  parallel = FALSE,
                                  modules = NULL)
       tmp = posterior::as_draws(jags_mod_full$samples)
       sumt = posterior::summarise_draws(tmp) %>% 
         as.data.frame
       save(jags_mod_full, file = paste0(model_dir, "/jags_mod_full.RData"))
       save(list = c("sumt","attempts"), file = paste0(model_dir, "/jags_mod_converge.RData"))
     }
     save(jags_mod_full, file = paste0(model_dir, "/jags_mod_full.RData"))
     save(c("sumt","attempts"), file = paste0(model_dir, "/jags_mod_converge.RData"))
     
    ##################### TRENDS AND TRAJECTORIES #################
     ### this trend calculation and plotting is not necessary, only for model-checking 
    dir.create(paste0(model_dir, "/plots"))
    
    # Stratum level
     # Stratum level
     inds = generate_indices(jags_mod = jags_mod_full,jags_data = jags_data)
     
     
     
     s_plots <- plot_indices(inds)
     
     for (i in 1:length(s_plots))
     {
       png(filename = paste0(model_dir, "/plots/", names(s_plots[i]), ".png"))
       print(s_plots[[i]])
       dev.off() 
     }
     
    
    }#end of full model parallel loop
    
stopCluster(cl = cluster)














for(model in models){
    
    ##################### CROSS VALIDATION ########################
  message(paste0("Beginning cross validation for model ", model, ", species ", species))
  model_dir <- paste0(sp.dir,
                      "/",
                      model)
  dir.create(paste0(model_dir, "/cv"))
    
  load(paste0(model_dir,"/jags_data.RData"))
  load(paste0(model_dir, "/jags_mod_full.RData"))
  
    inits <- get_final_values(model = jags_mod_full)
    model_file <- paste0("loo-models/", model, "-loo.jags")
    if(model == "firstdiff" & species == "Cooper's Hawk"){
      model_file <- paste0("loo-models/", model, "alt-loo.jags") #using the alternate first difference structure centered on the midyear    
    }
    # Set up parallel stuff
    # requires n_cores cores to run each cross-validation fold in parallel
    n_cores <- K
    cluster <- makeCluster(n_cores, type = "PSOCK")
    registerDoParallel(cluster)
    
    posterior <- foreach(kk = 1:K,
                         .packages = 'bbsBayes',
                         .inorder = FALSE,
                         .errorhandling = "pass",
                         .combine = 'cbind') %dopar%
      {
        
        indices_to_remove <- which(jags_data$ki == kk)
        true_count_k <- jags_data$count[indices_to_remove]
        
        n_remove <- as.integer(length(indices_to_remove))
        
        jags_data_loo <- jags_data
        jags_data_loo$count[indices_to_remove] <- NA
        jags_data_loo$I <- indices_to_remove
        jags_data_loo$Y <- true_count_k
        jags_data_loo$nRemove <- n_remove
        
        # Run model and track some LOOCV variables
        # Save the original models to disk using model_to_file() then
        # modify the models to add in the LOOCV variable tracking. 
        # Then give run_model() the path to that new model.
        jags_mod_loo <- run_model(jags_data = jags_data_loo,
                                  model_file_path = model_file,
                                  parameters_to_save = c("LambdaSubset"),
                                  track_n = FALSE,
                                  inits = inits,
                                  n_iter = n_iter,
                                  #n_adapt = n_adapt,
                                  n_burnin = 3000,
                                  n_chains = n_chains,
                                  n_thin = n_thin,
                                  parallel = FALSE,
                                  modules = NULL)
        
        # Just comment this line out if you don't want to save individual
        # model runs for each year left out
        save(jags_mod_loo, 
             file = paste0(model_dir, "/cv/k_", kk, " removed.RData"))
        
        jags_mod_loo$sims.list$LambdaSubset
      }
    
    stopCluster(cl = cluster)
    

  }

}#end species loop
