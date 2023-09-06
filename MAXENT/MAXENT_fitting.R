if (!require(dplyr)) install.packages("dplyr")

if (!require(dismo)) install.packages("dismo")

if (!require(biomod2)) install.packages("biomod2")

if (!require(sf)) install.packages("sf")

Sys.setenv(JAVA_HOME='C:/Program Files/Java/jre1.8')


Knames <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/Knames.rds")
#saveRDS(df_species, "E:/Madariaga/Documents/Phd/Manuscript2/df_species.rds")
df_species <-readRDS("E:/Madariaga/Documents/Phd/Manuscript2/df_species.rds")

formulas <- readRDS("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/Server/Kinzig_upd/Results/dredge_formulasK_upd.rds")
###Create function to unregister parallel 
cv_ids <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/cv_ids.rds")
#unregister_dopar <- function() {
#  env <- foreach:::.foreachGlobals
#  rm(list=ls(name=env), pos=env)
#}

###Create parallel 

#cl <- makePSOCKcluster(detectCores()-3, outfile="")
#parallel::mcaffinity(1:5)
#doParallel::registerDoParallel(cl)
#getDoParWorkers()
maxent_param <- function(data, occ, pred, k = 5){
  
  
  if (!require(dismo)) install.packages("dismo")
  if (!require(caret)) install.packages("caret")
  if (!require(precrec)) install.packages("precrec")
  if (!require(rJava)) install.packages("rJava")
  
  # generate balanced CV folds
  folds <- caret::createFolds(y = as.factor(data[,occ]), k = k)
  
  # regularisation multipliers
  ms <- c(0.5, 1, 2, 3, 4)
  grid <- expand.grid(
    regmult = paste0("betamultiplier=", ms),
    features = list(
      c("noautofeature", "nothreshold"), # LQHP
      c("noautofeature", "nothreshold", "noproduct"), # LQH
      c("noautofeature", "nothreshold", "nohinge", "noproduct"), # LQ
      c("noautofeature", "nothreshold", "nolinear", "noquadratic", "noproduct"), # H
      c("noautofeature", "nothreshold", "noquadratic", "nohinge", "noproduct")), # L
    stringsAsFactors = FALSE
  )
  AUCs <- c()
  for(n in seq_along(grid[,1])){
    full_pred <- data.frame()
    for(i in seq_len(length(folds))){
      trainSet <- unlist(folds[-i])
      testSet <- unlist(folds[i])
      if(inherits(try(
        maxmod <- dismo::maxent(x = data[trainSet, pred],
                                p = data[trainSet,occ],
                                removeDuplicates = FALSE,
                                args = as.character(unlist(grid[n, ]))
        )
      ), "try-error")){
        next
      }
      modpred <- predict(maxmod, data[testSet, pred], args = "outputformat=cloglog")
      pred_df <- data.frame(score = modpred, label = data[testSet,occ])
      full_pred <- rbind(full_pred, pred_df)
    }
    AUCs[n] <- precrec::auc(precrec::evalmod(scores = full_pred$score, 
                                             labels = full_pred$label))[1,4]
  }
  best_param <- as.character(unlist(grid[which.max(AUCs), ]))
  return(best_param)
}

source("E:/Madariaga/Documents/Phd/SSN/functions_modified.R")



set.seed(432)

maxent_list_cv <- list()

for(k in Knames) {
  for (i in 1:5){
  
  name <- k
  print(name)
  Fold <- paste0("Fold.", i)
  maxent_obj <- df_species[[name]]
  
  # Get the ssn object from the ssn_list
  id_list <- cv_ids[[name]][["splits"]][[i]][["in_id"]]
  
  maxent_obj[[paste0(name, Fold)]]<- NA
  
  subset_rows <- maxent_obj$id %in% id_list
  maxent_obj[[paste0(name,"Fold.", i)]][subset_rows]<- maxent_obj[[name]][subset_rows]
  print("check ssn_obj_df")
  

  Calib <- subset(maxent_obj,!is.na(maxent_obj[paste0(name,Fold)]))
  pred <- c("dh4", "fl1","ra4","th3","bio2","bio8","bio9","bio19") #predictor names
  
  print(pred)
  
    params <- maxent_param (data = Calib,
                            pred = pred,
                            occ = name,
                            k = 5)
    print(params)
   
    
    Maxent_model <- maxent(x = Calib[, pred], p = Calib[, paste0(name,Fold)],
                          args = params)
    
    
    #print(Maxent_model)
    
    Eval <- subset(maxent_obj,is.na(maxent_obj[paste0(name,Fold)]))
    PredProbEval <- predict(Maxent_model, Eval, args = "outputformat=cloglog")
    idsPresabsProb <- data.frame(Eval[,c("id",name)],PredProbEval)
    threshold <- eval_mod(idsPresabsProb)$threshold
    performance <- eval_mod(idsPresabsProb)$performance 
    
    print(performance)
    
    maxent_list_cv[[name]][[Fold]] <- list("model"= Maxent_model,
                                 "threshold"= threshold,
                                 "performance"=performance)
  }
}  

saveRDS(maxent_list_cv, "E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/maxent_models_cv.rds" )



#maxent_list <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/maxent_models.rds")
####


eval_maxent_cv_dflist <-list()
for (k in Knames){
  species_results <- data.frame()
  for (i in 1:5){
    Fold <- paste0("Fold.",i)
  if (is.list(maxent_list_cv[[k]][[Fold]])){
    iteration_df<- data.frame(
      Species= k,
      Fold= Fold,
      model= "maxent_cv",
      metric= maxent_list_cv[[k]][[Fold]][["performance"]][["Measures"]],
      values =maxent_list_cv[[k]][[Fold]][["performance"]][["Value"]]
    )
    species_results <- rbind(species_results, iteration_df)
  }  
  }
  eval_maxent_cv_dflist[[k]] <- species_results
}


maxent_cv_df<-do.call(rbind, eval_maxent_cv_dflist)

saveRDS(maxent_cv_df, "E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/maxent_cv_df.rds" )

mean_maxent_cv <- maxent_cv_df%>%
  group_by(Species, model, metric) %>%
  summarize(mean_value = mean(values))


saveRDS(mean_maxent_cv, "E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/mean_maxent_cv_df.rds" )


set.seed(432)

maxent_list_independent <- list()

for(k in 1:length(Knames)) {
  
    name <- Knames[[k]]
    print(name)
    maxent_obj <- df_species[[name]]
    maxent_obj[,paste0(name,".eva")]<- maxent_obj[name]
    
    maxent_obj[,paste0(name,".eva")]<- ifelse(maxent_obj$dataset == "Senckenberg", NA, maxent_obj[,paste0(name,".eva")])    
    
    print("check ssn_obj_df")
    
    
    Calib <- subset(maxent_obj,!is.na(maxent_obj[paste0(name,".eva")]))
    pred <- c("dh4", "fl1","ra4","th3","bio2","bio8","bio9","bio19") #predictor names
    
    print(pred)
    
    params <- maxent_param (data = Calib,
                            pred = pred,
                            occ = paste0(name,".eva"),
                            k = 5)
    print(params)
    
    
    Maxent_model <- maxent(x = Calib[, pred], p = Calib[, paste0(name,".eva")],
                           args = params)
    
    
    #print(Maxent_model)
    
    Eval <- subset(maxent_obj,is.na(maxent_obj[paste0(name,".eva")]))
    PredProbEval <- predict(Maxent_model, Eval, args = "outputformat=cloglog")
    idsPresabsProb <- data.frame(Eval[,c("id",name)],PredProbEval)
    threshold <- eval_mod(idsPresabsProb)$threshold
    performance <- eval_mod(idsPresabsProb)$performance 
    
    print(performance)
    
    maxent_list_independent[[name]] <- list("model"= Maxent_model,
                                           "threshold"= threshold,
                                           "performance"=performance)
  }
  

saveRDS(maxent_list_independent, "E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/maxent_list_independent_models.rds" )


maxent_independent_dflist <-list()
for (k in Knames){
  species_results <- data.frame()
  if (is.list(maxent_list_independent[[k]])){
    iteration_df<- data.frame(
      Species= k,
      model= "Maxent_independent",
      metric= maxent_list_independent[[k]][["performance"]][["Measures"]],
      values = maxent_list_independent[[k]][["performance"]][["Value"]]
    )
    species_results <- rbind(species_results, iteration_df)
  }  
  maxent_independent_dflist[[k]] <- species_results
}  

maxent_independent_result_df<-do.call(rbind, maxent_independent_dflist)
saveRDS(maxent_independent_result_df, "E:/Madariaga/Documents/Phd/Manuscript2/MAXENT/maxent_independent_df.rds" )
