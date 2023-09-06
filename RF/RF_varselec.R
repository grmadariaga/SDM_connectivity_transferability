library(tidymodels)
library(randomForest)
library(dplyr)
library(SSN)
library(SSN)
library(sirad)
library(parallel)
library(doParallel)
library(foreach)

Knames <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/Knames.rds")
#saveRDS(df_species, "E:/Madariaga/Documents/Phd/Manuscript2/df_species.rds")
df_species <-readRDS("E:/Madariaga/Documents/Phd/Manuscript2/df_species.rds")
formula_list <- readRDS("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/Server/Kinzig_upd/Results/dredge_formulasK_upd.rds")


cv_ids <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/cv_ids.rds")

source("E:/Madariaga/Documents/Phd/SSN/functions_modified.R")






####Training models in one dataset and assessing them in other

rf_list <- list()

for(k in 1:length(Knames)){
  
name <- Knames[[k]]
print(name)

# Get the ssn object from the ssn_list

rf_df <- df_species[[name]]


rf_df_sub_cal<- rf_df[rf_df$dataset=="Karan",]
print("check rf_df_sub_cal")

  prNum <- as.numeric(table(rf_df_sub_cal[[name]])["1"])          
  bgNum <- as.numeric(table(rf_df_sub_cal[[name]])["0"])
  samsize <- c("0" = bgNum, "1" = prNum)
  print(samsize)
  
  formula <- formula_list[[name]][["formula"]]
  frml <- as.formula(paste(paste0("as.factor(",name,")"), "~", paste0(" ", formula)[3]))
  
  print(frml)
    tryCatch(
      {
    
    rf_model <- randomForest(formula= frml, 
                             data = rf_df_sub_cal,
                             ntree = 1000,
                             sampsize = samsize,
                             replace = TRUE,
                             importance = TRUE,
                             na.action= na.omit)
    print(rf_model)
    
   
    Eval <- rf_df[rf_df$dataset=="Senckenberg",]
    PredProbEval <- as.numeric(predict(rf_model, Eval, type = "prob")[,"1"])
    idsPresabsProb <- data.frame(Eval[,c("id",name)],PredProbEval)
    threshold <- eval_mod(idsPresabsProb)$threshold
    performance <- eval_mod(idsPresabsProb)$performance 
    
    print(performance)
    
    rf_list[[name]] <- list("model"= rf_model,
                                 "threshold"= threshold,
                                 "performance"=performance)},
    error =function(e) {
      rf_list[[name]] <- list("model"= "failed",
                                    "threshold"= "failed",
                                    "performance"="failed")
    }
    )
}

  

####


rf_dflist <-list()
for (k in Knames){
  species_results <- data.frame()
    if (is.list(rf_list[[k]])){
      iteration_df<- data.frame(
        Species= k,
        model= "RF_independent",
        metric= rf_list[[k]][["performance"]][["Measures"]],
        values = rf_list[[k]][["performance"]][["Value"]]
      )
      species_results <- rbind(species_results, iteration_df)
    }  
    rf_dflist[[k]] <- species_results
}  

rf_result_list<-do.call(rbind, rf_dflist)


saveRDS(rf_result_list,"E:/Madariaga/Documents/Phd/Manuscript2/RF/rf_indpndnt.rds")


####Training models and assessing them with cv
rf_list2 <- list()
for(k in 1:length(Knames)){
  
  for (i in 1:5){
    tryCatch(
      {
  
  name <- Knames[[k]]
  print(name)
  Fold <- paste0("Fold.", i)
  rf_obj <- df_species[[name]]
  
  # Get the ssn object from the ssn_list
  id_list <- cv_ids[[name]][["splits"]][[i]][["in_id"]]
  
  rf_obj[[paste0(name, Fold)]]<- NA
  
  subset_rows <- rf_obj$id %in% id_list
  rf_obj[[paste0(name,"Fold.", i)]][subset_rows]<- rf_obj[[name]][subset_rows]
  calib<- subset(rf_obj, !is.na(rf_obj[[paste0(name,Fold)]]))
  print("check ssn_obj_df")
  prNum <- as.numeric(table(calib[[paste0(name,Fold)]])["1"])          
  bgNum <- as.numeric(table(calib[[paste0(name,Fold)]])["0"])
  samsize <- c("0" = bgNum, "1" = prNum)
  print(samsize)
  formula <- formula_list[[k]]$formula
  form_text<- as.formula(paste(paste0("as.factor(",name,Fold, ")"), "~", paste0(" ", formula)[3]))
  print(form_text)
     
        
        rf_model <- randomForest(formula= as.formula(form_text), 
                                 data = calib,
                                 ntree = 1000,
                                 sampsize = samsize,
                                 replace = TRUE,
                                 importance = TRUE,
                                 na.action= na.omit)
        print(rf_model)
        
        Eval <- subset(rf_obj,is.na(rf_obj[[paste0(name,Fold)]]))
        PredProbEval <- as.numeric(predict(rf_model, Eval, type = "prob")[,"1"])
        idsPresabsProb <- data.frame(Eval[,c("id",name)],PredProbEval)
        threshold <- eval_mod(idsPresabsProb)$threshold
        performance <- eval_mod(idsPresabsProb)$performance 
        
        print(performance)
        
        rf_list2[[name]][[Fold]] <- list("model"= rf_model,
                                      "threshold"= threshold,
                                      "performance"=performance)},
      error =function(e) {
        rf_list2[[name]][[Fold]] <- list("model"= "failed",
                                      "threshold"= "failed",
                                      "performance"="failed")
      }
    )
  }
}

####


rf_dflist2 <-list()
for (k in Knames){
  species_results <- data.frame()
  for (i in 1:5){
    Fold <- paste0("Fold.",i)
    if (is.list(rf_list2[[k]][[Fold]])){
      iteration_df<- data.frame(
        Species= k,
        Fold= Fold,
        model= "RF_cv",
        metric= rf_list2[[k]][[Fold]][["performance"]][["Measures"]],
        values = rf_list2[[k]][[Fold]][["performance"]][["Value"]]
      )
      species_results <- rbind(species_results, iteration_df)
    }  
  }
  rf_dflist2[[k]] <- species_results
}  

rf_df2<-do.call(rbind, rf_dflist2)


mean_rf2 <- rf_dflist2 %>%
  bind_rows() %>%
  group_by(Species, model, metric) %>%
  summarize(mean_value = mean(values))

saveRDS(mean_rf2,"E:/Madariaga/Documents/Phd/Manuscript2/RF/rf_cv.rds")
