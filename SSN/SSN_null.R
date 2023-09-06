library(tidymodels)
library(SSN)
library(sirad)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(caret)


###Load SSN object Kinzig


#null_AIC_list <-readRDS("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/server/Kinzig_upd/Results/dredge_formulasK_upd.rds")

options(na.action = "na.omit")


cv_ids <- ("E:/Madariaga/Documents/Phd/Manuscript2/cv_ids.rds")
#ssn_obj0 <- glmssn(Asellus.aquaticus.ev ~ ra6, ssn.object = ssn_list$Asellus.aquaticus,
#                  family = "binomial",
#                        CorModels = NULL)



glm_null_models <- list()

# Loop through each name in Knames
for (name in Knames) {
  for (i in 1:5) {
  
  tryCatch(
    {
      # Construct the column name with .ev suffix
      Fold <- paste0("Fold.", i)
      ssn_obj_i <- ssn_list[[name]][["obj"]]
      
      # Get the ssn object from the ssn_list
      id_list <- ssn_list[[name]][["folds"]][[1]][[i]][["in_id"]]
      
      ssn_obj_df <- getSSNdata.frame(ssn_obj_i, "Obs")
      ssn_obj_df[[paste0(name, Fold)]]<- NA
      
      subset_rows <- ssn_obj_df$id %in% id_list
      ssn_obj_df[[paste0(name,"Fold.", i)]][subset_rows]<- ssn_obj_df[[name]][subset_rows]
      print("check ssn_obj_df")
            
      ssn_obj_subset <- putSSNdata.frame(ssn_obj_df, K)
      
      frml <- null_AIC_list[[name]][["formula"]]%>%as.character()
      frml <- as.formula(paste(paste0(frml[2],Fold), frml[1], frml[3], collapse = " "))
      
      
      # Run the glmssn model with the .ev column
      model <- glmssn(
        formula = frml,
        ssn.object = ssn_obj_subset,
        family = "binomial",
        CorModels = NULL
      )
      
      # Save the model in the glm_null_models list with the name
      glm_null_models[[name]][[Fold]] <- model},  
    error =function(e) {
      # Store the error message in the list with the name from Knames
      glm_null_models[[name]][[Fold]] <<- paste("Failed model:", name, conditionMessage(e))
    }
  )
  }
}

source("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/functions_modified.R")

eval_results_null_glm <- list()

# Loop through each name in Knames
for (name in Knames) {
  for (i in 1:5) {
    

  tryCatch(
    {
      Fold <- paste0("Fold.", i)
      # Get the corresponding glmssn model from glm_null_models
      model <- glm_null_models[[name]][[Fold]]
      
      # Create the "_MissingObs_" evaluation dataframe
      Evaluation <- getSSNdata.frame(model, "_MissingObs_")
      
      # Create Eval_obj using the evaluation dataframe and ssn_obj0
      Eval_obj <- putSSNdata.frame(Evaluation, model, "_MissingObs_")
      
      # Predict using ssn_obj0 and the "_MissingObs_" evaluation dataframe
      Predeval <- predict(model, "_MissingObs_")
      
      # Calculate SSNProb for Prediceval
      SSN0_predev <- SSNProb(Predeval)
      
      # Subset pointdata and select relevant columns for evaluation
     specific_df <- glm_null_models[[name]][[Fold]][["ssn.object"]]@obspoints@SSNPoints[[1]]@point.data
     SSN0presabs <- specific_df[is.na(specific_df[[paste0(name,Fold)]]), c("id", name)] 
     
     # Combine SSN0presabs with SSN0_predev
      idssn0presabs <- data.frame(SSN0presabs, SSN0_predev)
      
      # Perform the evaluation using eval_mod on idssn0presabs
      result <- eval_mod(idssn0presabs)
      
      # Store the evaluation result in the list with the name from Knames
      eval_results_null_glm[[name]][[Fold]] <- result
    },  
    error =function(e) {
      # Store the error message in the list with the name from Knames
      eval_results_null_glm[[name]][[Fold]] <<- paste("Failed model:", name, conditionMessage(e))
    }
  )
  }
}


saveRDS(glm_null_models, "E:/Madariaga/Documents/Phd/Manuscript2/SSN/NULL/null_models.rds")
saveRDS(eval_results_null_glm, "E:/Madariaga/Documents/Phd/Manuscript2/SSN/NULL/results_null_models.rds")

glm_null_models<- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SSN/NULL/null_models.rds")
eval_results_null_glm <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SSN/NULL/results_null_models.rds")



null_df_list <- list()

for (k in Knames){
  species_results <- data.frame()
  for (i in 1:5){
    Fold <- paste0("Fold.", i)
    if (is.list(eval_results_null_glm[[k]][[Fold]])){
      iteration_df<- data.frame(
        Species= k,
        iteration= i,
        model= "glm_cv",
        metric= eval_results_null_glm[[k]][[Fold]][["performance"]][["Measures"]],
        values = eval_results_null_glm[[k]][[Fold]][["performance"]][["Value"]]
      )
      species_results <- rbind(species_results, iteration_df)
    }  
  }
  null_df_list[[k]] <- species_results
}  

null_dflist<-do.call(rbind, null_df_list)

mean_null <- null_dflist %>%
  bind_rows() %>%
  group_by(Species, model, metric) %>%
  summarize(mean_value = mean(values))


saveRDS(mean_null, "E:/Madariaga/Documents/Phd/Manuscript2/SSN/NULL/null_cv_df.rds")



##### Train and test in independent datasets


glm_null_independent <- list()

for (name in Knames) {
     tryCatch(
      {
        # Construct the column name with .ev suffix
        ssn_obj_subset <- ssn_objects[[name]]
        
        # Get the ssn object from the ssn_list
        print("check ssn_obj_df")
        
        frml <- null_AIC_list[[name]][["formula"]]%>%as.character()
        frml <- as.formula(paste(paste0(frml[2],".ev"), frml[1], frml[3], collapse = " "))
        
        
        # Run the glmssn model with the .ev column
        model <- glmssn(
          formula = frml,
          ssn.object = ssn_obj_subset,
          family = "binomial",
          CorModels = NULL
        )
        
        # Save the model in the glm_null_models list with the name
        glm_null_independent[[name]] <- model},  
      error =function(e) {
        # Store the error message in the list with the name from Knames
        glm_null_independent[[name]] <<- paste("Failed model:", name, conditionMessage(e))
      }
    )
  }


saveRDS(glm_null_independent,"E:/Madariaga/Documents/Phd/Manuscript2/glm_null_independent.rds")

source("E:/Madariaga/Documents/Phd/SSN/functions_modified.R")

eval_results_null_independent <- list()

# Loop through each name in Knames
for (name in Knames) {
    tryCatch(
      {
        # Get the corresponding glmssn model from glm_null_models
        model <- glm_null_independent[[name]]
        
        # Create the "_MissingObs_" evaluation dataframe
        Evaluation <- getSSNdata.frame(model, "_MissingObs_")
        
        # Create Eval_obj using the evaluation dataframe and ssn_obj0
        Eval_obj <- putSSNdata.frame(Evaluation, model, "_MissingObs_")
        
        # Predict using ssn_obj0 and the "_MissingObs_" evaluation dataframe
        Predeval <- predict(model, "_MissingObs_")
        
        # Calculate SSNProb for Prediceval
        SSN0_predev <- SSNProb(Predeval)
        
        # Subset pointdata and select relevant columns for evaluation
        specific_df <- glm_null_independent[[name]][["ssn.object"]]@obspoints@SSNPoints[[1]]@point.data
        SSN0presabs <- specific_df[is.na(specific_df[[paste0(name,".ev")]]), c("id", name)] 
        
        # Combine SSN0presabs with SSN0_predev
        idssn0presabs <- data.frame(SSN0presabs, SSN0_predev)
        
        # Perform the evaluation using eval_mod on idssn0presabs
        result <- eval_mod(idssn0presabs)
        
        # Store the evaluation result in the list with the name from Knames
        eval_results_null_independent[[name]] <- result
      },  
      error =function(e) {
        # Store the error message in the list with the name from Knames
        eval_results_null_independent[[name]] <<- paste("Failed model:", name, conditionMessage(e))
      }
    )
  }

eval_null_independent_list <- list()

for (k in Knames){
  species_results <- data.frame()
      if (is.list(eval_results_null_independent[[k]])){
      iteration_df<- data.frame(
        Species= k,
        model= "glm",
        metric= eval_results_null_independent[[k]][["performance"]][["Measures"]],
        values = eval_results_null_independent[[k]][["performance"]][["Value"]]
      )
      species_results <- rbind(species_results, iteration_df)
    }  
 
  eval_null_independent_list[[k]] <- species_results
} 


eval_null_independent_dflist <- eval_null_independent_list%>%
  bind_rows() %>%
  group_by(Species, model, metric) %>%
  summarize(mean_value = mean(values))


saveRDS(eval_null_independent_dflist, "E:/Madariaga/Documents/Phd/Manuscript2/glm_null_independent_df.rds")
