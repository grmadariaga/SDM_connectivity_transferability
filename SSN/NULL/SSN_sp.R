library(tidymodels)
library(SSN)
library(sirad)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(caret)


###Load SSN object Kinzig

K <- importSSN("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/Kin_upd.ssn", predpts = "Kinzig_loc")


#createDistMat(K, predpts="Kinzig_loc",  o.write = TRUE)
distPred1km <- getStreamDistMat(K, Name = "Kinzig_loc")

dmats <- getStreamDistMat(K)


K<-additive.function(K, "H2OArea", "computed.afv")

names(K)

###Get observation dataframe

obs_df <- getSSNdata.frame(K,"Obs")


###Load Kinzig dataset

Kinzig_pa <- read.csv("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/Kinzig_large_iha.csv")[,-1]
filtered_cols <- Kinzig_pa %>%
  select(7:98) %>%
  select_if(~ sum(.) >= 15)
Kinzig_df <- Kinzig_pa %>%
  select(1:6, 99:104, all_of(names(filtered_cols)))



###Load the bioclim dataset

Biovar <- read.csv("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/BiovarK_upd.csv", check.names = FALSE)[,-1]
Biovar_cl <- Biovar[!duplicated(Biovar),]
##Join with obs_df

obs_df <- left_join(obs_df,
                    Kinzig_df[,c("site_code", "year", "Channel",names(filtered_cols), "dh4", "fl1", "ra4", "th3")], 
                    by= c("site_code", "year", "Channel")) %>%
  filter(duplicated(id) == FALSE)%>%
  left_join(.,Biovar_cl, by= c("site_code", "year"))


###Check variable colinearity

vars <- obs_df[,114:136]

cor_matrix <- cor(vars)

threshold <- 0.65

highly_correlated <- which(abs(cor_matrix) > threshold & cor_matrix != 1, arr.ind = TRUE)
var_to_drop <- c()

for (i in 1:nrow(highly_correlated)) {
  row_index <- highly_correlated[i, "row"]
  col_index <- highly_correlated[i, "col"]
  var_to_drop <- c(var_to_drop, colnames(vars)[col_index])
}

obs_df <- obs_df[, !(colnames(obs_df) %in% var_to_drop)]

obs_ev <- obs_df

###Extract species names to create separate dataframes (evaluation sets)

Knames <- names(filtered_cols)


null_AIC_list <-readRDS("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/server/Kinzig_upd/Results/dredge_formulasK_upd.rds")
CorModels <-readRDS("E:/Madariaga/Documents/Phd/SSN/Kinzig_upd/server/Kinzig_upd/Results/ssn_selection_90.rds")


# Create an empty list to store the SSN objects
ssn_list <- list()

set.seed(432)

# Loop through each name in Knames
for (name in Knames) {
  # Subset the obs_df dataframe
  subset <- obs_df[, c(names(obs_df[,c(1:23,114:121)]), name)]
  subset_folded <- vfold_cv(subset,v=5, repeats = 1,
                            strata= name)
  
  # Create a new SSN object with the subset
  ssn_obj <- putSSNdata.frame(subset, K)  # Replace K with the appropriate value
  
  # Assign the SSN object to the list with the corresponding name
  ssn_list[[name]][["obj"]] <- ssn_obj
  ssn_list[[name]][["folds"]]<-subset_folded
}



options(na.action = "na.omit")


#ssn_obj0 <- glmssn(Asellus.aquaticus.ev ~ ra6, ssn.object = ssn_list$Asellus.aquaticus,
#                  family = "binomial",
#                        CorModels = NULL)



glm_ssn_models <- list()

# Loop through each name in Knames
for (name in Knames[19:90]) {
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
        cmdls <- CorModels[[name]]%>%as.character()
        
        # Run the glmssn model
        model <- glmssn(
          formula = frml,
          ssn.object = ssn_obj_subset,
          family = "binomial",
          addfunccol = "computed.afv",
          CorModels = cmdls
        )
        
        # Save the model in the glm_null_models list with the name
        glm_ssn_models[[name]][[Fold]] <- model},  
      error =function(e) {
        # Store the error message in the list with the name from Knames
        glm_ssn_models[[name]][[Fold]] <<- paste("Failed model:", name, conditionMessage(e))
      }
    )
  }
}

saveRDS(glm_ssn_models,"E:/Madariaga/Documents/Phd/Manuscript2/SSN/SSN_sp/ssn_models.rds" )
glm_ssn_models <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SSN/SSN_sp/ssn_models.rds")
source("E:/Madariaga/Documents/Phd/SSN/functions_modified.R")

eval_results_ssn <- list()

# Loop through each name in Knames
for (name in Knames) {
  for (i in 1:5) {
    
    
    tryCatch(
      {
        Fold <- paste0("Fold.", i)
        # Get the corresponding glmssn model from glm_null_models
        model <- glm_ssn_models[[name]][[Fold]]
        
        # Create the "_MissingObs_" evaluation dataframe
        Evaluation <- getSSNdata.frame(model, "_MissingObs_")
        
        # Create Eval_obj using the evaluation dataframe and ssn_obj0
        Eval_obj <- putSSNdata.frame(Evaluation, model, "_MissingObs_")
        
        # Predict using ssn_obj0 and the "_MissingObs_" evaluation dataframe
        Predeval <- predict(model, "_MissingObs_")
        
        # Calculate SSNProb for Prediceval
        SSN0_predev <- SSNProb(Predeval)
        
        # Subset pointdata and select relevant columns for evaluation
        specific_df <- glm_ssn_models[[name]][[Fold]][["ssn.object"]]@obspoints@SSNPoints[[1]]@point.data
        SSN0presabs <- specific_df[is.na(specific_df[[paste0(name,Fold)]]), c("id", name)] 
        
        # Combine SSN0presabs with SSN0_predev
        idssn0presabs <- data.frame(SSN0presabs, SSN0_predev)
        
        # Perform the evaluation using eval_mod on idssn0presabs
        result <- eval_mod(idssn0presabs)
        
        # Store the evaluation result in the list with the name from Knames
        eval_results_ssn[[name]][[Fold]] <- result
      },  
      error =function(e) {
        # Store the error message in the list with the name from Knames
        eval_results_ssn[[name]][[Fold]] <<- paste("Failed model:", name, conditionMessage(e))
      }
    )
  }
}

ssn_df_list <- list()

for (k in Knames){
  species_results <- data.frame()
  for (i in 1:5){
    Fold <- paste0("Fold.", i)
    if (is.list(eval_results_ssn[[k]][[Fold]])){
      iteration_df<- data.frame(
        Species= k,
        iteration= i,
        model= "glmssn_cv",
        metric= eval_results_ssn[[k]][[Fold]][["performance"]][["Measures"]],
        values = eval_results_ssn[[k]][[Fold]][["performance"]][["Value"]]
      )
      species_results <- rbind(species_results, iteration_df)
    }  
  }
  ssn_df_list[[k]] <- species_results
}  

ssn_dflist<-do.call(rbind, ssn_df_list)

mean_ssn_cv <- ssn_dflist %>%
  bind_rows() %>%
  group_by(Species, model, metric) %>%
  summarize(mean_value = mean(values))



saveRDS(mean_ssn_cv,"E:/Madariaga/Documents/Phd/Manuscript2/SSN/SSN_sp/mean_ssn_cv_df.rds" )



##### Train and test in independent datasets
###Extract species names to create separate dataframes (evaluation sets)

obs_df <- obs_df[, !(colnames(obs_df) %in% var_to_drop)]

obs_ev <- obs_df

Knames <- names(filtered_cols)
Knames_ev <- paste0(Knames, ".ev")

obs_ev <- obs_ev%>%
  rename_with(~ paste0(., ".ev"), any_of(Knames))

for (i in Knames_ev){
  obs_ev[,i]<- ifelse(obs_ev$dataset=="Senckenberg", NA,obs_ev[,i])
}


ssn_objects <- list()

# Loop through each name in Knames
for (name in Knames) {
  # Subset the obs_df dataframe
  subset <- obs_df[, c(names(obs_df[,c(1:23,114:121)]), name)]%>%
    left_join(.,obs_ev[,c("id", paste0(name,".ev"))], by= "id")
  
  
  # Create a new SSN object with the subset
  ssn_obj <- putSSNdata.frame(subset, K)  # Replace K with the appropriate value
  
  # Assign the SSN object to the list with the corresponding name
  ssn_objects[[name]] <- ssn_obj
}




glm_ssn_independent <- list()

for (name in Knames) {
  tryCatch(
    {
      # Construct the column name with .ev suffix
      ssn_obj_subset <- ssn_objects[[name]]
      
      # Get the ssn object from the ssn_list
      print("check ssn_obj_df")
      
      frml <- null_AIC_list[[name]][["formula"]]%>%as.character()
      frml <- as.formula(paste(paste0(frml[2],".ev"), frml[1], frml[3], collapse = " "))
      crmdl <- CorModels[[name]]%>%as.character()
      
      # Run the glmssn model with the .ev column
      model <- glmssn(
        formula = frml,
        ssn.object = ssn_obj_subset,
        family = "binomial",
        addfunccol = "computed.afv",
        CorModels = crmdl
      )
      
      # Save the model in the glm_null_models list with the name
      glm_ssn_independent[[name]] <- model},  
    error =function(e) {
      # Store the error message in the list with the name from Knames
      glm_ssn_independent[[name]] <<- paste("Failed model:", name, conditionMessage(e))
    }
  )
}


source("E:/Madariaga/Documents/Phd/SSN/functions_modified.R")

eval_results_ssn_independent <- list()

# Loop through each name in Knames
for (name in Knames) {
  tryCatch(
    {
      # Get the corresponding glmssn model from glm_null_models
      model <- glm_ssn_independent[[name]]
      
      # Create the "_MissingObs_" evaluation dataframe
      Evaluation <- getSSNdata.frame(model, "_MissingObs_")
      
      # Create Eval_obj using the evaluation dataframe and ssn_obj0
      Eval_obj <- putSSNdata.frame(Evaluation, model, "_MissingObs_")
      
      # Predict using ssn_obj0 and the "_MissingObs_" evaluation dataframe
      Predeval <- predict(model, "_MissingObs_")
      
      # Calculate SSNProb for Prediceval
      SSN0_predev <- SSNProb(Predeval)
      
      # Subset pointdata and select relevant columns for evaluation
      specific_df <- glm_ssn_independent[[name]][["ssn.object"]]@obspoints@SSNPoints[[1]]@point.data
      SSN0presabs <- specific_df[is.na(specific_df[[paste0(name,".ev")]]), c("id", name)] 
      
      # Combine SSN0presabs with SSN0_predev
      idssn0presabs <- data.frame(SSN0presabs, SSN0_predev)
      
      # Perform the evaluation using eval_mod on idssn0presabs
      result <- eval_mod(idssn0presabs)
      
      # Store the evaluation result in the list with the name from Knames
      eval_results_ssn_independent[[name]] <- result
    },  
    error =function(e) {
      # Store the error message in the list with the name from Knames
      eval_results_ssn_independent[[name]] <<- paste("Failed model:", name, conditionMessage(e))
    }
  )
}


eval_ssn_independent_list <- list()

for (k in Knames){
  species_results_ssn <- data.frame()
  if (is.list(eval_results_ssn_independent[[k]])){
    iteration_df_ssn<- data.frame(
      Species= k,
      model= "glm",
      metric= eval_results_ssn_independent[[k]][["performance"]][["Measures"]],
      values = eval_results_ssn_independent[[k]][["performance"]][["Value"]]
    )
    species_results_ssn <- rbind(species_results_ssn, iteration_df_ssn)
  }  
  
  eval_ssn_independent_list[[k]] <- species_results_ssn
} 


eval_ssn_independent_df <- eval_ssn_independent_list%>%
  bind_rows()


saveRDS(glm_ssn_independent, "E:/Madariaga/Documents/Phd/Manuscript2/glmssn_independentmodels.rds")
saveRDS(eval_ssn_independent_df, "E:/Madariaga/Documents/Phd/Manuscript2/glmssn_independenteval.rds")
