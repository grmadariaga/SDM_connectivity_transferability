library(tidymodels)
library(SSN)
library(sirad)
library(dplyr)
library(parallel)
library(doParallel)
library(foreach)
library(caret)


###Merge SSN object from STARS to dataframe on environmental variables

###Load SSN object Kinzig

K <- importSSN("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/Kinzig_upd/Kin_upd.ssn", predpts = "Kinzig_loc")


#createDistMat(K, predpts="Kinzig_loc",  o.write = TRUE)
distPred1km <- getStreamDistMat(K, Name = "Kinzig_loc")

dmats <- getStreamDistMat(K)


K<-additive.function(K, "H2OArea", "computed.afv")

names(K)

###Get observation dataframe

obs_df <- getSSNdata.frame(K,"Obs")


###Load Kinzig IHA dataset

Kinzig_pa <- read.csv("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/Kinzig_upd/Kinzig_large_iha.csv")[,-1]
filtered_cols <- Kinzig_pa %>%
  select(7:98) %>%
  select_if(~ sum(.) >= 15)
Kinzig_df <- Kinzig_pa %>%
  select(1:6, 99:104, all_of(names(filtered_cols)))



###Load the Kinzig bioclim dataset

Biovar <- read.csv("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/Kinzig_upd/BiovarK_upd.csv", check.names = FALSE)[,-1]
Biovar_cl <- Biovar[!duplicated(Biovar),]
##Join with obs_df

obs_df <- left_join(obs_df,
                    Kinzig_df[,c("site_code", "year", "Channel",names(filtered_cols), "dh4", "fl1", "ra4", "th3")], 
                    by= c("site_code", "year", "Channel"), relationship= "many-to-many") %>%
  filter(duplicated(id) == FALSE)%>%
  left_join(.,Biovar_cl, by= c("site_code", "year"), relationship= "many-to-many")


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

###Extract species names to create separate dataframes (evaluation sets)

Knames <- names(filtered_cols)


####Extract names to an object to avoid for rerunning this process

#saveRDS(Knames, "E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/Knames.rds")


# Create an empty list to store the SSN objects
ssn_list <- list()

set.seed(432)

# Create 5 balanced folds for each species according to specific prevalence, this folds are going to be used to train and test SDMs
for (name in Knames) {
  
  subset <- obs_df[, c(names(obs_df[,c(1:23,114:121)]), name)]
  subset_folded <- vfold_cv(subset,v=5, repeats = 1,
                            strata= name)
  
  
  ssn_obj <- putSSNdata.frame(subset, K)  
  
  # Assign the SSN object to the list with the corresponding name
  ssn_list[[name]][["obj"]] <- ssn_obj
  ssn_list[[name]][["folds"]]<-subset_folded
}



cv_ids <- list()

for (k in Knames){
  folds <- ssn_list[[k]][["folds"]]
  cv_ids[[k]] <- folds
}


####Extract cv_ids to an object

#saveRDS(cv_ids,"E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/cv_ids.rds")
#saveRDS(cv_ids,"E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/ssn_list.rds")



### Create object list for training in a dataset and evaluating in another
Knames_ev <- paste0(Knames, ".ev")

obs_ev <- obs_df

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


####SaveRDS(ssn_objects, "E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/ssn_objects.rds")

##### Now we select the best variables for each species
####Best model selection for null models

var_names <- colnames(obs_df[,114:121])

null_AIC_list <- list()

nCores <- detectCores() - 2
cluster <- makeCluster(nCores, type= "PSOCK")
doParallel::registerDoParallel(cl = cluster)

# Loop through each response variable
foreach(i= Knames, .combine = cbind) %do% {
  
  ssn_i <- ssn_objects[[i]]
  df_ssn_i<<- getSSNdata.frame(ssn_i,"Obs")
  
  # Initialize variables to keep track of the best model and its AIC
  best_model_formula <- NULL
  lowest_aic <- Inf
  
  # Create the formula for the current combination
  formula <- as.formula(paste(i, "~", paste(var_names, collapse = " + ")))
  
  # Fit the model
  model_1 <- try(glm(formula = formula, 
                     data = df_ssn_i,
                     family = binomial,
                     na.action = "na.fail")) 
  
  print("here")
  ###create cluster for dredge
  
  clusterExport(cluster, "df_ssn_i")
  
  ListModelNsp <- MuMIn::dredge(model_1, m.lim = c(0,5))
  
  BestModelNsp <- eval(attributes(ListModelNsp)$model.calls[[1]])
  
  forml <- BestModelNsp$formula
  
  
  
  # Save the best model formula with the lowest AIC for the current response variable
  null_AIC_list[[i]] <- list("formula" = forml)
}

stopCluster(cluster)


#saveRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/Kinzig_upd/server/Kinzig_upd/Results/dredge_formulasK_upd.rds")