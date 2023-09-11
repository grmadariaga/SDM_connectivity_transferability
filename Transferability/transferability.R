library(dplyr)
library(tidyverse)
library(rstatix)
library(rcompanion)



###Load maxent scores with CV  and independent validations

maxent_cv <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/MAXENT/mean_maxent_cv_df.rds")%>%
mutate(values= mean_value)%>%
  select(-mean_value)
maxent_independent <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/MAXENT/maxent_independent_df.rds")

null_independent <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/glm_null_independent_df.rds")

null_cv <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN/null_cv_df.rds")%>%
  mutate(values= mean_value)%>%
  select(-mean_value)

ssn_independent <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/glmssn_independenteval.rds")%>%
  mutate(model= "glmssn_independent")

ssn_cv <-readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN/SSN_sp/mean_ssn_cv_df.rds")%>%
mutate(values= mean_value)%>%
  select(-mean_value)

rf_cv <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/RF/rf_cv.rds")%>%
  mutate(values= mean_value)%>%
  select(-mean_value)

rf_independent <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/RF/rf_indpndnt.rds")%>%
  mutate(model="RF_independent")

Spmodels <- readRDS("E:/Madariaga/Documents/Phd/Manuscript2/SDM_connectivity_transferability/SSN_o/Kinzig_upd/Server/Kinzig_upd/Results/ssn_selection_90.rds")

sp_df <- df <- do.call(rbind, lapply(names(Spmodels), function(name) {
  data.frame(Species = name, model = Spmodels[[name]], stringsAsFactors = FALSE)
}))


model_counts <- table(sp_df$model)%>%as.data.frame()


model_counts <- model_counts%>%
  arrange(
    case_when(
      grepl("Euclid$", Var1) ~ 1,
      grepl("tailup$", Var1) ~ 2,
      grepl("taildown$", Var1) ~ 3,
      TRUE ~ 4
    )
  )
  

x_names <-model_counts$Var1%>%as.character()
x_names[15]<- "Failed"


model_counts$FillColor <- c(rep("blue",4), rep("cornflowerblue",5), rep("aquamarine4",5), "gray")

bar_colors <- c(rep("blue",4), rep("cornflowerblue",5), rep("aquamarine4",5), "gray")

ggplot(model_counts,
       aes(x = factor(Var1, levels= model_counts$Var1%>%as.character()%>%as.factor()),
           y = Freq, fill= FillColor)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = x_names)+
  scale_fill_identity()+
  labs(x = "Model", y = "Species using the model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Frequency of usage of spatial model")



All_scores <- rbind(rf_cv,maxent_cv,null_cv,ssn_cv,
                    rf_independent,maxent_independent,
                    null_independent, ssn_independent
                    )

ggplot(All_scores, aes(x = model, y = values, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(x = "Model", y = "Value", title = "Performance Measures by Model")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggplot(All_scores, aes(x = values)) + 
  geom_density(aes(fill = model), alpha = 0.5) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = 'Performance', y = 'Density') + 
  scale_fill_discrete(name = 'Model')

all_filt <- All_scores%>%
  group_by(Species) %>%
  filter(all(c(All_scores$model%>%unique()) %in% model))

ggplot(all_filt, aes(x = model, y = values, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(x = "Model", y = "Value", title = "Performance Measures by Model")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


independents <- all_filt[grepl("independent", all_filt$model, ignore.case = TRUE), ]
cv_datasets <- all_filt[grepl("_cv", all_filt$model, ignore.case = TRUE), ]


independent_large <- All_scores[grepl("independent", All_scores$model, ignore.case = TRUE), ]
cv_large <- All_scores[grepl("_cv", All_scores$model, ignore.case = TRUE), ]


AUC.kw <- dunn_test(values ~ model,
                    p.adjust.method = "bonferroni",
                    data= independents%>%ungroup()%>%
                      filter(., metric == "AUC_ROC"))
AUC.kw


sTSS.kw <- dunn_test(values ~ model,
                    p.adjust.method = "bonferroni",
                    data= independents%>%ungroup()%>%
                      filter(., metric == "sTSS"))
sTSS.kw


cv_kw <-dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= cv_datasets%>%ungroup()%>%
            filter(., metric == "AUC_ROC"))

dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= cv_datasets%>%ungroup()%>%
            filter(., metric == "sTSS"))


dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= independents%>%ungroup()%>%
            filter(., metric == "AUC_ROC"))


dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= independents%>%ungroup()%>%
            filter(., metric == "sTSS"))

dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= all_filt%>%ungroup()%>%
            filter(., metric == "sTSS"))

dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= all_filt%>%ungroup()%>%
            filter(., metric == "AUC_ROC"))

dunn_test(values ~ model,
          p.adjust.method = "bonferroni",
          data= cv_filt%>%ungroup()%>%
            filter(., metric == "sTSS"))


ggplot(independents, aes(x = values)) + 
  geom_density(aes(fill = model), alpha = 0.5) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = 'Performance', y = 'Density') + 
  scale_fill_discrete(name = 'Model')


ggplot(cv_datasets, aes(x = values)) + 
  geom_density(aes(fill = model), alpha = 0.5) +
  facet_wrap(~ metric, scales = "free") +
  labs(x = 'Performance', y = 'Density') + 
  scale_fill_discrete(name = 'Model')

ggplot(independents, aes(x = model, y = values, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(x = "Model", y = "Value", title = "Performance Measures by Model")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(cv_datasets, aes(x = model, y = values, fill = model)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  labs(x = "Model", y = "Value", title = "Performance Measures by Model")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
