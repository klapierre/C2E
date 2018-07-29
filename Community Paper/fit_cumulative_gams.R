################################################################################
##  fit_cumulative_gams.R: Script that fits GAMs to the cumulative RAC metrics
##  data. The goal is to compare models with and without a treatment effect.
##
##  Author: Andrew Tredennick (atredenn@gmail.com)
##  Date created: March 21, 2018
################################################################################

#########
### TODO: Make sure model structure is correct (i.e., do your homework)
### TODO: Add separate analysis where multivariate change is the response
########

# NOTES:
#  (1) for deviance and AIC, negative deltas indicate better models
#  (2) the p-value then says whether the deviance difference is significant
#  (3) some p-values will be NA -- this is OK and indicates the full the model is
#      DEFINITELY NOT BETTER than the null model. So, think of NA as p>0.05.


##  Clear the workspace
rm(list = ls(all.names = TRUE))



####
####  LOAD LIBRARIES -----------------------------------------------------------
####
library(tidyverse)
library(ggthemes)
library(mgcv)



####
####  SET WORKING DIRECTORIES AND FILENAMES ------------------------------------
####
work_dir  <- "~/Repos/C2E/Community Paper/" # change as needed
data_dir  <- "~/Dropbox/C2E/Products/CommunityChange/March2018 WG/"
results_dir <- "~/Dropbox/C2E/Products/CommunityChange/Summer2018_Results/"
data_file <- "CORRE_RACS_Subset_Perm.csv" # change as needed
setwd(work_dir)



####
####  DEFINE MODEL FITTING FUNCTION --------------------------------------------
####

fit_compare_gamms <- function(df, response){
  # Fits two GAMMs and compares them with AIC and LLR
  #
  # Args:
  #  data: a dataframe with necessary columns for fitting the GAMMs
  #  response: name of the response variable, must be a column in the dataframe
  #
  # Returns:
  #  A tibble with LLR delta deviance, LLR p-value, and delta AIC
  
  test_formula <- as.formula(
    paste(response, 
          "~ s(treatment_year, treatment, bs = 'fs', k = (num_years-1)) + 
          s(plot_id, bs='re')"
    )
  )
  
  null_formula <- as.formula(
    paste(response, 
          "~ s(treatment_year, bs = 'fs', k = (num_years-1)) + 
          s(plot_id, bs='re')"
    )
  )
  
  gam_test <- gam(
    test_formula, 
    data = df, 
    method = "REML"
  )
  
  gam_null <- gam(
    null_formula,
    data = df, 
    method = "REML"
  )
  
  # LLR tests
  pvalue <- anova(gam_null, gam_test, test="Chisq")$`Pr(>Chi)`[2]
  dev <- anova(gam_null, gam_test, test="Chisq")$`Resid. Dev`
  delta_div <- diff(dev)  # full - null
  
  # AIC tests
  aics <- AIC(gam_null, gam_test)$AIC
  delta_aic <- diff(aics)  # full - null
  
  return(
    tibble(
      response_var = response,
      p_value = pvalue,
      delta_deviance = delta_div,
      delta_aic = delta_aic
    )
  )
}  # end of model fit and comparison function



####
####  DEFINE FUNCTION TO FILL TIBBLE WHEN YEARS < 4
####

fill_empties <- function(...){
  return(
    tibble(
      response_var = c(
        "richness_change_abs", 
        "evenness_change_abs", 
        "rank_change", 
        "gains", 
        "losses"
      ),
      p_value = NA,
      delta_deviance = NA,
      delta_aic = NA
    )
  )
}



####
####  READ IN DATA AND CALCULATE CUMULATIVE CHANGE -----------------------------
####
change_metrics <- read.csv(paste0(data_dir,data_file))

##  Calculate cumulative sums of each metric (from Kevin)
change_cumsum <- change_metrics %>%
  group_by(site_project_comm, treatment, plot_id) %>%
  mutate(richness_change_abs = abs(richness_change)) %>%
  mutate(evenness_change_abs = abs(evenness_change)) %>%
  mutate_at(vars(richness_change, 
                 richness_change_abs, 
                 evenness_change,
                 evenness_change_abs, 
                 rank_change, 
                 gains, 
                 losses), 
            funs(cumsum)) %>%
  mutate(control = ifelse(plot_mani==0,"control","treatment"))



####
####  LOOP OVER SITE_PROJECT_COMMS AND COMPARE CONTROLS VS. TREATMENTS ---------
####
all_sites <- unique(change_cumsum$site_project_comm)
all_comparisons <- {} # empty object for storage

for(do_site in all_sites){
  site_data <- filter(change_cumsum, site_project_comm == do_site)
  site_controls <- filter(site_data, plot_mani == 0)
  site_treatments <- filter(site_data, plot_mani != 0)
  all_treatments <- unique(site_treatments$treatment)
  
  for(do_treatment in all_treatments){
    treatment_data <- filter(site_treatments, treatment == do_treatment)
    model_data <- rbind(site_controls, treatment_data)
    num_years <- length(unique(model_data$treatment_year))
    
    # Skip data with less than four years
    if(num_years < 4){
      tmp_out <- fill_empties() %>%
        mutate(
          site_proj_comm = do_site,
          treatment = do_treatment
        ) %>%
        dplyr::select(
          site_proj_comm,
          treatment,
          response_var,
          p_value,
          delta_deviance,
          delta_aic
        )
    }
    
    # Compare models for data with more than four years
    if(num_years > 3){
      
      # Richness
      rich_test <- fit_compare_gamms(
        df = model_data,
        response = "richness_change_abs"
      )
      
      # Evennes
      even_test <- fit_compare_gamms(
        df = model_data,
        response = "evenness_change_abs"
      )
      
      # Rank change
      rank_test <- fit_compare_gamms(
        df = model_data,
        response = "rank_change"
      )
      
      # Gains
      gain_test <- fit_compare_gamms(
        df = model_data,
        response = "gains"
      )
      
      # Losses
      loss_test <- fit_compare_gamms(
        df = model_data,
        response = "losses"
      )
      
      tmp_out <- bind_rows(
        rich_test,
        even_test,
        rank_test,
        gain_test,
        loss_test
      ) %>%
        mutate(
          site_proj_comm = do_site,
          treatment = do_treatment
        ) %>%
        dplyr::select(
          site_proj_comm,
          treatment,
          response_var,
          p_value,
          delta_deviance,
          delta_aic
        )
      
    } # end if for num_years
    
    all_comparisons <- all_comparisons %>%
      bind_rows(tmp_out)
      
  } # end treatment loop
  
  print(paste("Done with site:", do_site))
  
} # end site loop



####
####  SAVE DELTA_AIC TABLE -----------------------------------------------------
####
save_comparisons <- all_comparisons %>%
  filter(is.na(delta_deviance) == FALSE) %>%
  mutate(
    sig_diff_cntrl_trt = ifelse(
      p_value <= 0.05 & sign(delta_deviance) == -1,
      "yes",
      "no"
    )
  ) %>%
  mutate(
    sig_diff_cntrl_trt = ifelse(is.na(sig_diff_cntrl_trt) == TRUE, "no", sig_diff_cntrl_trt)
  )

write_csv(
  x = save_comparisons, 
  path = paste0(results_dir, "gam_comparison_table.csv")
)



####
####  VISUALIZE THE DLETA AIC TABLE --------------------------------------------
####
# delta_aics <- read.csv(paste0(data_dir,"gam_delta_aic_table.csv"), 
#                        row.names = 1) %>%
#   gather(key = metric, value = delta_aic, rich_delta_aic:loss_delta_aic) %>%
#   mutate(site_treatment = paste(site_project_comm, treatment, sep = "::"),
#          different = ifelse(delta_aic < -10, "yes", "no"))
# 
# ggplot(delta_aics, aes(y = site_treatment, x = metric))+
#   geom_tile(aes(fill = different))+
#   scale_fill_brewer(type = "qual", 
#                     labels = c("C and T not different", "C and T different", "NA"), 
#                     name = NULL)+
#   scale_x_discrete(labels = c("Evenness","Gains","Losses","Rank Change", "Richness"))+
#   xlab("Metric")+
#   ylab("Site and Treatment")
# ggsave(filename = paste0(data_dir,"figures/delta_aic_figure.pdf"),
#        height = 14, 
#        width = 7, 
#        units = "in")


####
####  TEST FIT USING JUST PPLOTS -----------------------------------------------
####
##  Subset controls and one treatment for pplots
# do_site <- "CDR_e001_D"
# do_treatment <- 6
# pplots_cumsum  <- filter(change_cumsum, site_project_comm == do_site)
# pplot_controls <- filter(pplots_cumsum, treatment == 9 | treatment == do_treatment)
# 
# ##  Fit nested GAMs
# ### TODO: CHECK THESE MODEL STRUCTURES!!! Email Gavin?
# gam_test <- gam(evenness_change_abs ~ s(treatment_year, treatment, bs = "fs") + 
#                   s(plot_id, bs="re"), data = pplot_controls)
# gam_null <- gam(evenness_change_abs ~ s(treatment_year) + s(plot_id, bs="re"),
#                 data = pplot_controls)
# 
# ##  Compare nested models
# AIC(gam_test, gam_null) # less conservative
# BIC(gam_test, gam_null) # more conservative
# deltaAIC <- diff(AIC(gam_test, gam_null)$AIC)
# 
# ##  Make example plot of predictions from the complex model
# newdata <- data.frame(treatment_year = rep(0:11,2),
#                       treatment = rep(c("N1P0","N2P0"), each = 12),
#                       plot_id = rep(factor(27),times = 12*2))
# gam_preds <- predict(gam_test, 
#                      newdata = newdata, 
#                      type = "response", 
#                      exclude = "s(plot_id)")
# plot_preds <- data.frame(predictions = gam_preds,
#                          treatment_year = rep(0:11,2),
#                          treatment = rep(c("N1P0","N2P0"), each = 12))
# plot_data <- pplot_controls %>%
#   ungroup() %>%
#   select(treatment_year, rank_change, treatment)
# 
# ggplot()+
#   geom_point(data = plot_data, 
#              aes(x = treatment_year, y = rank_change, color = treatment))+
#   geom_line(data = plot_preds, 
#             aes(x = treatment_year, y = predictions, color = treatment),
#             size = 1.2)+
#   scale_color_brewer(type = "qual", name = "Treatment")+
#   xlab("Treatment year")+
#   ylab("Rank change")+
#   ggtitle("KNZ PPlots Example", subtitle = "Control (N1P0) vs. N2P0")
# ggsave(filename = paste0(data_dir,"figures/pplot_gam_example.pdf"),
#        height = 4,
#        width = 5,
#        units = "in")
# 
