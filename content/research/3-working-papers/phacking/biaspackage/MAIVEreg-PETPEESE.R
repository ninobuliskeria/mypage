library("tidyverse")
library("writexl")
library("readxl")

rm(list = ls())

# Set Working Directory and Import Data
#setwd("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/GitHub/maive")
#setwd("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/DATA")
  # general <- read_excel("H:/My Drive/BIAS/DATA/MAIVEreg.xlsx") 

general <- read_excel("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/DATA/MAIVEreg.xlsx")
all_objects <- ls()
rm(list = all_objects[all_objects != "general"])



                                                      # POOLED EFFECT ####
# gen <- general %>%
  # select("E1", "SE1", "N1", "E2", "SE2", "N2", "E5", "SE5", "N5","studyID", "metaID", "studyPublishD", "np" )
#select("bs", "sebs", "Ns", "studyid", "meta_id", "studyPublishD", "np")

gen <- general %>%
  select("E1", "SE1", "N1","studyID", "metaID", "studyPublishD", "np" )

# Unique meta_ids
numbers <- gen$metaID %>% unique()

object<-c("E","SE","F","Hausman","CV_of_Chi2", "Obs")
# Function to Split and Name Data
split_and_name <- function(data, condition, prefix) {
  split_data <- split(data[condition, ], data$meta_id[condition])
  split_data <- lapply(split_data, function(df) {
    df$meta_id <- NULL
    return(df)
  })
  names(split_data) <- paste0(prefix, names(split_data))
  return(split_data)
}

# FUNCTION TO HANDLE THE REPETITIVE TASKS 
run_analysis <- function(list_of_meta, y, w) {
  source("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/GitHub/maive/maive-PETPEESE/maivefunction.R")
  # OPTIONS: 
  method <- 3 # method: PET:1, PEESE:2, PET-PEESE:3, EK:4 (default 3)
  weight <- 0 #no weight: 0; weights: 1, adjusted weights: 2 
  instrument <- 1 # instrumenting 
  studylevel <-0  # none: 0, fixed effects: 1, cluster: 2
  # default options are method=3; weight=0; instrument=1; studylevel=0 
  
  MAIVEresults <- data.frame(object)
  for (i in numbers){
    dat <- list_of_meta[[paste0("Meta_", i)]]
    MAIVE <- tryCatch(
      {
        maive(dat=dat, method=method, weight=weight, instrument=instrument, studylevel=studylevel)
      },
      error = function(e) {
        return(NULL)
      }
    )
    # If MAIVE is not NULL (meaning no error), extract results and save them
    if (!is.null(MAIVE)) {
      MAIVEresults[[paste0("Meta_", i)]] <- c(MAIVE$beta, MAIVE$SE, MAIVE$`F-test`, MAIVE$Hausman, MAIVE$Chi2, nrow(dat))
    } else {
      MAIVEresults[[paste0("Meta_", i)]] <- c(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    }
    print(i) # To monitor the progress
  }
  # If you want to view the results after the loop
  print(MAIVEresults)
  
  
  results <- as.data.frame(t(MAIVEresults))
  colnames(results) <- results[1, ]
  results <- results[-1, ] %>%
    mutate(across(everything(), as.numeric))
  
  # Add row names as a new column in the data frame
  results_with_row_names <- results %>%
    tibble::add_column(Row_Names = row.names(results), .before = 1)
  
  # Save the modified data frame to an Excel file
  write_xlsx(results_with_row_names, paste0("MAIVE_", y, "pooled_" , w, ".xlsx"))
}

## WINSORIZATION at 1% ####
w <- 1
# Renaming and Preprocessing Data
df <- general %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E1", "SE1", "N1", "studyID", "metaID")) #, "n_E_m")) 
# Splitting Data Based on Conditions
assign(paste0("all_", w), split_and_name(df, , "Meta_"))
assign(paste0("all_", w), lapply(all_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("p_", w),   lapply(p_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("wp_", w),  lapply(wp_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("wps_", w), lapply(wps_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("np_", w),  lapply(np_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS
run_analysis(get(paste0("all_", w)), "all", w)
run_analysis(get(paste0("p_", w)), "p", w)
run_analysis(get(paste0("wp_", w)), "wp", w)
run_analysis(get(paste0("wps_", w)), "wps", w)
run_analysis(get(paste0("np_", w)), "np", w)

## WINSORIZATION at 2% ####
w <- 2
# Renaming and Preprocessing Data
df <- general %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E2", "SE2", "N2", "studyID", "metaID")) #, "n_E_m")) 

# Splitting Data Based on Conditions
assign(paste0("all_", w), split_and_name(df, , "Meta_"))
assign(paste0("all_", w), lapply(all_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("p_", w),   lapply(p_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("wp_", w),  lapply(wp_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("wps_", w), lapply(wps_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("np_", w),  lapply(np_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS
run_analysis(get(paste0("all_", w)), "all", w)
run_analysis(get(paste0("p_", w)), "p", w)
run_analysis(get(paste0("wp_", w)), "wp", w)
run_analysis(get(paste0("wps_", w)), "wps", w)
run_analysis(get(paste0("np_", w)), "np", w)

## WINSORIZATION at 5% ####
w <- 5
# Renaming and Preprocessing Data
df <- general %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E5", "SE5", "N5", "studyID", "metaID")) #, "n_E_m")) 

# Splitting Data Based on Conditions
# Splitting Data Based on Conditions
assign(paste0("all_", w), split_and_name(df, , "Meta_"))
assign(paste0("all_", w), lapply(all_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("p_", w),   lapply(p_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("wp_", w),  lapply(wp_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("wps_", w), lapply(wps_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("np_", w),  lapply(np_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS
run_analysis(get(paste0("all_", w)), "all", w)
run_analysis(get(paste0("p_", w)), "p", w)
run_analysis(get(paste0("wp_", w)), "wp", w)
run_analysis(get(paste0("wps_", w)), "wps", w)
run_analysis(get(paste0("np_", w)), "np", w)



                                                # FIXED EFFECT ####
gen <- general %>%
  select("E1", "SE1", "N1", "E2", "SE2", "N2", "E5", "SE5", "N5","studyID", "metaID", "studyPublishD", "np" )
  #select("bs", "sebs", "Ns", "studyid", "meta_id", "studyPublishD", "np")


# Unique meta_ids
numbers <- gen$metaID %>% unique()

object<-c("E","SE","F","Hausman","CV_of_Chi2", "Obs")
# Function to Split and Name Data
split_and_name <- function(data, condition, prefix) {
  split_data <- split(data[condition, ], data$meta_id[condition])
  split_data <- lapply(split_data, function(df) {
    df$meta_id <- NULL
    return(df)
  })
  names(split_data) <- paste0(prefix, names(split_data))
  return(split_data)
}

# FUNCTION TO HANDLE THE REPETITIVE TASKS 
 run_analysis_fe <- function(list_of_meta, y, w) {
   source("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/GitHub/maive/maive-PETPEESE/maivefunction.R")
   # OPTIONS: 
   method <- 3 # method: PET:1, PEESE:2, PET-PEESE:3, EK:4 (default 3)
   weight <- 0 #no weight: 0; weights: 1, adjusted weights: 2 
   instrument <- 1 # instrumenting 
   studylevel <-1  # none: 0, fixed effects: 1, cluster: 2
   # default options are method=3; weight=0; instrument=1; studylevel=0 
   
   MAIVEresults <- data.frame(object)
   for (i in numbers){
     dat <- list_of_meta[[paste0("Meta_", i)]]
     MAIVE <- tryCatch(
       {
         maive(dat=dat, method=method, weight=weight, instrument=instrument, studylevel=studylevel)
       },
       error = function(e) {
         return(NULL)
       }
     )
     # If MAIVE is not NULL (meaning no error), extract results and save them
     if (!is.null(MAIVE)) {
       MAIVEresults[[paste0("Meta_", i)]] <- c(MAIVE$beta, MAIVE$SE, MAIVE$`F-test`, MAIVE$Hausman, MAIVE$Chi2, nrow(dat))
     } else {
       MAIVEresults[[paste0("Meta_", i)]] <- c(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
     }
     print(i) # To monitor the progress
   }
   # If you want to view the results after the loop
   print(MAIVEresults)
   

  results <- as.data.frame(t(MAIVEresults))
  colnames(results) <- results[1, ]
  results <- results[-1, ] %>%
    mutate(across(everything(), as.numeric))

  # Add row names as a new column in the data frame
  results_with_row_names <- results %>%
    tibble::add_column(Row_Names = row.names(results), .before = 1)

  # Save the modified data frame to an Excel file
  write_xlsx(results_with_row_names, paste0("MAIVE_", y, "FE_" , w, ".xlsx"))
 }

                                                          ## WINSORIZATION at 1% ####
w <- 1
# Renaming and Preprocessing Data
df <- gen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E1", "SE1", "N1", "studyID", "metaID")) #, "n_E_m")) 
# Splitting Data Based on Conditions
assign(paste0("fe_all_", w), split_and_name(df, , "Meta_"))
assign(paste0("fe_all_", w), lapply(fe_all_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("fe_p_", w),   lapply(fe_p_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("fe_wp_", w),  lapply(fe_wp_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("fe_wps_", w), lapply(fe_wps_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("fe_np_", w),  lapply(fe_np_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS
run_analysis_fe(get(paste0("fe_all_", w)), "all", w)
run_analysis_fe(get(paste0("fe_p_", w)), "p", w)
run_analysis_fe(get(paste0("fe_wp_", w)), "wp", w)
run_analysis_fe(get(paste0("fe_wps_", w)), "wps", w)
run_analysis_fe(get(paste0("fe_np_", w)), "np", w)

                                                    ## WINSORIZATION at 2% ####
w <- 2
# Renaming and Preprocessing Data
df <- gen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E2", "SE2", "N2", "studyID", "metaID")) #, "n_E_m")) 

# Splitting Data Based on Conditions
assign(paste0("fe_all_", w), split_and_name(df, , "Meta_"))
assign(paste0("fe_all_", w), lapply(fe_all_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("fe_p_", w),   lapply(fe_p_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("fe_wp_", w),  lapply(fe_wp_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("fe_wps_", w), lapply(fe_wps_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("fe_np_", w),  lapply(fe_np_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS
run_analysis_fe(get(paste0("fe_all_", w)), "all", w)
run_analysis_fe(get(paste0("fe_p_", w)), "p", w)
run_analysis_fe(get(paste0("fe_wp_", w)), "wp", w)
run_analysis_fe(get(paste0("fe_wps_", w)), "wps", w)
run_analysis_fe(get(paste0("fe_np_", w)), "np", w)

                                                      ## WINSORIZATION at 5% ####
w <- 5
# Renaming and Preprocessing Data
df <- gen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E5", "SE5", "N5", "studyID", "metaID")) #, "n_E_m")) 

# Splitting Data Based on Conditions
assign(paste0("fe_all_", w), split_and_name(df, , "Meta_"))
assign(paste0("fe_all_", w), lapply(fe_all_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_p_", w),   split_and_name(df, df$studyPublishD == 1, "Meta_"))
assign(paste0("fe_p_", w),   lapply(fe_p_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wp_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0, "Meta_"))
assign(paste0("fe_wp_", w),  lapply(fe_wp_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_wps_", w), split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 0, "Meta_"))
assign(paste0("fe_wps_", w), lapply(fe_wps_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("fe_np_", w),  split_and_name(df%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), df$studyPublishD == 0 & df$np == 1, "Meta_"))
assign(paste0("fe_np_", w),  lapply(fe_np_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS 
run_analysis_fe(get(paste0("fe_all_", w)), "all", w)
run_analysis_fe(get(paste0("fe_p_", w)), "p", w)
run_analysis_fe(get(paste0("fe_wp_", w)), "wp", w)
run_analysis_fe(get(paste0("fe_wps_", w)), "wps", w)
run_analysis_fe(get(paste0("fe_np_", w)), "np", w)
















                                                            # BETWEEN EFFECTS ####
gen <- general %>%
  select("E1", "SE1", "N1", "E2", "SE2", "N2", "E5", "SE5", "N5","studyID", "metaID", "studyPublishD", "np" )

avgen <- gen %>%
  group_by(metaID, studyID) %>%
  summarise(
    E1  = mean(E1, na.rm = TRUE),
    SE1 = mean(SE1, na.rm = TRUE),
    N1  = mean(N1, na.rm = TRUE),
    E2  = mean(E2, na.rm = TRUE),
    SE2 = mean(SE2, na.rm = TRUE),
    N2  = mean(N2, na.rm = TRUE),
    E5  = mean(E5, na.rm = TRUE),
    SE5 = mean(SE5, na.rm = TRUE),
    N5  = mean(N5, na.rm = TRUE),
    
    studyPublishD = mean(studyPublishD, na.rm = TRUE),
    np  = mean(np, na.rm = TRUE)
  )


# Unique meta_ids
numbers <- avgen$metaID %>% unique()
object<-c("E","SE","F","Hausman","CV_of_Chi2", "Obs")
# Function to Split and Name Data
split_and_name <- function(data, condition, prefix) {
  split_data <- split(data[condition, ], data$meta_id[condition])
  split_data <- lapply(split_data, function(df) {
    df$meta_id <- NULL
    return(df)
  })
  names(split_data) <- paste0(prefix, names(split_data))
  return(split_data)
}

# FUNCTION TO HANDLE THE REPETITIVE TASKS 
run_analysis_be <- function(list_of_meta, y, w) {
  source("/Users/ninobuliskeria/Library/CloudStorage/GoogleDrive-40603931@fsv.cuni.cz/My Drive/BIAS/GitHub/maive/maive-PETPEESE/maivefunction.R")
  # OPTIONS: 
  method <- 3 # method: PET:1, PEESE:2, PET-PEESE:3, EK:4 (default 3)
  weight <- 0 #no weight: 0; weights: 1, adjusted weights: 2 
  instrument <- 1 # instrumenting 
  studylevel <- 0  # none: 0, fixed effects: 1, cluster: 2
  # default options are method=3; weight=0; instrument=1; studylevel=0 
  
  MAIVEresults <- data.frame(object)
  for (i in numbers){
    dat <- list_of_meta[[paste0("Meta_", i)]]
    MAIVE <- tryCatch(
      {
        maive(dat=dat, method=method, weight=weight, instrument=instrument, studylevel=studylevel)
      },
      error = function(e) {
        return(NULL)
      }
    )
    # If MAIVE is not NULL (meaning no error), extract results and save them
    if (!is.null(MAIVE)) {
      MAIVEresults[[paste0("Meta_", i)]] <- c(MAIVE$beta, MAIVE$SE, MAIVE$`F-test`, MAIVE$Hausman, MAIVE$Chi2, nrow(dat))
    } else {
      MAIVEresults[[paste0("Meta_", i)]] <- c(NULL, NULL, NULL, NULL, NULL, NULL, NULL)
    }
    print(i) # To monitor the progress
  }
  # If you want to view the results after the loop
  print(MAIVEresults)
  
  results <- as.data.frame(t(MAIVEresults))
  colnames(results) <- results[1, ]
  results <- results[-1, ] %>%
    mutate(across(everything(), as.numeric))
  
  # Add row names as a new column in the data frame
  results_with_row_names <- results %>%
    tibble::add_column(Row_Names = row.names(results), .before = 1)
  
  # Save the modified data frame to an Excel file
  write_xlsx(results_with_row_names, paste0("MAIVE_", y, "BE_" , w, ".xlsx"))
}


## WINSORIZATION at 1% ####
w <- 1
# Renaming and Preprocessing Data
avdf <- avgen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E1", "SE1", "N1", "studyID", "metaID")) #, "n_E_m")) 


# Splitting Data Based on Conditions
assign(paste0("be_all_", w), split_and_name(avdf, , "Meta_"))
assign(paste0("be_all_", w), lapply(be_all_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_p_", w),   split_and_name(avdf, avdf$studyPublishD == 1, "Meta_"))
assign(paste0("be_p_", w),   lapply(be_p_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wp_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0, "Meta_"))
assign(paste0("be_wp_", w),  lapply(be_wp_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wps_", w), split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 0, "Meta_"))
assign(paste0("be_wps_", w), lapply(be_wps_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_np_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 1, "Meta_"))
assign(paste0("be_np_", w),  lapply(be_np_1, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS 
run_analysis_be(get(paste0("be_all_", w)), "all", w)
run_analysis_be(get(paste0("be_p_", w)), "p", w)
run_analysis_be(get(paste0("be_wp_", w)), "wp", w)
run_analysis_be(get(paste0("be_wps_", w)), "wps", w)
run_analysis_be(get(paste0("be_np_", w)), "np", w)



## WINSORIZATION at 2% ####
w <- 2
# Renaming and Preprocessing Data
avdf <- avgen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E2", "SE2", "N2", "studyID", "metaID")) #, "n_E_m")) 


# Splitting Data Based on Conditions
assign(paste0("be_all_", w), split_and_name(avdf, , "Meta_"))
assign(paste0("be_all_", w), lapply(be_all_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_p_", w),   split_and_name(avdf, avdf$studyPublishD == 1, "Meta_"))
assign(paste0("be_p_", w),   lapply(be_p_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wp_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0, "Meta_"))
assign(paste0("be_wp_", w),  lapply(be_wp_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wps_", w), split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 0, "Meta_"))
assign(paste0("be_wps_", w), lapply(be_wps_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_np_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 1, "Meta_"))
assign(paste0("be_np_", w),  lapply(be_np_2, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS 
run_analysis_be(get(paste0("be_all_", w)), "all", w)
run_analysis_be(get(paste0("be_p_", w)), "p", w)
run_analysis_be(get(paste0("be_wp_", w)), "wp", w)
run_analysis_be(get(paste0("be_wps_", w)), "wps", w)
run_analysis_be(get(paste0("be_np_", w)), "np", w)



## WINSORIZATION at 5% ####
w <- 5
# Renaming and Preprocessing Data
avdf <- avgen %>%
  rename_with(~c("bs", "sebs", "Ns", "studyid", "meta_id"), #, "ncoefm"), 
              .cols = c("E5", "SE5", "N5", "studyID", "metaID")) #, "n_E_m")) 

# Splitting Data Based on Conditions
assign(paste0("be_all_", w), split_and_name(avdf, , "Meta_"))
assign(paste0("be_all_", w), lapply(be_all_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_p_", w),   split_and_name(avdf, avdf$studyPublishD == 1, "Meta_"))
assign(paste0("be_p_", w),   lapply(be_p_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wp_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0, "Meta_"))
assign(paste0("be_wp_", w),  lapply(be_wp_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_wps_", w), split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 0, "Meta_"))
assign(paste0("be_wps_", w), lapply(be_wps_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

assign(paste0("be_np_", w),  split_and_name(avdf%>% select("bs", "sebs", "Ns", "studyid", "meta_id"), avdf$studyPublishD == 0 & avdf$np == 1, "Meta_"))
assign(paste0("be_np_", w),  lapply(be_np_5, function(x) x %>% select("bs", "sebs", "Ns", "studyid")))

# RUNNING THE FUNCTION FOR DIFFERENT SUBSETS 
run_analysis_be(get(paste0("be_all_", w)), "all", w)
run_analysis_be(get(paste0("be_p_", w)), "p", w)
run_analysis_be(get(paste0("be_wp_", w)), "wp", w)
run_analysis_be(get(paste0("be_wps_", w)), "wps", w)
run_analysis_be(get(paste0("be_np_", w)), "np", w)

























 
