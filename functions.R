#---------------- About ------------------
# Project: Recall error in the Patient-Oriented Eczema Measure (POEM) 

# This script contains common functions used throughout the project.

#-----------------------------------------

library(tidyverse)

#' Load in daily symptom data.
#'
#' @param file File path to dataset.
#' @param remove_NA Bool indicating whether to remove observations with NA or not.
#'
#' @return A dataframe containing the daily symptom scores of the clinical dataset.
#' @export
#'
#' @examples
load_daily <- function(file = "/homedir/data.RData", remove_NA = TRUE) {
  
  stopifnot(is_scalar_logical(remove_NA))
  
  load(file = file)
  
  if (remove_NA) {
    daily_df <- drop_na(daily_df)
  } else { # Drop rows with more than 1 missing of the 12 symptoms.
    twelve_symptoms <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep", "burning", "painful", "redness", "scaly", "stinging")
    
    # Count the number of missing symptoms for each observation.
    num_NA_12 <- rowSums(is.na(daily_df %>% select(twelve_symptoms))) 
    
    # Drop rows with more than 1 missing symptom.
    i_drop <- as.numeric(rownames(daily_df[num_NA_12 > 1, ])) 
    daily_df <- daily_df[-i_drop, ]
  }
  
  return(daily_df)
  
}


#' Load in observed POEM (r-POEM) data.
#'
#' @param file File path to dataset.
#' @param remove_NA Bool indicating whether to remove observations with NA or not.
#'
#' @return A dataframe containing the POEM scores of the clinical dataset.
#' @export
#'
#' @examples
load_POEM <- function(file = "/homedir/data.RData") {
  
  load(file = file)
  
  drop_cols <- c("Difficult to fall asleep after waking up", "How difficult to fall asleep", "How many nights wake up at least once")
  
  POEM_df <- POEM_df %>%
    select(!all_of(drop_cols)) %>%
    drop_na()
  
  return(POEM_df)
  
}


#' Helper function to add a column to front of dataframe.
#'
#' @param df Dataframe we want to add column to.
#' @param row_daily Data for column we want to add.
#' @param missing_day Name of column.
#'
#' @return df Dataframe with added column.
#' @export
#'
#' @examples
add_column_to_df <- function(df, row_daily, missing_day){
  
  q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")
  
  missing_day <- paste0("day_", missing_day)
  impute <- as.numeric(row_daily %>% select(subset = q_labels))
  
  df <- add_column(df, "{missing_day}" := impute, .before = 1) #library(tibble)
  
  return(df)
}


#' Helper function that converts, for a symptom, the number of present days to the POEM item score (0-4).
#'
#' @param num_days Numeric or vector of number of days.
#'
#' @return Numeric or vector of corresponding POEM item scores.
#' @export
#'
#' @examples
num_days_to_POEM_item <- function(num_days){
  
  dict <- c(0, 1, 1, 2, 2, 3, 3, 4)
  
  POEM_item_scores <- unlist(lapply(num_days, function(n_days){dict[n_days + 1]}))
  
  return(POEM_item_scores)
}


#' Helper function that generates a dataframe of the daily scores observed in the week corresponding to the POEM observation.
#'
#' @param ID_rows_daily Rows of daily_df that belong to the patient we are currently interested in.
#' @param day_POEM Day in study of POEM observation.
#'
#' @return Dataframe with rows representing the 7 measured items and columns the scores of the days observed in the week previous to observed POEM.
#' Content are the scores for each item.
#' The dataframe will have 7 columns if all 7 days in the week are measured. If only 1 day measured then only 1 column.
#' e.g. of a dataframe returned:
#'           day_24 day_25 day_26 day_27 day_28 day_29 day_30
#'bleeding      0      0      0      0      1      0      0     
#'cracked       0      1      0      0      0      0      0    
#'dry           5      4      5      7      4      4      5     
#'flaky         3      3      3      3      2      3      4    
#'itchy         4      4      5      5      6      5      4    
#'oozing        0      0      0      0      1      0      0     
#'sleep         3      4      3      2      2      2      2     
#' @export
#'
#' @examples
generate_week_of_daily <- function(ID_rows_daily, day_POEM){
  
  q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")
  df <- data.frame(item = q_labels)
  
  # For each day in the previous week, add its daily observations of the 7 items to the dataframe if measured.
  # Count day of observation as one of the 7 days in the week.
  for (day_before in (day_POEM-6):day_POEM){ 
    row_daily <- ID_rows_daily %>% filter(Analysis.Relative.Day == day_before)
    
    if(dim(row_daily)[1] != 0){ # If we have a daily observation for that day.
      df <- cbind(df, as.numeric(row_daily %>% select(subset = all_of(q_labels))))
      colnames(df)[ncol(df)] <- paste0("day_", day_before)
    }
  }
  rownames(df) <- q_labels
  df <- df %>% select(subset = -c("item"))
  
  return(df)
  
}


#' Calculate POEM score from the scores of the past 7 days' symptoms.
#'
#' @param row_input Row from POEM or daily dataframe for corresponding week we want to calculate POEM.
#' @param type_row "POEM" or "daily" to indicate from which dataset the row comes from.
#' @param threshold Threshold to convert daily symptom score (0-10) to binary presence of symptoms for POEM. 
#' Can be NA if scores have already been thresholded. Days with scores >= threshold are considered as having symptoms,
#' otherwise they are considered absent (score < threshold).
#'
#' @return A dataframe with the same patient information as the row_input. Plus for each symptom and the total score (total 8 rows), 
#' the observed_POEM (r-POEM) score, the calculated_POEM (d-POEM) score, and the daily scores for the corresponding week.
#'
#' @examples
calculate_POEM <- function(row_input, type_row, threshold = 2, daily_df) {
  
  stopifnot(isTRUE(type_row %in% c("POEM", "daily")))
  
  q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")
  
  day_obs <- row_input$Analysis.Relative.Day
  
  ID_rows_daily <- daily_df %>% filter(ID == row_input$ID)
  
  df <- generate_week_of_daily(ID_rows_daily = ID_rows_daily, day_POEM = day_obs)
  
  # Calculate POEM if we have all 7 days.
  if (ncol(df) == 7) { 
    if(is.na(threshold)){ # Threshold is NA if the scores are already dichotomized (thresholded).
      # For each item, count the number of days in the week with score present.
      df$count <- rowSums(df) 
    } else{
      # For each item, count the number of days in the week with score >= threshold.
      df$count <- rowSums(df >= threshold) 
    }
    
    # Convert number of days to POEM item score 0-4.
    df <- df %>% mutate(POEM_score_calculated = num_days_to_POEM_item(count))
    
    # POEM score from 0-28 (sum of item scores from 0-4).
    total_POEM <- sum(df$POEM_score_calculated)
    
    # Total daily scores (sum the symptom scores for each of the 7 days in the week).
    total_daily_scores <- as.numeric(colSums(df[, 1:7]))
    total_daily_df <- t(as.data.frame(total_daily_scores))
    colnames(total_daily_df) <- paste0("day", seq(1,7), "_observed")
    
    # Daily scores.
    daily7_df <- as.data.frame(df[, 1:7])
    colnames(daily7_df) <- paste0("day", seq(1,7), "_observed")
    
    # Add all information into a dataframe.
    # First, add information for the total POEM score.
    rows_7days <- row_input %>%
      select(ID, Patient, Age, Country, Race, Sex, Treatment,  Analysis.Relative.Day) %>%
      rename(Day = Analysis.Relative.Day) %>%
      mutate(Item = "total",
             calculated_POEM = total_POEM) %>%
      bind_cols(total_daily_df) 
    
    # Second, add information for the POEM item scores.
    d <- bind_cols(rows_7days %>% select(ID, Patient, Age, Country, Race, Sex, Treatment,  Day), data.frame(Item = q_labels, calculated_POEM = df$POEM_score_calculated), daily7_df)
    rows_7days <- rows_7days %>% bind_rows(d)
    
    rownames(rows_7days) <- seq(1, nrow(rows_7days))
    
    if(type_row == "POEM"){ # Add information about observed POEM scores.
      rows_7days <- rows_7days %>% add_column(observed_POEM = as.numeric(row_input %>% select(subset = c("Total POEM Score", q_labels))), .after = "Day")
    }
    
    out <- rows_7days
  } else {
    out <- NULL
  }
  
  return(out)
}


#' Generate observed (r-POEM) and calculated (d-POEM) POEM dataframe.
#'
#' @param threshold Threshold to convert daily symptom score (0-10) to binary presence of symptoms for POEM.
#' Days with scores >= threshold are considered as having symptoms otherwise they are considered absent (score < threshold).
#' #' Can be NA if scores have already been thresholded.
#' @param daily_df dataframe of daily observations.
#' @param POEM_df dataframe of POEM observations.  
#'
#' @return A dataframe with the observed (r-POEM) and calculated (d-POEM) POEM scores for each item and their total score 
#' for all comparisons for which we are able to calculate a corresponding POEM score. As well as the daily observed scores 
#' for the week measured by POEM.
#' @export
#'
#' @examples
generate_observed_calculated <- function(threshold = 2, daily_df, POEM_df) {
  
  df <- which(POEM_df$Analysis.Relative.Day != 1) %>%
    lapply(function(i) {
      calculate_POEM(POEM_df[i, ], type_row = "POEM", threshold = threshold, daily_df = daily_df)
    }) %>%
    bind_rows()
  
  cols_numeric = c("Age", "Day", "observed_POEM", "calculated_POEM", "day1_observed", "day2_observed", "day3_observed", "day4_observed", "day5_observed", "day6_observed", "day7_observed")
  df[, cols_numeric] = apply(df[, cols_numeric], 2, function(x) as.numeric(x))
  
  return(df)
}


#' Generate observed (r-POEM) and calculated (d-POEM) POEM dataframe with a range of thresholds.
#' @param threshold_low Lowest value of threshold to generate calculated with.
#' @param threshold_high Highest value of threshold to generate calculated with.
#' @param daily_df dataframe of daily observations.
#' @param POEM_df dataframe of POEM observations.
#'
#' @return Dataframe that compiles outputs of generate_observed_calculated() for a range of thresholds.
#' A dataframe with the observed (r-POEM) and calculated (d-POEM) POEM scores for each item and their total score 
#' for all comparisons for which we are able to calculate a corresponding POEM score. As well as the daily observed scores 
#' for the week measured by POEM. For a range of thresholds.
#' @export
#'
#' @examples 
generate_varying_thresholds <- function(threshold_low = 1, threshold_high = 10, daily_df, POEM_df) {
  
  stopifnot(between(threshold_low, 1, 10),
            between(threshold_high, 1, 10),
            threshold_low <= threshold_high)
  
  # Use the generate_observed_calculated() function to generate calculated POEM for a range of thresholds.
  seq(threshold_low:threshold_high) %>%
    lapply(function(i) {
      generate_observed_calculated(threshold = i, daily_df, POEM_df) %>%
        bind_cols(data.frame(threshold = i))
    }) %>%
    bind_rows()
}


#' Visualize frequency of scores in the different datasets.
#'
#' @param df Dataframe including items to plot ("observed_calculated_df", "daily_df" or "POEM_(full_)df").
#' @param to_plot Name of severity measurement want to plot ("observed_POEM", "calculated_POEM", "daily" or "POEM")
#'
#' @return Bar plot with frequencies of the scores for each severity item in the dataframe.
#' @export
#'
#' @examples
visualize_freq <- function(df, to_plot){
  
  stopifnot(to_plot %in% c("observed_POEM", "calculated_POEM", "daily", "POEM"))
  
  q_labels <- q_labels <- c("bleeding", "cracked", "dry", "flaky", "itchy", "oozing", "sleep")
  
  # Two ways of plotting since the full dataframes (daily_df and POEM_df) have a different
  # format than the observed_calculated_df.
  # Option 1: df is observed_calculated_df
  if(to_plot %in% c("observed_POEM", "calculated_POEM")){
    pl <- lapply(c("total", q_labels),
                 function(item) {
                   title <- ifelse(item == "total", yes = "POEM", no = item)
                   
                   p <- ggplot(data = df %>% filter(Item == item), aes(x = .data[[to_plot]])) + 
                     geom_bar(color = "black", fill = "grey") + 
                     scale_y_continuous(expand = c(0, 0), sec.axis = sec_axis(~ . / (nrow(df %>% filter(Item == item)))))+
                     labs(x = "Score", y = "Num. Observations", title = title) +
                     theme_bw(base_size = 15) +
                     theme(panel.grid.minor.x = element_blank())
                   
                   if(item == "total"){
                     p = p + scale_x_continuous(breaks = seq(0, 28, 4))
                   }else{
                     p = p + scale_x_continuous(breaks = seq(0, 4))
                   }
                 }
    )
    
  # Option 2: df is daily_df or POEM_df
  }else{ # to_plot == "daily" or "POEM"
    pl <- lapply(q_labels,
                 function(item) {
                   p <- ggplot(data = df, aes(x = .data[[item]])) + 
                     geom_bar(color = "black", fill = "grey") + 
                     scale_y_continuous(expand = c(0, 0), sec.axis = sec_axis(~ . / (length(df[,item]))))+
                     labs(x = "Score", y = "Num. Observations", title = item) +
                     theme_bw(base_size = 15) +
                     theme(panel.grid.minor.x = element_blank())
                   
                   if(to_plot == "daily"){
                     p = p + scale_x_continuous(breaks = seq(0, 10))
                   }else{ # to_plot == "POEM"
                     p = p + scale_x_continuous(breaks = seq(0, 4))
                   }
                 })
  }
  return(plot_grid(plotlist = pl, ncol = 3))
}


#' Obtain (r-POEM) observed POEM dataset with added column of the number of days with daily scores recorded for the corresponding week.
#'
#' @param daily_df dataframe of daily observations.
#' @param POEM_df dataframe of POEM observations.
#'
#' @return POEM_df with an added column of the number of days in the corresponding POEM week we have daily symptom scores recorded.
#' All rows with Analysis.Relative.Day < 7 are removed, as we cannot calculate a POEM score with 7 days of data.
#' @export
#'
#' @examples
generate_numdays_present <- function(daily_df, POEM_df) {
  
  # We cannot calculate a POEM score from the daily symptom scores without 7 days of data, so discard rows with study day < 7.
  df <- which(POEM_df$Analysis.Relative.Day > 6) %>% 
    lapply(function(i) {
      row_input <- POEM_df[i, ]
      day_obs <- row_input$Analysis.Relative.Day
      
      ID_rows_daily <- daily_df[daily_df$ID == row_input$ID, ]
      df <- generate_week_of_daily(ID_rows_daily = ID_rows_daily, day_POEM = day_obs)
      
      # Number of days we have daily score for the week of POEM
      return(row_input %>% mutate(num_days_daily_scores_present = ncol(df)))
    }) %>%
    bind_rows()
  
  return(df)
}


#' Obtain dataframe listing missing days (no daily symptom scores recorded) that we need to calculate POEM (d-POEM).
#'
#' @param daily_df dataframe of daily observations.
#' @param POEM_df dataframe of POEM observations.
#'
#' @return Dataframe with 3 columns: Patient, ID, and Analysis.Relative.Day (missing days).
#' @export
#'
#' @examples
missing_daily_days <- function(daily_df, POEM_df) {
  
  # We cannot calculate a POEM score from the daily symptom scores without 7 days of data, so discard rows with study day < 7.
  df <- which(POEM_df$Analysis.Relative.Day > 6) %>% 
    lapply(function(i) {
      row_input <- POEM_df[i, ]
      day_obs <- row_input$Analysis.Relative.Day
      
      ID_rows_daily <- daily_df %>% filter(ID == row_input$ID)
      df <- generate_week_of_daily(ID_rows_daily = ID_rows_daily, day_POEM = day_obs)
      
      days_recorded <- as.numeric(regmatches(colnames(df), gregexpr("[[:digit:]]+", colnames(df))))
      missing_days <- setdiff(seq(day_obs - 6,day_obs), days_recorded)
      
      if(length(missing_days) > 0){
        return(cbind(ID_rows_daily[1, ] %>% select(c("Patient", "ID")), data.frame(Analysis.Relative.Day = missing_days), row.names = NULL))
      }
    }) %>%
    bind_rows()
  
  return(df)
}


#' Predict probability of daily presence or absence from fitted model.
#' For a threshold of 2 (predicted daily score >= 2 means presence of symptoms, absence otherwise).
#'
#' @param fit Stan fit object (obtained from fitting training set on Markov chain model).
#' @param predict_df_MC Dataframe with the scores we want to predict. Has columns with Patient, Time, y0 (previous score), and dt (time since y0).
#'
#' @return predict_df_MC with added column of the probability of abs/pres (0 to 1) and dichotomized probability (0 or 1).
#' @export
#'
#' @examples
predict_daily_presabs <- function(fit, predict_df_MC){
  # Extract the transition matrix P for all possible dt (delayed time).
  # a has dimension [number samples in chains from MCMC, max dt in training set, max score + 1 (= 11), max score + 1 (= 11)]
  a <- rstan::extract(fit, pars = "P")[[1]]
  
  # Take the mean across the samples, so we have a transition probability for each possible state transition for each dt.
  b <- apply(a, c(2, 3, 4), mean)
  
  # We can only predict scores with dt <= max dt in the training dataset (which is already very large, e.g. 70). 
  predict_df_MC <- predict_df_MC %>%
    filter(dt <= dim(b)[1])
  
  # Obtain probability of presence.
  # For threshold 2, to get the probability of symptom being present, sum probabilities of state being 3 or higher
  # States in Markov Chain model range from 1 to 11, instead of 0 to 10, so state of 1 or 2 means absence and 3 or higher means presence for threshold of 2.
  predict_df_MC[["prob_symp"]] <- vapply(1:nrow(predict_df_MC),
                                         function(i) {
                                           with(predict_df_MC,
                                                sum(b[dt[i], y0[i], -c(1, 2)])) 
                                         }, numeric(1))
  
  # Dichotomize probability into presence (1) or absence (0). 
  # 0 if probability < 0.5, 1 otherwise.
  predict_df_MC <- predict_df_MC %>%
    mutate(impute = as.numeric(prob_symp > 0.5))
  
  return(predict_df_MC)
}
