#Create table of percentage of hematopoietic-related populations relative
#to all cell types on each day of differentiation

library(gt)
library(tidyverse)
library(dplyr)

save_dir <- "~/save_dir/"

cds <- readRDS("~/sci_cds_BBI_preprocessed.RDS")
cds_hemato_only <- readRDS("~/sci_cds_BBI_preprocessed_hemato-only.RDS")

#Determine number of cells on each day of differentiation for Line 1
sci_d7_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 7"])
sci_d8_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 8"])
sci_d11_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 11"])
sci_d14_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 14"])
sci_d18_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 18"])
sci_d21_line1 <- ncol(cds[,colData(cds)$line_and_time.point == "45-iPS Day 21"])

#Determine number of cells on each day of differentiation for Line 2
sci_d7_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 7"])
sci_d8_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 8"])
sci_d11_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 11"])
sci_d14_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 14"])
sci_d18_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 18"])
sci_d21_line2 <- ncol(cds[,colData(cds)$line_and_time.point == "MSC-iPS Day 21"])

#Determine number of Hematopoietic-related cells on each day of differentiation for Line 1
hemato_d7_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 7"])
hemato_d8_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 8"])
hemato_d11_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 11"])
hemato_d14_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 14"])
hemato_d18_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 18"])
hemato_d21_line1 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "45-iPS Day 21"])

#Determine number of Hematopoietic-related cells on each day of differentiation for Line 2
hemato_d7_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 7"])
hemato_d8_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 8"])
hemato_d11_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 11"])
hemato_d14_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 14"])
hemato_d18_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 18"])
hemato_d21_line2 <- ncol(sci_hemato[,colData(sci_hemato)$line_and_time.point == "MSC-iPS Day 21"])

#Create vectors of calculated cell numbers separated by line and all or hematopoietic-related
line1_days_all <- c(sci_d7_line1, sci_d8_line1, sci_d11_line1, sci_d14_line1, 
                    sci_d18_line1, sci_d21_line1)
line2_days_all <- c(sci_d7_line2, sci_d8_line2, sci_d11_line2, sci_d14_line2, 
                    sci_d18_line2, sci_d21_line2)
line1_days_hemato <- c(hemato_d7_line1, hemato_d8_line1, hemato_d11_line1, 
                       hemato_d14_line1, hemato_d18_line1, hemato_d21_line1)
line2_days_hemato <- c(hemato_d7_line2, hemato_d8_line2, hemato_d11_line2, 
                       hemato_d14_line2, hemato_d18_line2, hemato_d21_line2)

#Calculate percentage of hematopoietic-related cells for each line and day of
#differentiation
line1_perc_hemato <- round((line1_days_hemato/line1_days_all)*100, digits = 2)
line2_perc_hemato <- round((line2_days_hemato/line2_days_all)*100, digits = 2)

#Create dataframes of the percentage data
line_1_df <- data.frame("Hematopoietic" = line1_days_hemato, 
                        "All" = line1_days_all, 
                        "Percent" = line1_perc_hemato)
line_1_df$day <- c("Day 7", "Day 8", "Day 11", "Day 14", "Day 18", "Day 21")
line_1_df <- line_1_df %>% select("day", "Hematopoietic", "All", "Percent")
colnames(line_1_df) <- c("Day of Differentiation",
                         "Hematopoietic-Related (HR)", 
                         "All Cells", 
                         "% HR")

line_2_df <- data.frame("Hematopoietic" = line2_days_hemato, 
                        "All" = line2_days_all, 
                        "Percent" = line2_perc_hemato)
line_2_df$day <- c("Day 7", "Day 8", "Day 11", "Day 14", "Day 18", "Day 21")
line_2_df <- line_2_df %>% select("day", "Hematopoietic", "All", "Percent")
colnames(line_2_df) <- c("Day of Differentiation",
                         "Hematopoietic-Related (HR)", 
                         "All Cells",
                         "% HR")

#Create tibbles of the dataframes for gt
line_1_tib <- as_tibble(line_1_df)
line_2_tib <- as_tibble(line_2_df)

#create a separate table for each line
line_1_tbl <- gt(line_1_tib, rowname_col = "Day of Differentiation") |>
  tab_header(title = "Day of Differentiation Cell Numbers (Line 1)") |> 
  gtsave(filename = paste0(save_dir, "line1_cell_numbers.png"))

line_2_tbl <- gt(line_2_tib, rowname_col = "Day of Differentiation") |>
  tab_header(title = "Day of Differentiation Cell Numbers (Line 2)") |> 
  gtsave(filename = paste0(save_dir, "line2_cell_numbers.png"))

#Create a single table containing the percentages for both lines 
days_line1 <- c("day_7_1", "day_8_1", "day_11_1", "day_14_1", "day_18_1", "day_21_1")
for (i in 1:6) {
  assign(days_line1[i], paste0(as.character(line1_days_hemato[i]), 
                               " / ", 
                               as.character(line1_days_all[i]), 
                               " (",
                               as.character(line1_perc_hemato[i]),
                               ")"))
}

days_line2 <- c("day_7_2", "day_8_2", "day_11_2", "day_14_2", "day_18_2", "day_21_2")
for (i in 1:6) {
  assign(days_line2[i], paste0(as.character(line2_days_hemato[i]), 
                               " / ", 
                               as.character(line2_days_all[i]), 
                               " (",
                               as.character(line2_perc_hemato[i]),
                               ")"))
}

cell_percents_df <- data.frame("Day of Differentiation" = c("Day 7", "Day 8", "Day 11", "Day 14", "Day 18", "Day 21"),
                               "Line 1" = c(day_7_1, day_8_1, day_11_1, day_14_1, day_18_1, day_21_1),
                               "Line 2" = c(day_7_2, day_8_2, day_11_2, day_14_2, day_18_2, day_21_2))

colnames(cell_percents_df) <- c("Day of Differentiation", "Line 1", "Line 2")
cell_percents_tib <- as_tibble(cell_percents_df)

tbl <- gt(cell_percents_tib, rowname_col = "Day of Differentiation") |>
  tab_stubhead(label = md("*Day of Differentiation*")) |>
  tab_header(title = md("**Hematopoietic-Related Population Percentages**"),
             subtitle = md("**Hematopoietic-Related Cells / All Cells (%)**")) |>
  cols_label("Line 1" = md("*Line 1*"),
             "Line 2" = md("*Line 2*")) |>
  tab_options(heading.padding.horizontal = px(40)) |>
  gtsave(filename = paste0(save_dir, "all_line_cell_numbers.png"))
