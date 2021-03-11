rm(list = ls())
#----------------------------------------------------------------------------
# Process the COVID19 Data (last day is 2020-12-02)
#The Data include county-wise fatality, case count,population , mobility
library(tidyverse)
library(dplyr)
library(reshape2)
#----------------------------------------------
#1) Tidy: Fatality count data
Tx_COVID_Fat <- read_csv("data/Texas COVID-19 Fatality Count Data by County.csv")

dat_f = Tx_COVID_Fat %>% pivot_longer(
  cols = colnames(Tx_COVID_Fat)[-1], 
  names_to = "date", 
  values_to = "fat_total") %>%
  dplyr::mutate(date = as.Date(date, "%m/%d/%Y"))%>% 
  rename(County = "County Name") %>%
  group_by(County) %>% 
  mutate(fat_new = c(NA, diff(fat_total))) #new fatalities

#2) Tidy: Case count data
Tx_COVID_Case <- read_csv("data/Texas COVID-19 Fatality Count Data by County.csv")

dat_c = Tx_COVID_Case %>% pivot_longer(
  cols = colnames(Tx_COVID_Case)[-1], 
  names_to = "date", 
  values_to = "case_total") %>%
  dplyr::mutate(date = as.Date(date, "%m/%d/%Y"))%>% 
  rename(County = "County Name") %>%
  group_by(County) %>% 
  mutate(case_new = c(NA, diff(case_total))) #new cases

#3) Tidy: Population data
Tx_Pop <- read_delim("data/Texas County Population.csv", ";", escape_double = FALSE, col_types = cols(X1 = col_skip(),Population = col_double()), trim_ws = TRUE)

dat_p = Tx_Pop[complete.cases(Tx_Pop$Population), ]%>%
  dplyr::rename(County = `County Name`)%>%
  distinct() %>%
  filter(County != "Total")%>% 
  mutate(Population = (Population/max(Population,na.rm = TRUE))) #scale population to be < 1

#4) Tidy: Mobility data
Tx_Mob <- read_csv("data/GMR_TX.csv", col_types = cols(X1 = col_skip(),date = col_date(format = "%Y-%m-%d")))

#covariates: percentage -> decimals
dat_m = Tx_Mob %>% 
  rename(County = county)%>%
  dplyr::select(-c("parks"))%>% 
  group_by(County) %>%
  summarise(retail = mean(as.numeric(retail)/100,na.rm = TRUE), 
            grocery = mean(as.numeric(grocery)/100,na.rm = TRUE), 
            transit = mean(as.numeric(transit)/100,na.rm = TRUE), 
            workplaces = mean(as.numeric(workplaces)/100,na.rm = TRUE), 
            residential = mean(as.numeric(residential)/100,na.rm = TRUE)) 

#-----------------------------------------------------
#4):Combine dat_f,dat_c, dat_p,dat_m to get X,Y
XY_df = dat_f %>% 
  inner_join(dat_c, on = c("County","date"))%>%
  inner_join(dat_p, on = c("County"))%>%
  inner_join(dat_m, on = c("County"))%>%
  filter(County != "Total")%>% 
  filter(case_total >= 10)%>% #remove county with <10 total cases
  mutate(date_ind = seq(length(fat_new))) %>% 
  group_by(County) %>% mutate(n_obs = n()) %>%filter(n_obs >= 10)%>% ungroup() #remove county with <10 observations
  

Y = XY_df%>%
  dplyr::select("County","date_ind","fat_new") %>%
  dcast(County ~ date_ind) %>% 
  column_to_rownames(var = "County")%>%
  as.matrix()

X = XY_df %>% 
  dplyr::select(c("Population","retail","grocery","transit","workplaces","residential")) %>%
  distinct() %>%
  as.matrix()

#------------------------------------------------------
# Save the X, Y, Ytot, tau, and the survey weights:
save(X, Y, file = paste('data/pa_covid19.RData', sep=''))

