# Function for the kidney data preprocessing with the same approach 
# of Thei et al., 2024

# @return the cleaned dataset, also save the corresponding file in the fold of .\Data
Data_proprocessing1 <- function(){

  
  ## locate the file on your computer
  NOTR_DGF.df <- read_excel("../Data/NOTR2.xlsx", sheet = 'NOTR data alle centra_20230509_', guess_max=50000)
  
  
  
  # create patient survival variables
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(
      PT_surv_time = difftime(Recipientdateofdeath, Dateoftransplant, units = "days"),
      PT_Death = if_else(is.na(Recipientdateofdeath), 0,1),
      PT_Earlydeath = if_else(PT_surv_time < 10 & PT_Death == 1, 1, 0),
      Year_of_death = format(Recipientdateofdeath, format = "%Y")) 
  
  
  # rename some variables 
  NOTR_DGF.df<- NOTR_DGF.df%>%
    dplyr::rename("Number_of_transplantations" = "TransplantKidneyordernrorgan",
                  "Recipient_weight" = "Initialweight",
                  "Recipient_height" = "Initialheight",
                  "Recipient_EBV" = "InitialEBVIGG",
                  "Recipient_CMV" = "InitialCMV",
                  'Recipient_Vascular_event' = "InitialVasculaireventindicator",
                  "Recipient_Diabetes" = "InitialDiabeteseverindicator",
                  "Recipient_CVA" = "InitialCVAindicator",
                  "Recipient_Cardiac_event" = "InitialCardiaceventindicator",
                  "Warm_ischaemic_period_1" = "Warmischaemicperiod1",
                  "Warm_ischaemic_period_2" = "InitialWarmischaemicperiod2",
                  "Cold_ischaemic_period" = "InitialColdischaemicperiod",
    )
  
  
  # remove some variables
  NOTR_DGF.df<- NOTR_DGF.df%>%
    select(-Transplantmethod, -Lowestcreatininemgdl, -Highestcreatininemgdl)
  
  # replace "Unknown" values with NA
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(
      Recipient_Vascular_event = na_if(Recipient_Vascular_event, "Unknown"),
      Recipient_Diabetes = na_if(Recipient_Diabetes, "Unknown"),
      Recipient_CVA = na_if(Recipient_CVA, "Unknown"),
      Recipient_Cardiac_event = na_if(Recipient_Cardiac_event, "Unknown"),
      Warm_ischaemic_period_1 = na_if(Warm_ischaemic_period_1, 0),
      Warm_ischaemic_period_2 = na_if(Warm_ischaemic_period_2, 0), 
      Lowestcreatininemumoll = case_when(Lowestcreatinine¦Ìmoll <1001 ~ Lowestcreatinine¦Ìmoll,
                                         Lowestcreatinine¦Ìmoll >1000 ~NA),
      Cold_ischaemic_period = ifelse(Cold_ischaemic_period <0, 0,Cold_ischaemic_period )
    )%>%
    mutate(across(where(is.numeric), ~ifelse(. == -1, NA, .)))
  
  # Create some calculated variables
  NOTR_DGF.df$DonorBMI <- NOTR_DGF.df$Donorweight / ((NOTR_DGF.df$Donorheight/100)^2)
  NOTR_DGF.df$Recipient_BMI <- NOTR_DGF.df$Recipient_weight/((NOTR_DGF.df$Recipient_height/100)^2)
  NOTR_DGF.df$Recipient_BMI <- ifelse(NOTR_DGF.df$Recipient_BMI > 16 & NOTR_DGF.df$Recipient_BMI < 60, NOTR_DGF.df$Recipient_BMI, NA)
  NOTR_DGF.df$DonorBMI <- ifelse(NOTR_DGF.df$DonorBMI > 16 & NOTR_DGF.df$DonorBMI < 60, NOTR_DGF.df$DonorBMI, NA)
  NOTR_DGF.df$Year_Of_Tx = NOTR_DGF.df$Dateoftransplant
  NOTR_DGF.df$Year_Of_Tx <- as.numeric(format(NOTR_DGF.df$Year_Of_Tx, format = "%Y"))
  
  # Add an ID column which will help us to identify recipients who have received a transplant before.
  NOTR_DGF.df<-
    NOTR_DGF.df%>% group_by(fictiefrecipientnummer) %>% mutate(ID = paste0(fictiefrecipientnummer, ".", 1:n())) 
  
  # Create some categorical variables for the cox models later
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(Recipient_age.cat = cut(Recipientage, breaks = c(18, 25, 45, 65, Inf))) %>%
    mutate(Donor_age.cat = cut(Donorage, breaks = c(18, 25, 45, 65, Inf))) %>%
    mutate(Donor_BMI.cat = cut(DonorBMI, breaks = c(0, 35,  Inf)))
  
  ## here we will categorize the causes of death that have been listed in the data set
  # Create a vector of conditions for each cause of death group
  circulatory_conditions <- c("Circulational: Acute Myocard Infarct",
                              "Circulational: Cardiac Arrest",
                              "Circulational: Not Otherwise Specified",
                              "Circulational: Recidief Myocard Infarct")
  icva_conditions <- c("CVA: Cerebral Ischemia",
                       "CVA: Cerebro Vascular Accident Not Otherwise Specified")
  
  bcva_conditions <- c(
    "CVA: Intra Cerebral Bleeding",
    "SAB: Sub Arachnoidal Bleeding")
  
  other_conditions <- c("Medical complication: Diagnostic treatment",
                        "Medical complication: Surgical / Medical treatment",
                        "Non-accident: Not Otherwise Specified",
                        "Not Otherwise Specified",
                        "rej Epilepticus",
                        "Suicide: Not Otherwise Specified",
                        "Meningitis: Bacterial",
                        "Brain Tumor: Astrocytoma grade 1 or 2", 
                        "Brain Tumor: Astrocytoma grade 3",
                        "Brain Tumor: Benign",
                        "Brain Tumor: Malignant")
  
  respiratory_conditions <- c("Respirational: Epiglotitis / Laryngitis",
                              "Respirational: Not Otherwise Specified",
                              "Respirational: rej Asthmaticus",
                              "Suicide: Respiratory",
                              "Trauma: Suffocation",
                              "Trauma: Drowning")
  trauma_conditions <- c("Suicide: Head injury",
                         "Suicide: Jump",
                         "Trauma: Capitis",
                         "Trauma: Falling",
                         "Trauma: Mechanical",
                         "Trauma: Travel accident on land",
                         "Trauma: Not Otherwise Specified",
                         "EDH: Epi Dural Hematoma",
                         "SDH: Sub Dural Hematoma")
  
  # NOTR_sex.df$DonorDeathGroup.dbd <-
  #   case_when(
  #     is.element(NOTR_sex.df$d_cause_of_death_ET, circulatory_conditions) ~ "Cardiac",
  #     is.element(NOTR_sex.df$d_cause_of_death_ET , icva_conditions) |
  #       is.element(NOTR_sex.df$d_cause_of_death_ET , bcva_conditions) ~ "CVA",
  #     is.element(NOTR_sex.df$d_cause_of_death_ET , respiratory_conditions) |
  #       is.element(NOTR_sex.df$d_cause_of_death_ET , other_conditions) ~ "Other",
  #     is.element(NOTR_sex.df$d_cause_of_death_ET , trauma_conditions) ~ "Trauma",
  #   ) # HAVING NO IDEA ABOUT WHAT THIS FILE IS
  
  
  ## Creating graft survival time
  # Setting end date for living recipient with a functioning graft to the last day of follow-up (day of data-extraction)
  caseDay <- ymd("2023-05-09")
  
  # The formatting of dates is a little tricky from excel to R so that is what I'm doing here. It does change for every version of excel so check if this also aplies to your version
  NOTR_DGF.df<-
    NOTR_DGF.df%>% mutate(InitialGraftfaildate1 = ymd(InitialGraftfaildate)) %>%
    mutate(Dateoftransplant1 = ymd(Dateoftransplant))%>%
    mutate(last_FU_date = case_when(!is.na(Recipientdateofdeath) ~ ymd(Recipientdateofdeath),
                                    is.na(Dateseen.END) & is.na(Recipientdateofdeath) ~ caseDay,
                                    is.na(Recipientdateofdeath) & !is.na(Dateseen.END) ~ ymd(Dateseen.END)
                                    #                                 
                                    # mutate(last_FU_date = case_when(!is.na(r_date_of_death) ~ ymd(r_date_of_death),
                                    #                                 is.na(date_last_seen) & is.na(r_date_of_death)  ~ caseDay,
                                    #                                 is.na(r_date_of_death) & !is.na(date_last_seen) ~ ymd(date_last_seen)
    ))  #Dateseen.END == date_last_seen?
  
  
  # Now calculate the times that we're interested in
  NOTR_DGF.df<-
    NOTR_DGF.df%>%
    mutate(Follow_up_time = as.numeric(difftime(last_FU_date, Dateoftransplant1,  units = 'weeks')) /
             52.25) %>%
    mutate(DGF_duration = difftime(ymd(DGFLastdialysisdate), ymd(Dateoftransplant), units = 'days')) %>%
    mutate(GraftLoss = if_else(is.na(InitialGraftfaildate)&
                                 is.na(InitialGraftfailcause), 0, 1)) 
  
  
  
  #### HLA matching (package transplantr)
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(HLA_mismatch = hla_mm_level(a = MismatchA, b = MismatchB, dr = MismatchDR))
  
  
  # Here we create the PNF variable. Rationale/definition in supplements 
  NOTR_DGF.df<-
    NOTR_DGF.df %>% 
    mutate(
      PNF = case_when(
        Delayedgraftfunction == 'Never' |
          InitialGraftfailcause == 'Permanent Non-Function' |
          InitialGraftfailcause == 'Non-Viable Kidney' |
          DGF_duration > 90|
          InitialGraftfaildate1 < 7 + Dateoftransplant1 ~ 1,
        TRUE ~ 0  
      )
    ) %>%
    mutate(
      PNF = case_when( InitialGraftfailcause ==
                         "Patient died with functioning transplant" ~ 0 , TRUE ~ PNF)) 
  
  #Make it a factor
  NOTR_DGF.df$PNF <- as.factor(NOTR_DGF.df$PNF)
  
  
  
  # Impute a graft fail date for the recipients with PNF and a missing graft fail date
  NOTR_DGF.df$InitialGraftfaildate <-
    case_when(
      NOTR_DGF.df$PNF == 1 &
        is.na(NOTR_DGF.df$InitialGraftfaildate1) ~
        NOTR_DGF.df$Dateoftransplant1 + 90,
      TRUE ~ NOTR_DGF.df$InitialGraftfaildate1
    )
  
  
  # Create the delayed graft function variable
  NOTR_DGF.df$dDGF <- as.factor(
    case_when(
      NOTR_DGF.df$Delayedgraftfunction == "Direct" |
        NOTR_DGF.df$Delayedgraftfunction == "Never" |
        NOTR_DGF.df$PNF == 1 ~ 0, 
      NOTR_DGF.df$Delayedgraftfunction == "Delayed" |
        !is.na(NOTR_DGF.df$DGFLastdialysisdate) ~ 1,
      TRUE ~ NA
    )
  )
  
  
  # In case of missingness in dDGF, check if creatinine is available to judge if subjected to dialysis 
  # to avoid unnecessary work we only do this for the cases we'll include
  
  ## INCLUSION & EXCLUSION CRITERIA
  # check.creat <- subset(NOTR_DGF.df,Typeofdonor=='Deceased' 
  #                       & Recipientage >17
  #                       & Donorage >17
  #                       & Number_of_transplantations==1
  #                       & (Combinedtransplants=='LKi'
  #                          | Combinedtransplants=='RKi') 
  #                       & PNF== 0
  #                       # & is.na(Transatl) WHAT DOES THIS MEAN??
  #                       & PT_Earlydeath ==0
  #                       & !is.na(dDGF)
  #                       & Year_Of_Tx > 2013
  #                       & Year_Of_Tx < 2023)
  # 
  # check.creat2 <- gather(check.creat, day, creatinine, CreatinineD1:CreatinineD7, factor_key = TRUE)
  # check.creat2 <- subset(check.creat2, !is.na(creatinine))
  # check.creat2 <- check.creat2%>% arrange(ID)
  # check.creat3 <- check.creat2[1:500,]
  
  # ggplot(check.creat3,
  #        aes(
  #          x = day,
  #          y = creatinine,
  #          color = factor(ID),
  #          group = factor(ID)
  #        )
  # ) +
  #   geom_point() +
  #   geom_line() +
  #   facet_wrap( ~ factor(ID), nrow = 4)
  # 
  #saveRDS(check.creat2, "check.creat2.rds")
  
  
  ### The following ID's were either clearly subjected to dialysis or clearly not subjected to dialysis
  Missing.No.dDGF <- c(
    "394849.1",
    "401860.1",
    "544225.1",
    "545601.1",
    "560034.1",
    "567893.1",
    "572340.1",
    "573015.1",
    "573287.1",
    "578053.1",
    "578271.1",
    "583026.1"
  )
  
  Missing.dDGF <- c("541652.1",
                    "541993.1",
                    "561811.1")
  
  
  
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(dDGF = if_else(ID %in% Missing.No.dDGF, factor(0), dDGF)) %>%
    mutate(dDGF = if_else(ID %in%  Missing.dDGF, factor(1), dDGF))
  
  NOTR_DGF.df$Donorsex <- as.factor(NOTR_DGF.df$Donorsex)
  
  NOTR_DGF.df<- NOTR_DGF.df%>% mutate(dDGF = relevel(as.factor(dDGF), ref = "0"))%>% 
    mutate(Donorsex = relevel(Donorsex, ref = "Male"))
  
  ########## Create groups for graft loss 
  NOTR_DGF.df$Donor_type_deceased <- NOTR_DGF.df$TypecadavericDBDDCD 
  
  NOTR_DGF.df<- NOTR_DGF.df%>%
    mutate(Donor_type_deceased = case_when(Donor_type_deceased == 'Donation after brain death' ~ "DBD",
                                           Donor_type_deceased == 'Donation after circulatory death' ~ "DCD"))
  
  #### Create variables for competing risk analysis this is the death-censored graft loss definition we'll use for outcome analyses
  NOTR_DGF.df$time_tx_death <- round(difftime(NOTR_DGF.df$Recipientdateofdeath,NOTR_DGF.df$Dateoftransplant,units="days"))
  
  NOTR_DGF.df$time_tx_graftloss <- round(difftime(NOTR_DGF.df$InitialGraftfaildate,NOTR_DGF.df$Dateoftransplant,units="days"))
  
  NOTR_DGF.df$time_tx_date_max <- round(difftime(NOTR_DGF.df$last_FU_date,NOTR_DGF.df$Dateoftransplant,units="days"))
  
  NOTR_DGF.df$time_tx_date_max2 <- round(difftime(date("2023-05-09"),NOTR_DGF.df$Dateoftransplant,units="days"))
  
  
  NOTR_DGF.df$diff_tx_death <- NOTR_DGF.df$time_tx_graftloss - NOTR_DGF.df$time_tx_death 
  
  NOTR_DGF.df$status <- NA
  NOTR_DGF.df$status <- ifelse(!is.na(NOTR_DGF.df$InitialGraftfaildate) | NOTR_DGF.df$PNF==1, "graftloss", NA )
  NOTR_DGF.df$status <- ifelse(!is.na(NOTR_DGF.df$diff_tx_death) & NOTR_DGF.df$diff_tx_death >=0 , "death", NOTR_DGF.df$status)
  NOTR_DGF.df$status <- ifelse(NOTR_DGF.df$InitialGraftfailcause == "Patient died with functioning transplant","death", NOTR_DGF.df$status)
  NOTR_DGF.df$status <- ifelse(is.na(NOTR_DGF.df$status), "alive", NOTR_DGF.df$status)
  
  NOTR_DGF.df$time <- ifelse(NOTR_DGF.df$status=="graftloss", NOTR_DGF.df$time_tx_graftloss, NA)
  NOTR_DGF.df$time <- ifelse(NOTR_DGF.df$status=="death", NOTR_DGF.df$time_tx_death, NOTR_DGF.df$time)
  NOTR_DGF.df$time <- ifelse(NOTR_DGF.df$status=="alive", NOTR_DGF.df$time_tx_date_max2,  NOTR_DGF.df$time)
  
  
  ## create factor variable for other analyses
  NOTR_DGF.df$status <- as.factor(NOTR_DGF.df$status)
  NOTR_DGF.df$status2 <- ifelse(NOTR_DGF.df$status=="death", 2, NA)
  NOTR_DGF.df$status2 <- ifelse(NOTR_DGF.df$status=="graftloss", 1, NOTR_DGF.df$status2)
  NOTR_DGF.df$status2 <- ifelse(NOTR_DGF.df$status=="alive", 0, NOTR_DGF.df$status2)
  # NOTR_DGF.df$status2 <- as.factor(NOTR_DGF.df$status2)
  
  
  # Check if all recipients with graft loss date are labeled correctly
  table(is.na(NOTR_DGF.df$InitialGraftfaildate), NOTR_DGF.df$status, exclude =NULL)
  
  
  # only plasma creatinine is registered in the data set so we'll calculate eGFR with the BMI at inclusion
  columns_eGFR <- names(NOTR_DGF.df)[grepl("Creatinine", names(NOTR_DGF.df), ignore.case = FALSE)]
  
  NOTR_DGF.df.egfr<- NOTR_DGF.df%>%
    select(ID, Recipientage, Recipientsex, columns_eGFR) %>%
    gather(day, creatinine, columns_eGFR, factor_key = TRUE)
  
  # Using ethnicity is not recomended in Europe
  NOTR_DGF.df.egfr$ethnicity <- 0
  
  
  # we use the 2009 formula in the ckd-epi package
  NOTR_DGF.df.egfr <- NOTR_DGF.df.egfr %>%
    mutate(eGFR = ckd_epi(creat = creatinine, age = Recipientage, sex = Recipientsex, ethnicity = ethnicity))
  
  
  # check if it looks ok
  summary(NOTR_DGF.df.egfr$eGFR)
  
  # now for the variable/column names 
  NOTR_DGF.df.egfr <- select(NOTR_DGF.df.egfr, ID, eGFR,day)
  NOTR_DGF.df.egfr$day <-
    factor(
      NOTR_DGF.df.egfr$day,
      levels = c( 'CreatinineD1', 'CreatinineD2',
                  'CreatinineD3', 'CreatinineD4',
                  'CreatinineD5', 'CreatinineD6',
                  'CreatinineD7',
                  "Creatinine.M3" ,
                  "Creatinine.Y01",
                  "Creatinine.Y02",
                  "Creatinine.Y03" ,
                  "Creatinine.Y04",
                  "Creatinine.Y05",
                  "Creatinine.Y06",
                  "Creatinine.Y07" ,
                  "Creatinine.Y08",
                  "Creatinine.Y09" ,
                  "Creatinine.Y10",
                  "Creatinine.Y11",
                  "Creatinine.Y12",
                  "Creatinine.Y13" ,
                  "Creatinine.Y14",
                  "Creatinine.Y15",
                  "Creatinine.Y16" ,
                  "Creatinine.Y17",
                  "Creatinine.Y18" ,
                  "Creatinine.Y19" ,
                  "Creatinine.Y20",
                  "Creatinine.Y21",
                  "Creatinine.Y22",
                  "Creatinine.Y23",
                  "Creatinine.END"
      ),
      labels = c( 'eGFRD1', 'eGFRD2',
                  'eGFRD3', 'eGFRD4',
                  'eGFRD5', 'eGFRD6',
                  'eGFRD7',
                  'eGFR.M3',
                  'eGFR.Y01',
                  'eGFR.Y02',
                  'eGFR.Y03',
                  'eGFR.Y04',
                  'eGFR.Y05',
                  'eGFR.Y06',
                  'eGFR.Y07',
                  'eGFR.Y08',
                  'eGFR.Y09',
                  'eGFR.Y10',
                  'eGFR.Y11',
                  'eGFR.Y12',
                  'eGFR.Y13',
                  'eGFR.Y14',
                  'eGFR.Y15',
                  'eGFR.Y16',
                  'eGFR.Y17',
                  'eGFR.Y18',
                  'eGFR.Y19',
                  'eGFR.Y20',
                  'eGFR.Y21',
                  'eGFR.Y22',
                  'eGFR.Y23',
                  'eGFR.END'
      )
    )
  
  # back to wide format
  NOTR_DGF.df.egfr.short <- spread(NOTR_DGF.df.egfr[,!names(NOTR_DGF.df.egfr) %in% c("fictiefrecipientnummer")], "day", "eGFR")
  
  # add the eGFR's to the dataset
  NOTR_DGF.df<- merge(NOTR_DGF.df, NOTR_DGF.df.egfr.short, by.x = "ID", by.y = "ID")
  
  
  ## Done with the data cleaning. on to the creatinine imputation
  saveRDS(NOTR_DGF.df, "./Data/NOTR_DGF.Rds")
  
  return (NOTR_DGF.df)
}


# Function for the kidney data preprocessing
# @param dataset: a dataframe, which was cleaned 
#                 with the same approach of  Thei et al., 2024
# @return: the reorganized dataset, in a dataframe class
Real_world_data_reorganization <- function(dataset){
  ## INCLUSION & EXCLUSION CRITERIA
  df_using <- subset(dataset,Typeofdonor=='Deceased' 
                     & Recipientage >17
                     & Donorage >17
                     & Number_of_transplantations==1
                     & (Combinedtransplants=='LKi'
                        | Combinedtransplants=='RKi') 
                     & Year_Of_Tx > 2013
                     & Year_Of_Tx < 2023)
  
  
  # select the necessary variables
  variableName_df = read.csv(file = "../Data/codeboek2_.csv")
  
  variableName = variableName_df$Name_NOTR 
  variableName = if_else(variableName == "Initialweight", "Recipient_weight", variableName) # use the renamed name
  variableName = if_else(variableName == "Initialheight", "Recipient_height", variableName) # use the renamed name
  variableName = if_else(variableName == "Warmischaemicperiod1", "Warm_ischaemic_period_1", variableName) 
  variableName = if_else(variableName == "InitialColdischaemicperiod", "Cold_ischaemic_period", variableName) 
  variableName = if_else(variableName == "LowestcreatinineÂµmoll", "Lowestcreatininemumoll", variableName)
  variableName_using = variableName[nzchar(variableName)]
  variableName_using = variableName_using[!variableName_using %in% c("InitialPrimaryDisease", 
                                                                     "InitialPrimaryDiseaserenine",
                                                                     "DonorHLAstring")]
  
  variableName_using <- c(variableName_using, c("time", "status", "status2", "HLA_mismatch")) # including survival related variables
  df_using_ <- df_using |> select(variableName_using)
  df_using_$Lowestcreatininemumoll <- df_using$Lowestcreatininemumoll
  df_using <- df_using_
  head(df_using)
  str(df_using)
  
  print(variableName_using)
  
  
  
  
  ### re-label/factorization variables which have 2 or three unique values
  df_using = df_using |>
    mutate(Recipientsex = factor(Recipientsex, levels = c("Male","Female"), labels  = c("Male","Female")),
           Donorsex = factor(Donorsex, levels = c("Female", "Male"), labels  = c("Female", "Male")),
           Smoking = factor(Smoking, levels = c("Yes", "No"), labels  = c("Yes", "No")),
           TypecadavericDBDDCD = factor(TypecadavericDBDDCD, levels = c("Donation after circulatory death", "Donation after brain death"), labels  = c("DCD", "DBD"))
    )
  
  
  
  # # group the InitialPrimaryDiseaseET
  InitialPrimaryDiseaseET_group = read.csv("../Data/tx_all_primary disease.csv")
  print(head(InitialPrimaryDiseaseET_group))
  
  completed_df_ <- dataset


  Num_overlap <- 0
  for (InitialPrimaryDiseaseET in sort(unique(as.character(completed_df_$InitialPrimaryDiseaseET)))){
    if (InitialPrimaryDiseaseET %in% colnames(InitialPrimaryDiseaseET_group)){
      Num_overlap = Num_overlap + 1
      print(InitialPrimaryDiseaseET)
    }
    else{
      
    }
  }
  Num_overlap
  library(dplyr)
  completed_df_ <- completed_df_ |> mutate(InitialPrimaryDiseaseET_regroup = case_when(
    InitialPrimaryDiseaseET == "Amyloidosis" ~ "Other multisystem diseases",
    InitialPrimaryDiseaseET == "Chronic renal failure - etiology uncertain"~ "Chronic renal failure; aetiology uncertain", #"Chronic.renal.failure...etiology.uncertain",
    InitialPrimaryDiseaseET == "Congenital renal dysplasia with/without urinary tract malformation" ~ "Other congenital and hereditary kidney diseases", #"Congenital.renal.dysplasia.with.without.urinary.tract.malformation",
    InitialPrimaryDiseaseET == "Congenital renal hypoplasia - Type unspecified" ~ "Other congenital and hereditary kidney diseases", #"Congenital.renal.hypoplasia...Type.unspecified",
    InitialPrimaryDiseaseET == "Cortical or Tubular necrosis" ~  "Others", #"Cortical.or.Tubular.necrosis",
    
    
    InitialPrimaryDiseaseET == "Cryoglobulinemic glomerulonephritis" ~ "Other multisystem diseases", #"Cryoglobulinemic.glomerulonephritis",
    InitialPrimaryDiseaseET == "Cystic kidney disease - Other specified type" ~ "Cystic kidney diseases", #"Cystic.kidney.disease...Other.specified.type",
    InitialPrimaryDiseaseET == "Cystic kidney disease - Type unspecified" ~ "Cystic kidney diseases", #"Cystic.kidney.disease...Type.unspecified", 
    InitialPrimaryDiseaseET == "Cystinosis" ~ "Other congenital and hereditary kidney diseases", #"Cystinosis",
    InitialPrimaryDiseaseET == "Diabetes Type I" ~ "Diabetes Mellitus", #"Diabetes.Type.I",
    
    
    InitialPrimaryDiseaseET == "Diabetes Type II" ~ "Diabetes Mellitus", # "Diabetes.Type.II", 
    InitialPrimaryDiseaseET == "Focal segmental glomerulosclerosis with nephrotic syndrome in adults" ~ "Glomerulonephritis", #"Focal.segmental.glomerulosclerosis.with.nephrotic.syndrome.in.adults", 
    InitialPrimaryDiseaseET == "Focal segmental glomerulosclerosis with nephrotic syndrome in children" ~ "Glomerulonephritis", #"Focal.segmental.glomerulosclerosis.with.nephrotic.syndrome.in.children", 
    InitialPrimaryDiseaseET == "Glomerulonephritis - histologically examined" ~ "Glomerulonephritis", 
    #"Glomerulonephritis...histologically.examined", 
    InitialPrimaryDiseaseET == "Glomerulonephritis - histologically not examined" ~ "Glomerulonephritis", #"Glomerulonephritis...histologically.not.examined",
    
    
    InitialPrimaryDiseaseET == "Goodpasture's syndrome" ~ "Other multisystem diseases", #"Goodpasture.s.syndrome", 
    InitialPrimaryDiseaseET == "Gout nephropathy (urate)" ~ "Others", #"Gout.nephropathy..urate.", 
    InitialPrimaryDiseaseET == "Hemolytic Uremic Syndrome including Moschcowitz syndrome" ~ "Other multisystem diseases", #"Hemolytic.Uremic.Syndrome.including.Moschcowitz.syndrome", 
    InitialPrimaryDiseaseET == "Henoch-Schonlein Purpura" ~ "Other multisystem diseases", #"Henoch.Schonlein.Purpura", 
    InitialPrimaryDiseaseET == "Hereditary nephritis with nerve deafness (Alport's syndrome)" ~ "Other congenital and hereditary kidney diseases", #"Hereditary.nephritis.with.nerve.deafness..Alport.s.syndrome.", 
    
    
    InitialPrimaryDiseaseET == "Hereditary Nephropathy - Other" ~ "Other congenital and hereditary kidney diseases", #"Hereditary.Nephropathy...Other", 
    InitialPrimaryDiseaseET == "Hereditary/Familial nephropathy - Type unspecified" ~ "Other congenital and hereditary kidney diseases", #"Hereditary.Familial.nephropathy...Type.unspecified",
    InitialPrimaryDiseaseET ==  "IgA nephropathy (proven by immunofluorescence)" ~ "Glomerulonephritis", #"IgA.nephropathy..proven.by.immunofluorescence." , 
    InitialPrimaryDiseaseET == "Iscemic renal disease / cholesterol embolism" ~ "Renal vascular disease, excluding vasculitis", #"Iscemic.renal.disease...cholesterol.embolism", 
    InitialPrimaryDiseaseET ==  "Kidney tumor" ~ "Others", # "Kidney.tumor", 
    
    
    InitialPrimaryDiseaseET == "Lupus Erythematosus" ~ "Other multisystem diseases", #"Lupus.Erythematosus", 
    InitialPrimaryDiseaseET == "Membrano-proliferative glomerulonephritis - Type I" ~ "Glomerulonephritis", #"Membrano.proliferative.glomerulonephritis...Type.I", 
    InitialPrimaryDiseaseET == "Membranous nephropathy" ~ "Glomerulonephritis", #"Membranous.nephropathy", 
    InitialPrimaryDiseaseET == "Multisystem disease - Other specified type" ~ "Other multisystem diseases", #"Multisystem.disease...Other.specified.type", 
    InitialPrimaryDiseaseET == "Myelomatosis - Light chain deposit disease" ~ "Other multisystem diseases", #"Myelomatosis...Light.chain.deposit.disease", #"Myelomatosis...Light.chain.deposit.disease", 
    
    
    InitialPrimaryDiseaseET == "Nephropathy caused by other specific drug" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", # "Nephropathy.caused.by.other.specific.drug",
    InitialPrimaryDiseaseET == "Nephropathy due to Cisplatinum" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Nephropathy.due.to.Cisplatinum", 
    InitialPrimaryDiseaseET == "Nephropathy due to Cyclosporin A" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Nephropathy.due.to.Cyclosporin.A", 
    InitialPrimaryDiseaseET == "Oligomeganephronic hypoplasia" ~ "Other congenital and hereditary kidney diseases", #"Oligomeganephronic.hypoplasia", 
    InitialPrimaryDiseaseET == "Other identified renal disorders - Specify" ~ "Others", #"Other.identified.renal.disorders...Specify", #*****
    
    
    InitialPrimaryDiseaseET == "Polycystic Kidneys - Adult type (dominant)" ~ "Cystic kidney diseases", # "Polycystic.Kidneys...Adult.type..dominant.", 
    InitialPrimaryDiseaseET == "Primary oxalosis" ~ "Other congenital and hereditary kidney diseases", # "Primary.oxalosis", 
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Acquired obstructive uropathy" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...Acquired.obstructive.uropathy", 
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Cause not specified" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...Cause.not.specified",  
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Congenital obstructive uropathy" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...Congenital.obstructive.uropathy",
    
    
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Neurogenic bladder" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...Neurogenic.bladder",
    
    
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Other cause" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", # "Pyelonephritis.Interstitial.nephritis...Other.cause",
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - Urolithiasis" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...Urolithiasis", 
    InitialPrimaryDiseaseET == "Pyelonephritis/Interstitial nephritis - VU reflux without obstruction" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", #"Pyelonephritis.Interstitial.nephritis...VU.reflux.without.obstruction", 
    InitialPrimaryDiseaseET == "Rapidly progressive Glomerulonephritis without systemic disease" ~ "Glomerulonephritis", #"Rapidly.progressive.Glomerulonephritis.without.systemic.disease",
    InitialPrimaryDiseaseET == "Renal vascular disease - Classified" ~ "Renal vascular disease, excluding vasculitis", # "Renal.vascular.disease...Classified", 
    
    
    InitialPrimaryDiseaseET == "Renal vascular disease - Type unspecified" ~  "Renal vascular disease, excluding vasculitis",#"Renal.vascular.disease...Type.unspecified", 
    InitialPrimaryDiseaseET == "Renal vascular disease due to hypertension" ~ "Renal vascular disease, excluding vasculitis", #"Renal.vascular.disease.due.to.hypertension", 
    InitialPrimaryDiseaseET == "Renal vascular disease due to malignant hypertension" ~ "Renal vascular disease, excluding vasculitis", #"Renal.vascular.disease.due.to.malignant.hypertension", 
    InitialPrimaryDiseaseET == "Renal vascular disease due to polyarteritis" ~ "Other multisystem diseases", # "Renal.vascular.disease.due.to.polyarteritis", 
    InitialPrimaryDiseaseET == "Systemic sclerosis (scleroderma)" ~ "Other multisystem diseases", #"Systemic.sclerosis..scleroderma.", 
    
    
    InitialPrimaryDiseaseET == "Traumatic or Surgical loss of kidney" ~ "Others", #"Traumatic.or.Surgical.loss.of.kidney", 
    InitialPrimaryDiseaseET == "Tuberculosis" ~ "Others", #"Tuberculosis",
    InitialPrimaryDiseaseET == "Tubulo-interstitial nephritis (not Pyelonephritis)" ~ "Interstitial nephritis, including pyelonephritis, drug induced nephropathy and urolithiasis", # "Tubulo.interstitial.nephritis..not.Pyelonephritis.", 
    InitialPrimaryDiseaseET == "Wegener's granulomatosis" ~ "Other multisystem diseases",#"Wegener.s.granulomatosis",
    InitialPrimaryDiseaseET == "Nephrocalcinosis and Hypercalcemic Nephropathy" ~ "Others"#"Nephrocalcinosis.and.Hypercalcemic.Nephropathy"
    
  ))
  
  
  df_using <- subset(completed_df_, select = -InitialPrimaryDiseaseET)
  df_using$InitialPrimaryDiseaseET_regroup = factor(df_using$InitialPrimaryDiseaseET_regroup)
  return(df_using)
}
