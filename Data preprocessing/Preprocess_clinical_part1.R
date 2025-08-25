# Loading required packages.
library("factoextra")
library("FactoMineR")
library("R.matlab")
library("matlabr")
library("dplyr")

# Getting file directory.
getwd()
# Setting the directory.
setwd("C:/Users/Data")

# Uploading the data set. #GTEX = healthy #TCGA = with cancer.

BreastCancer <-
  read.csv(file = "Clinical_Original.csv",
           stringsAsFactors = FALSE,
           header = TRUE)

# Delete unnecessary data. 
BreastCancer <-subset(BreastCancer, type=="BRCA")

# Resetting the row names so they start from 1.
rownames(BreastCancer) <- NULL

colnames(BreastCancer)

BreastCancer = as.data.frame(BreastCancer)

# Making sure data is deleted and checking new observations.
dim(BreastCancer)

# Finding the data type for all attributes
summary(BreastCancer)

sapply(BreastCancer, function(x) sum(is.na(x)))

BreastCancer <-
  BreastCancer[-c(1047, 1059),]


sapply(BreastCancer,class)

# Resetting the row names so they start from 1.
rownames(BreastCancer) <- NULL


################## ATTRIBUTES TO DELETE ##################

table(BreastCancer$clinical_stage, exclude=NULL) 
which(BreastCancer$clinical_stage != '[Not Applicable]', arr.ind=T)


table(BreastCancer$histological_grade, exclude=NULL) 
which(BreastCancer$histological_grade != '[Not Available]', arr.ind=T)


BreastCancer$birth_days_to <- as.double(BreastCancer$birth_days_to)
summary(BreastCancer$birth_days_to)
which(is.na(BreastCancer$birth_days_to) | NA)  


BreastCancer$last_contact_days_to <- as.double(BreastCancer$last_contact_days_to)
summary(BreastCancer$last_contact_days_to)
which(is.na(BreastCancer$last_contact_days_to) | NA)      


table(BreastCancer$death_days_to, exclude=NULL) 


table(BreastCancer$cause_of_death, exclude=NULL) 


table(BreastCancer$new_tumor_event_type, exclude=NULL) 
which(BreastCancer$new_tumor_event_type != '#N/A')


table(BreastCancer$new_tumor_event_site, exclude=NULL) 
which(!BreastCancer$new_tumor_event_site == '#N/A' 
      | BreastCancer$new_tumor_event_site == 'Other, specify')



table(BreastCancer$new_tumor_event_site_other, exclude=NULL) 
which(BreastCancer$new_tumor_event_site_other != '#N/A' )


BreastCancer$new_tumor_event_dx_days_to <- as.double(BreastCancer$new_tumor_event_dx_days_to)
table(BreastCancer$new_tumor_event_dx_days_to, exclude=NULL) 
summary(BreastCancer$new_tumor_event_dx_days_to)
which(!is.na(BreastCancer$new_tumor_event_dx_days_to) | NA)      


table(BreastCancer$treatment_outcome_first_course, exclude=NULL) 


table(BreastCancer$residual_tumor, exclude=NULL) 


table(BreastCancer$DSS, exclude=NULL) 
which(BreastCancer$DSS == '#N/A')


which(BreastCancer$DSS.time == '#N/A')


table(BreastCancer$DFI, exclude=NULL) 
which(BreastCancer$DFI == '#N/A')


table(BreastCancer$DFI.time, exclude=NULL) 
which(BreastCancer$DFI.time == '#N/A')


table(BreastCancer$PFI, exclude=NULL) 
which(BreastCancer$PFI.time == '#N/A')


table(BreastCancer$Redaction, exclude=NULL) 
which(BreastCancer$Redaction == 'Redacted')



BreastCancer[ , c(7,9,12,15,16,17,18,19,20,21,22,24,27,28,29,30,31,32,33)] <- list(NULL)




rownames(BreastCancer) <- NULL


################################### Pre processing ############################

table(BreastCancer$type, exclude = NULL)


BreastCancer$age_at_initial_pathologic_diagnosis[BreastCancer$age_at_initial_pathologic_diagnosis == '#N/A' ] <- '58'
BreastCancer$age_at_initial_pathologic_diagnosis <- as.double(BreastCancer$age_at_initial_pathologic_diagnosis)
summary(BreastCancer$age_at_initial_pathologic_diagnosis)
table(BreastCancer$age_at_initial_pathologic_diagnosis, exclude=NULL) 
which(is.na(BreastCancer$age_at_initial_pathologic_diagnosis) | NA)  


# Finding the rows that do not contain the gender as female.
table(BreastCancer$gender, exclude=NULL) 
which(!grepl('FEMALE',BreastCancer$gender))
BreastCancer <- BreastCancer[-c(18, 199, 286, 381, 392, 424, 530, 557, 768, 818, 976, 990),]
rownames(BreastCancer) <- NULL



table(BreastCancer$race, exclude=NULL) 
which(!grepl('^[A-Z]',BreastCancer$race))
BreastCancer$race[BreastCancer$race == "[Not Evaluated]" ] <- "[Not Available]"


table(BreastCancer$ajcc_pathologic_tumor_stage, exclude=NULL) 
which(!grepl('^[A-Z]',BreastCancer$ajcc_pathologic_tumor_stage))
BreastCancer$ajcc_pathologic_tumor_stage[BreastCancer$ajcc_pathologic_tumor_stage == "[Discrepancy]" ] <- "[Not Available]"


table(BreastCancer$histological_type, exclude=NULL) 
which(BreastCancer$histological_type == '[Not Available]' | 
        BreastCancer$histological_type == 'Other, specify')
BreastCancer$histological_type[BreastCancer$histological_type == "[Not Available]" |
                                 BreastCancer$histological_type == "Mixed Histology (please specify)"] <- "Other, specify"


BreastCancer$initial_pathologic_dx_year[BreastCancer$initial_pathologic_dx_year == '#N/A'] <- 2009
BreastCancer$initial_pathologic_dx_year <- as.double(BreastCancer$initial_pathologic_dx_year)
table(BreastCancer$initial_pathologic_dx_year, exclude=NULL) 
summary(BreastCancer$initial_pathologic_dx_year)
which(is.na(BreastCancer$initial_pathologic_dx_year) | NA)      
                                        

table(BreastCancer$menopause_status, exclude=NULL) 
which(BreastCancer$menopause_status == '[Not Available]' | 
        BreastCancer$menopause_status == '[Unknown]' |
        BreastCancer$menopause_status == '[Not Evaluated ]' |  
        BreastCancer$menopause_status == 'Other, specify')
BreastCancer$menopause_status[BreastCancer$menopause_status == "[Not Evaluated]" |
                                 BreastCancer$menopause_status == "[Unknown]"] <- "[Not Available]"


table(BreastCancer$vital_status, exclude=NULL) 


table(BreastCancer$tumor_status, exclude=NULL) 
which(!grepl('^[A-Z]',BreastCancer$tumor_status))
BreastCancer$tumor_status[BreastCancer$tumor_status == "#N/A"] <- "[Not Available]"
BreastCancer$tumor_status[BreastCancer$tumor_status == "[Not Available]"] <- "WITH TUMOR"


table(BreastCancer$margin_status, exclude=NULL) 
which(BreastCancer$margin_status == '[Not Available]'
      | BreastCancer$margin_status == '[Unknown]', arr.ind=T)
BreastCancer$margin_status[BreastCancer$margin_status == "[Unknown]"] <- "[Not Available]"


table(BreastCancer$OS, exclude=NULL) 


BreastCancer$OS.time <- as.numeric(BreastCancer$OS.time)
table(BreastCancer$OS.time, exclude=NULL) 

rownames(BreastCancer) <- NULL




write.csv(BreastCancer,'Survival(Pre-Processed).csv', row.names=FALSE)
