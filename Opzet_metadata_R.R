
#===================================================== Imports and sub selection ==================================================
set.seed(456)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(writexl)
library(tibble)


#File that holds all metadata on how and where swabs were taken
SwabMetaData <- read_excel("Masters BiBc/Research stage NIOZ/NIOZ365_Analyse/MetaDataJasper.xlsx")

#File that holds all metadata on how DNA extractions were done 
ExtrMetaData <- read_excel("Masters BiBc/Research stage NIOZ/Resultaten DNA-extracties/Extractie_metadata.xlsx")
ExtrMetaData$Extraction_batch <- as.factor(ExtrMetaData$Extraction_batch)
ExtrMetaData$Amount_of_material = as.numeric(ExtrMetaData$Amount_of_material)

#Subsetting the swab metadata to only those samples on which I did DNA extractions
SwabMetaData <- subset(SwabMetaData, SwabMetaData$RINGNR %in% ExtrMetaData$Ringnumber)
SwabMetaData <- subset(SwabMetaData, Condition != "RECAPTURE_OWN") #Removing an extra entry of Z101082 that was a recapture and did not have measurements other than weight

#Checking which samples I performed DNA extractions on do not have swab collection metadata (should only be my negative controls)
absent_ringnr_counter = 0
for (i in ExtrMetaData$Ringnumber) {
  if (!(i %in% SwabMetaData$RINGNR)) {
    absent_ringnr_counter = absent_ringnr_counter + 1
    print(i)
  } 
}
absent_ringnr_counter
#write.csv(SwabMetaData, "Masters BiBc/Research stage NIOZ/SwabMetaData.csv")

#Adding the catch date as a single column
SwabMetaData <- SwabMetaData %>% add_column("CATCH_DATE(YMD)" = NA)
for (row in 1:nrow(SwabMetaData)){
  year = SwabMetaData[row,"CATCH_YEAR"]
  month = SwabMetaData[row, "CATCH_MONTH"]
  day = SwabMetaData[row, "CATCH_DAY"]
  
  datum = paste(year[1,1],"-",month[1,1],"-",day[1,1], sep = "")
  #print(date)
  
  SwabMetaData[row, "CATCH_DATE(YMD)"] = datum
}

#write.csv(SwabMetaData, "Masters BiBc/Research stage NIOZ/SwabMetaData2.csv")

# ================================================= General data exploration =========================================================
#Checking if the amount of entries in both dataframes makes sense 
if (!((length(SwabMetaData$RINGNR) + absent_ringnr_counter) == length(ExtrMetaData$Ringnumber))){
  print("WARNING: Number of entries unequal!")
} else {
  print("Number of entries is equal")
}

#Obtaining counts of sex and age in the data
SwabMetaData %>% count(AGE)
SwabMetaData %>% count(SEX)
SwabMetaData %>% count(CATCH_LOCATION)
SwabMetaData %>% count(PLUM)
SwabMetaData %>% count(CATCH_YEAR)

# ================================================ Plotting swab collection metadata =================================================

#Sex per age group
ggplot(data=SwabMetaData, aes(x=AGE, fill=SEX)) +geom_bar()
chisq.test(SwabMetaData$SEX, SwabMetaData$AGE, correct = T)

#Sex vs mass
ggplot(data = SwabMetaData, aes(x=SEX, y=MASS)) +geom_boxplot()

#Catch location counts
ggplot(data=SwabMetaData, aes(y=CATCH_LOCATION, fill=as.factor(CATCH_MONTH))) + geom_bar()

#Bill length per age group
ggplot(data=SwabMetaData, aes(x=as.factor(AGE), y=BILL)) +geom_boxplot()
t.test(BILL~AGE, data = SwabMetaData)

#Preparedness for migration per age group
ggplot(data = SwabMetaData, aes(x=as.factor(AGE), y=MASS)) +geom_boxplot()
t.test(MASS~AGE, data = SwabMetaData)
ggplot(data = SwabMetaData, aes(x=as.factor(AGE), y=PLUM)) +geom_violin()
chisq.test(SwabMetaData$PLUM, SwabMetaData$AGE, correct = T)

#Mass per plumage
ggplot(data= SwabMetaData, aes(x=as.factor(PLUM), y=MASS)) +geom_boxplot()
ggplot(data= SwabMetaData, aes(x=as.factor(PLUM), y=MASS)) +geom_violin()
summary(aov(MASS ~PLUM, data=SwabMetaData))

#Bill length per plumage
ggplot(data= SwabMetaData, aes(x=as.factor(PLUM), y=BILL)) +geom_boxplot()
ggplot(data= SwabMetaData, aes(x=as.factor(PLUM), y=BILL)) + geom_violin()
summary(aov(BILL ~PLUM, data=SwabMetaData))


#Bill length vs mass
ggplot(data= SwabMetaData, aes(x=BILL, y=MASS)) + geom_point()
cor.test(SwabMetaData$BILL, SwabMetaData$MASS, method="pearson")

#Catch years and months
ggplot(data=SwabMetaData, aes(x=CATCH_YEAR, fill=as.factor(CATCH_MONTH))) +geom_bar()

#Catch year and mass
ggplot(data=SwabMetaData, aes(x=as.factor(CATCH_YEAR), y=MASS)) +geom_boxplot()
t.test(MASS~CATCH_YEAR, data = SwabMetaData[SwabMetaData$CATCH_YEAR %in% c(2019, 2021),]) # <-- should not be a t-test because we have >2 groups!!


#=================================================== Plotting Extraction metadata ===================================================

# Amount of material spread over extraction batches 
ggplot(data = ExtrMetaData, aes(x=Amount_of_material, fill=Extraction_batch)) + geom_bar() + scale_fill_brewer(palette="Set2")
summary(aov(Amount_of_material ~ Extraction_batch, data = ExtrMetaData[complete.cases(ExtrMetaData["Amount_of_material"]),]))  #complete.cases removes rows that have a NA for the amount of material, since these cannot be used in the anova

# Amount of material and being frozen or not
ggplot(data = ExtrMetaData, aes(x=Amount_of_material, fill=Swab_frozen_during_transfer)) + geom_bar()
summary(aov(Amount_of_material ~ Swab_frozen_during_transfer, data = ExtrMetaData[complete.cases(ExtrMetaData["Amount_of_material"]),]))  #complete.cases removes rows that have a NA for the amount of material, since these cannot be used in the anova




# ============================================== Free coding ================================================================

#Doing a t-test on the plumages and weights of catchyears 2019 and 2021, since they were during the same season

t.test(MASS~CATCH_YEAR, data = SwabMetaData[SwabMetaData$CATCH_YEAR %in% c(2019, 2021),])
chisq.test(SwabMetaData$PLUM, as.factor(SwabMetaData$CATCH_YEAR), correct = T)#, data=SwabMetaData[SwabMetaData$CATCH_YEAR %in% c(2019, 2021),])

print(SwabMetaData[SwabMetaData$CATCH_YEAR %in% c(2019, 2021),]) #Blijkbaar is opeens iedereen catch_year 2022

print(SwabMetaData[1,2])

str(SwabMetaData)
