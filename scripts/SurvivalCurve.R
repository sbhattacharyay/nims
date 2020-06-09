setwd("C:/Users/Matt/OneDrive/OneDrive - Johns Hopkins University/Desktop/BIMS/Survival Curve/");
data <- read.csv(file = "patient_clinical_data.csv")

death12months <- data[, "Death12Months"]
deathdate <- data[, "DeathDate"]

missing_death12months <- is.nan(death12months)
missing_deathdate <- deathdate == "NaT"
havedeath12months_missingdate <- ((death12months == 1) + missing_deathdate) == 2
havedeath12months_missingdate[is.na(havedeath12months_missingdate)] <- FALSE

index = (missing_death12months + havedeath12months_missingdate) == 0
data <- data[index,c("HospitalDischargeDate", "DiedDuringThisHospitalStay_", "DeathDate", "Death12Months")]


event_time <- c()
for(i in seq(from = 1, to = nrow(data), by = 1)) {
    a <- as.character(data$HospitalDischargeDate[i])
    a <- strsplit(a,"-")
    a[[1]][2] <- match(a[[1]][2],month.abb)
    a[[1]][3] <- paste("20",a[[1]][3], sep = "")
    a <- paste(a[[1]][1], a[[1]][2], a[[1]][3], sep = "/")
    a <- as.numeric(as.Date(a, "%d/%m/%Y"))
    
    
    b <- as.character(data$DeathDate[i])
    if (b == "NaT") {
        b <- NA
        diff <- 365
    } else {
        b = strsplit(b,"-")
        b[[1]][2] <- match(b[[1]][2],month.abb)
        b[[1]][3] <- paste("20",b[[1]][3], sep = "")
        b <- paste(b[[1]][1], b[[1]][2], b[[1]][3], sep = "/")
        b <- as.numeric(as.Date(b, "%d/%m/%Y"))
        diff <- b-a
        if (diff > 365) {
            diff <- 365
        }
    }
    event_time <- c(event_time,diff)
}

data$event_time <- event_time


library(survival)
library(survminer)
library(dplyr)

surv_object <- Surv(time = data$event_time, event = data$Death12Months)
fit1 <- survfit(surv_object ~ 1, data = data)

ggsurvplot(fit1, data = data)

