#### Master Script 1: Clinical Data Extraction and Examination ####
#
# Shubhayu Bhattacharyay, Matthew Wang
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

library(tidyverse)
library(table1)
library(devtools)
library(plotly)
library(naniar)
library(survival)
library(survminer)

source('./functions/load_patient_clinical_data.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source('./functions/viz_lol_2D.R')
source("./functions/generateRootDir.R")

# Load patient clinical data (sorts by PY numbering and corrects variable types)
data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv')
data$fav_mort_dis <- factor(data$fav_mort_dis, levels=c("Fav","Unfav",2), labels=c("Alive at Discharge", "Died During Hospital Stay", "P-Value"))


# Recode variables for clinical characteristics chart
data$gender <- factor(data$Sex, levels=c("M", "F"), labels=c("Male", "Female"))
data$stroke <- factor(data$CVA, levels=c(0, 1), labels=c("0", "1"))
data$ICH <- factor(data$ICH, levels=c(0, 1), labels=c("0", "1"))
data$SAH <- factor(data$SAH, levels=c(0, 1), labels=c("0", "1"))
data$BrainTumorOrLesion <- factor(data$BrainTumorOrLesion, levels=c(0, 1), labels=c("0", "1"))
data$SDHOrEDH <- factor(data$SDHOrEDH, levels=c(0, 1), labels=c("0", "1"))
data$TBI <- factor(data$TBI, levels=c(0, 1), labels=c("0", "1"))
data$fav_mRS_dis <- factor(data$fav_mRS_dis, levels=c("Unfav", "Fav"), labels=c("0", "1"))
data$fav_GOSE_dis <- factor(data$fav_GOSE_dis, levels=c("Unfav", "Fav"), labels=c("0", "1"))


# Add labels to dataframe variable names
label(data$Sex) <- "Sex"
label(data$CVA) <- "CVA"
label(data$ICH) <- "ICH"
label(data$SAH) <- "SAH"
label(data$BrainTumorOrLesion) <- "BrainTumorOrLesion"
label(data$SDHOrEDH) <- "SDH"
label(data$TBI) <- "TBI"
label(data$fav_GOSE_dis) <- "Favorable GOSE @Discharge"
label(data$fav_mRS_dis) <- "Favorable mRS @Discharge"
label(data$DaysInNCCU) <- "Length of Stay in NCCU"
units(data$DaysInNCCU) <- "d"

# Establish render functions for clinical chart
rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- data[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ data$fav_mort_dis)$p.value
        } else {
            p <- chisq.test(table(y, droplevels(data$fav_mort_dis)))$p.value
        }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
    } else {
        render.default(x=x, name=name, ...)
    }
}
rndr.strat <- function(label, n, ...) {
    ifelse(n==0, label, render.strat.default(label, n, ...))
}

# Produce clinical characteristics chart
table1(~ Sex + CVA + ICH + SAH + BrainTumorOrLesion + SDHOrEDH + TBI + fav_GOSE_dis + fav_mRS_dis + DaysInNCCU | fav_mort_dis, 
       data = data, overall = F, droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass = "Rtable1-grid")

# Format death labels for survival plot
death12months <- data[, "fav_mort_12m"]
death12months <- factor(death12months, levels = c("Fav","Unfav"), labels = c(0, 1))
death12months2 <- as.numeric(levels(death12months))[death12months]
data$fav_mort_12m <- death12months2
deathdate <- data[, "DeathDate"]

missing_death12months <- is.na(death12months)
missing_deathdate <- is.na(deathdate)
havedeath12months_missingdate <- ((death12months == 1) + missing_deathdate) == 2
havedeath12months_missingdate[is.na(havedeath12months_missingdate)] <- FALSE

index = (missing_death12months + havedeath12months_missingdate) == 0
data2 <- data[index,c("HospitalDischargeDate", "DiedDuringThisHospitalStay_", "DeathDate", "fav_mort_12m")]

event_time <- c()
for(i in seq(from = 1, to = nrow(data2), by = 1)) {
    a <- as.character(data2$HospitalDischargeDate[i])
    a <- strsplit(a,"-")
    a <- paste(a[[1]][3], a[[1]][2], a[[1]][1], sep = "/")
    a <- as.numeric(as.Date(a, "%d/%m/%Y"))
    b <- as.character(data2$DeathDate[i])
    if (is.na(b)) {
        b <- NA
        diff <- 365
    } else {
        b = strsplit(b,"-")
        b <- paste(b[[1]][3], b[[1]][2], b[[1]][1], sep = "/")
        b <- as.numeric(as.Date(b, "%d/%m/%Y"))
        diff <- b-a
        if (diff > 365) {
            diff <- 365
        }
    }
    event_time <- c(event_time,diff)
}
data2$event_time <- event_time

# Plot survival curve
surv_object <- Surv(time = data2$event_time, event = data2$fav_mort_12m)
fit1 <- survfit(surv_object ~ 1, data = data2)
ggsurvplot(fit1, data = data2)