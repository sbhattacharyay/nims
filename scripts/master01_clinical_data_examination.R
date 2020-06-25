#### Master Script 1: Clinical Data Extraction and Examination ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# Source: https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html
require(table1)

patient_clinical_data <- load_patient_clinical_data('../clinical_data/patient_clinical_data.csv',gose_thresh,mrs_thresh)

data$death <- factor(data$death, levels=c(0,1,2), labels=c("Alive at Discharge", "Died During Hospital Stay", "P-Value"))


###
data$gender <- factor(data$gender, levels=c("M", "F"), labels=c("Male", "Female"))
data$stroke <- factor(data$stroke, levels=c(0, 1), labels=c("0", "1"))
data$ich <- factor(data$ich, levels=c(0, 1), labels=c("0", "1"))
data$sah <- factor(data$sah, levels=c(0, 1), labels=c("0", "1"))
data$bt <- factor(data$bt, levels=c(0, 1), labels=c("0", "1"))
data$sdh <- factor(data$sdh, levels=c(0, 1), labels=c("0", "1"))
data$tbi <- factor(data$tbi, levels=c(0, 1), labels=c("0", "1"))


###
label(data$gender) <- "Sex"
label(data$stroke) <- "CVA"
label(data$ich) <- "ICH"
label(data$sah) <- "SAH"
label(data$bt) <- "BrainTumorOrLesion"
label(data$sdh) <- "SDH"
label(data$tbi) <- "TBI"
label(data$gose) <- "GOSE Discharge"
label(data$hlm_en) <- "HLM Enrollment"
label(data$hlm_dis) <- "HLM Discharge"
label(data$gcs_en) <- "GCS Enrollment"
label(data$gcs_dis) <- "GCS Discharge"
label(data$gcs_motor_en) <- "GCS Motor Enrollment"
label(data$gcs_motor_dis) <- "GCS Motor Discharge"

label(data$los) <- "Length of Stay in NCCU"
units(data$los) <- "d"


rndr <- function(x, name, ...) {
    if (length(x) == 0) {
        y <- data[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        if (is.numeric(y)) {
            p <- t.test(y ~ data$death)$p.value
        } else {
            p <- chisq.test(table(y, droplevels(data$death)))$p.value
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


table1(~ gender + stroke + ich + sah + sdh + tbi + + gose + gcs_en + gcs_dis + gcs_motor_en + gcs_motor_dis + hlm_en + hlm_dis + los| death, 
       data = data, overall = F, droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass = "Rtable1-grid")

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

library(devtools)
if (!require(lolR)) install_github('neurodata/lol', build_vignettes=TRUE, force=TRUE)
library(lolR)
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(plotly)
library(naniar)
library(MASS)

source('./functions/load_patient_clinical_data.R')
source('./functions/update_clinicalVariableList.R')
source('./functions/get_motion_features.R')
source('./functions/lol_project_motion_features.R')
source('./functions/viz_lol_2D.R')
source("./functions/generateRootDir.R")
source("./functions/plot_descriptive_figs.R")
    
# Load patient clinical data (sorts by PY numbering and corrects variable types)
patient_clinical_data<-load_patient_clinical_data()

# Load and update clinical variable list
clinicalVariableList<-update_clinicalVariableList()

plotDescriptiveFigs(patient_clinical_data,clinicalVariableList)