#### Master Script 1: Clinical Data Extraction and Examination ####
# Decoding Quantitative Motor Features for Classification and Prediction
# in Severe Acquired Brain Injury
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# Department of Biomedical Engineering
# Department of Applied Mathematics and Statistics
# Whiting School of Engineering, Johns Hopkins University
# email address: shubhayu@jhu.edu

# https://cran.r-project.org/web/packages/table1/vignettes/table1-examples.html
require(table1)

data <- read.csv(file = "patient_clinical_data.csv")
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

