#### Train regression models ####
#
# Shubhayu Bhattacharyay, Matthew Wang, Eshan Joshi
# University of Cambridge
# Johns Hopkins University
# email address: sb2406@cam.ac.uk


# Load necessary packages
library(caret)
library(rms)
library(tidyverse)
library(nnet)
library(MASS)
library(UBL)
library(ROCR)
if(.Platform$OS.type == "unix") {
  library(doMC)
} else {
  library(doParallel)
}

# Set the number of parallel cores
no.parallel.cores <- floor(2 * detectCores() / 3)
if(.Platform$OS.type == "unix") {
  registerDoMC(cores = no.parallel.cores)
} else {
  registerDoParallel(cores = no.parallel.cores)
}

# Directories of imputations
imp.dirs <- list.files('../all_motion_feature_data/04_formatted_predictor_matrices',pattern = 'imp*',include.dirs = TRUE, full.names = TRUE)

# Load detection labels and partitions
load('../validation_resampling/detection_labels.RData')
load('../validation_resampling/detection_partitions.RData')

# Dimension choices
r <- c(2,3,4,5,10,20,30,40,50)
source('./functions/macro_average_AUROC.R')

# Initialize empty dataframe to store compiled tuning results
compiled.tuning.df <- as.data.frame(matrix(ncol = 5, nrow = 0))
for (curr.r in r){
  # Iterate through imputations and identify detection folders
  for (i in 1:length(imp.dirs)){
    print(paste("Imputation no.",i,"out of",length(imp.dirs),"started."))
    curr.imp.dirs <- imp.dirs[i]
    detection.folders <- list.files(path = curr.imp.dirs,pattern = 'detection_*',include.dirs = TRUE, full.names = TRUE)
    
    # Create directory for imputation
    dir.create(path = file.path("../results/detection_results/regression_results",paste0("imp",i)),showWarnings = FALSE,recursive = TRUE) 
    
    # Iterate through detection folders and tune, train, and validate models
    for (j in 1:length(detection.folders)){
      
      # Extract current detection folder information
      curr.window.size <- as.numeric(sub(".*window_", "", detection.folders[j]))
      curr.window.idx <- which(det_parameters$obs_windows == curr.window.size)
      curr.label.set <- det_gcs_labels[[curr.window.idx]]
      
      # Create directory for current detection paradigm
      dir.create(paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size),
                 showWarnings = FALSE,
                 recursive = TRUE)
      
      print(paste("Detection folder no.",j,"out of",length(detection.folders),"started."))
      
      # Load current training matrix,scale columns, and rename variables to match formula
      curr.motor.train.matrix <- as.data.frame(readRDS(file.path(detection.folders[j],'motor_train_matrix_lol2.rds')))[,1:curr.r] %>%
        scale() %>%
        as.data.frame()
      names(curr.motor.train.matrix) <- paste(paste0(rep("MF.",ncol(curr.motor.train.matrix)),1:ncol(curr.motor.train.matrix)))
      curr.motor.train.matrix$GCSm <- as.factor(curr.label.set$Best.Motor.Response[det_motor_train_idx[[curr.window.idx]]])
      
      # Perform SMOTE on training test
      curr.motor.smote.matrix <- SmoteClassif(GCSm ~., dat = curr.motor.train.matrix) 
      curr.motor.smote.matrix[,1:curr.r] <- scale(curr.motor.smote.matrix[,1:curr.r]) 
      
      # Load current testing matrix and rename variables to match formula
      curr.motor.test.matrix <- as.data.frame(readRDS(file.path(detection.folders[j],'motor_test_matrix_lol2.rds')))[,1:curr.r] %>%
        scale() %>%
        as.data.frame()        
      names(curr.motor.test.matrix) <- paste(paste0(rep("MF.",ncol(curr.motor.test.matrix)),1:ncol(curr.motor.test.matrix)))
      curr.motor.test.matrix$GCSm <- as.factor(curr.label.set$Best.Motor.Response[det_motor_test_idx[[curr.window.idx]]])
      
      # Define regression formula
      regression.formula <- as.formula(paste("GCSm","~",paste(paste0(rep("MF.",curr.r),1:(curr.r)),collapse = " + ")))
      
      ## ORDINAL REGRESSION
      # Train ordinal regression model without SMOTE
      curr.polr.mdl <- orm(regression.formula, 
                           data = curr.motor.train.matrix, 
                           model=TRUE,
                           linear.predictors=TRUE,
                           se.fit=TRUE,
                           maxit=500,
                           scale=TRUE)
      saveRDS(object = curr.polr.mdl, 
              file = paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/polr_mdl_',curr.r,'.rds'))
      
      # Test model on validation set
      curr.polr.mdl.preds <- cbind(data.frame(true.labels = curr.label.set$Best.Motor.Response[det_motor_test_idx[[curr.window.idx]]],
                                              pred.labels = max.col(predict(curr.polr.mdl,curr.motor.test.matrix,type='fitted.ind'),ties.method="first")),
                                   predict(curr.polr.mdl,curr.motor.test.matrix,type='fitted.ind'))
      write.csv(curr.polr.mdl.preds,
                paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/polr_mdl_',curr.r,'_preds.csv'),
                row.names = FALSE)
      
      compiled.tuning.df <- rbind(compiled.tuning.df,
                                  data.frame(model = 'polr',
                                             r = curr.r,
                                             aurocs = macro.average.AUROC(curr.polr.mdl.preds),
                                             obs.window = curr.window.size,
                                             imp = i))
      
      # Train ordinal regression model with SMOTE
      curr.polr.smote.mdl <- orm(regression.formula, 
                                 data = curr.motor.smote.matrix, 
                                 model=TRUE,
                                 linear.predictors=TRUE,
                                 se.fit=TRUE,
                                 maxit=500,
                                 scale=TRUE)
      saveRDS(object = curr.polr.smote.mdl, 
              file = paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/SMOTE_polr_mdl_',curr.r,'.rds'))
      
      # Test model on validation set
      curr.polr.smote.mdl.preds <- cbind(data.frame(true.labels = curr.label.set$Best.Motor.Response[det_motor_test_idx[[curr.window.idx]]],
                                                    pred.labels = max.col(predict(curr.polr.smote.mdl,curr.motor.test.matrix,type='fitted.ind'),ties.method="first")),
                                         predict(curr.polr.smote.mdl,curr.motor.test.matrix,type='fitted.ind'))
      write.csv(curr.polr.smote.mdl.preds,
                paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/SMOTE_polr_mdl_',curr.r,'_preds.csv'),
                row.names = FALSE)
      
      compiled.tuning.df <- rbind(compiled.tuning.df,
                                  data.frame(model = 'polr.smote',
                                             r = curr.r,
                                             aurocs = macro.average.AUROC(curr.polr.smote.mdl.preds),
                                             obs.window = curr.window.size,
                                             imp = i))
      
      # Demarcate baseline outcome for multinomial fit
      curr.motor.train.matrix$GCSm <- relevel(curr.motor.train.matrix$GCSm, ref = "1")
      curr.motor.smote.matrix$GCSm <- relevel(curr.motor.smote.matrix$GCSm, ref = "1")
      
      # Train multinomial regression model without SMOTE
      curr.mnlr.mdl <- multinom(regression.formula,
                                data = curr.motor.train.matrix,
                                Hess = TRUE,
                                model = TRUE,
                                trace = FALSE)
      curr.mnlr.mdl$z.scores <- summary(curr.mnlr.mdl)$coefficients/summary(curr.mnlr.mdl)$standard.errors
      curr.mnlr.mdl$p.values <- (1 - pnorm(abs(curr.mnlr.mdl$z.scores), 0, 1)) * 2
      saveRDS(object = curr.mnlr.mdl, 
              file = paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/mnlr_mdl_',curr.r,'.rds'))
      
      # Test model on validation set
      curr.mnlr.mdl.preds <- cbind(data.frame(true.labels = curr.label.set$Best.Motor.Response[det_motor_test_idx[[curr.window.idx]]],
                                              pred.labels = predict(curr.mnlr.mdl, curr.motor.test.matrix)),
                                   as.data.frame(predict(curr.mnlr.mdl, curr.motor.test.matrix,type = 'p')) %>% 
                                     rename('GCSm=1' = '1','GCSm=2' = '2','GCSm=3' = '3','GCSm=4' = '4','GCSm=5' = '5','GCSm=6' = '6'))
      write.csv(curr.mnlr.mdl.preds,
                paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/mnlr_mdl_',curr.r,'_preds.csv'),
                row.names = FALSE)
      
      compiled.tuning.df <- rbind(compiled.tuning.df,
                                  data.frame(model = 'mnlr',
                                             r = curr.r,
                                             aurocs = macro.average.AUROC(curr.mnlr.mdl.preds),
                                             obs.window = curr.window.size,
                                             imp = i))
      
      # Train multinomial regression model with SMOTE
      curr.mnlr.smote.mdl <- multinom(regression.formula,
                                      data = curr.motor.smote.matrix,
                                      Hess = TRUE,
                                      model = TRUE,
                                      trace = FALSE)
      curr.mnlr.smote.mdl$z.scores <- summary(curr.mnlr.smote.mdl)$coefficients/summary(curr.mnlr.smote.mdl)$standard.errors
      curr.mnlr.smote.mdl$p.values <- (1 - pnorm(abs(curr.mnlr.smote.mdl$z.scores), 0, 1)) * 2
      saveRDS(object = curr.mnlr.smote.mdl, 
              file = paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/SMOTE_mnlr_mdl_',curr.r,'.rds'))
      
      # Test model on validation set
      curr.mnlr.smote.mdl.preds <- cbind(data.frame(true.labels = curr.label.set$Best.Motor.Response[det_motor_test_idx[[curr.window.idx]]],
                                                    pred.labels = predict(curr.mnlr.smote.mdl, curr.motor.test.matrix)),
                                         as.data.frame(predict(curr.mnlr.smote.mdl, curr.motor.test.matrix,type = 'p')) %>% 
                                           rename('GCSm=1' = '1','GCSm=2' = '2','GCSm=3' = '3','GCSm=4' = '4','GCSm=5' = '5','GCSm=6' = '6'))
      write.csv(curr.mnlr.smote.mdl.preds,
                paste0('../results/detection_results/regression_results/imp',i,'/detection_window_',curr.window.size,'/SMOTE_mnlr_mdl_',curr.r,'_preds.csv'),
                row.names = FALSE)
      
      compiled.tuning.df <- rbind(compiled.tuning.df,
                                  data.frame(model = 'mnlr.smote',
                                             r = curr.r,
                                             aurocs = macro.average.AUROC(curr.mnlr.smote.mdl.preds),                                             
                                             obs.window = curr.window.size,
                                             imp = i))      
    }
  }
}

# Save compiled tuning results
write.csv(compiled.tuning.df, '../results/detection_results/regression_results/compiled_tuning_results.csv')

opt.r.df <- compiled.tuning.df %>%
  group_by(model,r,obs.window) %>%
  summarise(meanAUROCS = mean(aurocs),sdAUROCS = sd(aurocs),n = n()) %>%
  group_by(model,obs.window) %>%
  summarise(maxAUROCdim = r[which.max(meanAUROCS)],maxAUROC = max(meanAUROCS))

# Save optimal tuning results as a CSV
write.csv(opt.r.df,'../results/detection_results/regression_results/optimal_tuning_results.csv')

# Imputation directories
imp.dirs <- list.files('../results/detection_results/regression_results',pattern = '^imp*',include.dirs = TRUE,full.names = TRUE)

# Create new directory to store tuned model metrics
dir.create('../results/detection_results/regression_results/tuned_model_metrics',showWarnings = FALSE)

# Loop through optimally tuned model results and calculate a bevy of classification metrics
super.compiled.metrics.df <- as.data.frame(matrix(ncol= 6, nrow = 0))
for (i in 1:nrow(opt.r.df)){
  
  if (opt.r.df$model[i] == 'mnlr'){
    model.file.name <- paste0('mnlr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'mnlr.smote') {
    model.file.name <- paste0('SMOTE_mnlr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'polr') {
    model.file.name <- paste0('polr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'polr.smote') {
    model.file.name <- paste0('SMOTE_polr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  }
   
  compiled.metrics.df <- as.data.frame(matrix(ncol=74,nrow = 0))
   
  for (j in 1:length(imp.dirs)){
    
    curr.mdl.results <-
      read.csv(file.path(
        '../results/detection_results/regression_results',
        paste0('imp', j),
        paste0('detection_window_',opt.r.df$obs.window[i]),
        model.file.name
      )) %>%
      mutate(true.labels = factor(true.labels, levels = 1:6),
             pred.labels = factor(pred.labels, levels = 1:6))
    
      curr.mdl.cm <- caret::confusionMatrix(curr.mdl.results$pred.labels,curr.mdl.results$true.labels)
      no.overall <- length(curr.mdl.cm$overall)
      dim.byClass <- dim(curr.mdl.cm$byClass)
      
      curr.metrics.df <- as.data.frame(matrix(nrow = 1, ncol = no.overall + (dim.byClass[1]*dim.byClass[2])))
      names(curr.metrics.df)[1:no.overall] <- names(curr.mdl.cm$overall)
      names(curr.metrics.df)[(no.overall+1):ncol(curr.metrics.df)] <- paste(rep(colnames(curr.mdl.cm$byClass), each = length(rownames(curr.mdl.cm$byClass))), rownames(curr.mdl.cm$byClass), sep = ".")
      
      curr.metrics.df[1,names(curr.mdl.cm$overall)] <- curr.mdl.cm$overall[names(curr.mdl.cm$overall)]
      curr.metrics.df[1,(no.overall+1):ncol(curr.metrics.df)] <- as.vector(curr.mdl.cm$byClass)
      curr.metrics.df$imp <- j
      
      compiled.metrics.df <- rbind(compiled.metrics.df,curr.metrics.df)
  }
  
  compiled.metrics.df <- compiled.metrics.df %>% 
    pivot_longer(cols = -imp, names_to = 'metric') %>%
    mutate(model = opt.r.df$model[i],
           obs.window = opt.r.df$obs.window[i],
           opt.r = opt.r.df$maxAUROCdim[i])

  write.csv(
    compiled.metrics.df,
    file.path(
      '../results/detection_results/regression_results/tuned_model_metrics',
      paste0(opt.r.df$model[i], '_', opt.r.df$obs.window[i], '.csv')
    ),
    row.names = FALSE
  )
  super.compiled.metrics.df <- rbind(super.compiled.metrics.df, compiled.metrics.df)
}

# Save super compiled metrics dataframe
write.csv(super.compiled.metrics.df, '../results/detection_results/regression_results/tuned_model_metrics/compiled_metrics.csv')

grouped.metrics <- super.compiled.metrics.df %>% group_by(obs.window,model,metric) %>%
  summarise(meanValue = mean(value), sdValue = sd(value))

# Loop through optimally tuned model results and calculate a bevy of probability metrics and axes

# Axes for interpolation
xq <- seq(from = 0, to = 1, length.out = 4000)
compiled.curve.axes <- data.frame(matrix(ncol = 9, nrow = 0))
compiled.aucs <- data.frame(matrix(ncol = 7, nrow = 0))
  
for (i in 1:nrow(opt.r.df)){
  
  print(paste("Model no.",i,"out of",nrow(opt.r.df),"started."))
  
  if (opt.r.df$model[i] == 'mnlr'){
    model.file.name <- paste0('mnlr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'mnlr.smote') {
    model.file.name <- paste0('SMOTE_mnlr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'polr') {
    model.file.name <- paste0('polr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  } else if (opt.r.df$model[i] == 'polr.smote') {
    model.file.name <- paste0('SMOTE_polr_mdl_',opt.r.df$maxAUROCdim[i],'_preds.csv')
  }
  for (j in 1:length(imp.dirs)){
    print(paste("Imputation no.",j,"out of",length(imp.dirs),"started."))
    curr.mdl.results <-
      read.csv(file.path(
        '../results/detection_results/regression_results',
        paste0('imp', j),
        paste0('detection_window_',opt.r.df$obs.window[i]),
        model.file.name
      )) %>%
      mutate(true.labels = factor(true.labels, levels = 1:6),
             pred.labels = factor(pred.labels, levels = 1:6))
    
    unique.labels <- sort(unique(curr.mdl.results$true.labels))
    
    for (k in 1:length(unique.labels)){
      print(paste("GCSm",k,"started."))
      curr.label <- unique.labels[k]
      temp.results.df <- curr.mdl.results %>% 
        mutate(temp.label = as.numeric(true.labels == curr.label))
      curr.label.prob.name <- names(temp.results.df)[grep(paste('GCSm',curr.label,sep = "."),names(temp.results.df))]
      
      curr.pred <- prediction(temp.results.df[,curr.label.prob.name],temp.results.df$temp.label)
      
      curr.perf.roc <- performance(curr.pred, measure = "tpr", x.measure = "fpr")
      curr.perf.prc <- performance(curr.pred, measure = "prec", x.measure = "rec")
      
      curr.perf.auroc <- performance(curr.pred, measure = "auc")
      curr.perf.auprc <- performance(curr.pred, measure = "aucpr")
      
      interp.roc <- approx(x = curr.perf.roc@x.values[[1]], y = curr.perf.roc@y.values[[1]],xout = xq)
      interp.prc <- approx(x = curr.perf.prc@x.values[[1]], y = curr.perf.prc@y.values[[1]],xout = xq)
      
      compiled.curve.axes <- rbind(compiled.curve.axes,
                                   data.frame(fpr = interp.roc$x, 
                                              tpr = interp.roc$y, 
                                              rec = interp.prc$x, 
                                              prec = interp.prc$y,
                                              GCSm = curr.label,
                                              model = opt.r.df$model[i], 
                                              obs.window = opt.r.df$obs.window[i],
                                              opt.r = opt.r.df$maxAUROCdim[i],
                                              imp = j))

      compiled.aucs <- rbind(compiled.aucs,
                             data.frame(auroc = curr.perf.auroc@y.values[[1]], 
                                        auprc = curr.perf.auprc@y.values[[1]], 
                                        GCSm = curr.label,
                                        model = opt.r.df$model[i], 
                                        obs.window = opt.r.df$obs.window[i],
                                        opt.r = opt.r.df$maxAUROCdim[i],
                                        imp = j))
    }
  }
}

# Save compiled dataframes
write.csv(compiled.curve.axes,'../results/detection_results/regression_results/tuned_model_metrics/compiled_curve_axes.csv',row.names = FALSE)
write.csv(compiled.aucs,'../results/detection_results/regression_results/tuned_model_metrics/compiled_aucs.csv',row.names = FALSE)

# Average compiled axes dataframe across imputations
group.curve.axes <- compiled.curve.axes %>%
  group_by(model,obs.window,GCSm,fpr,rec) %>%
  summarise(meanTPR = mean(tpr,na.rm = TRUE), meanPREC = mean(prec,na.rm = TRUE), sdTPR = sd(tpr,na.rm = TRUE), sdPREC = sd(prec,na.rm = TRUE)) %>%
  mutate(upperTPR = min(meanTPR + sdTPR, 1),lowerTPR = max(meanTPR - sdTPR, 0),upperPREC = min(meanPREC + sdPREC, 1),lowerPREC = max(meanPREC - sdPREC, 0))

# Save averaged curve information
write.csv(group.curve.axes,'../results/detection_results/regression_results/tuned_model_metrics/averaged_curve_axes.csv',row.names = FALSE)

# Average compiled AUC dataframe across imputations
group.aucs <- compiled.aucs %>%
  group_by(model, obs.window, GCSm) %>%
  summarise(
    meanAUROC = mean(auroc, na.rm = TRUE),
    meanAUPRC = mean(auprc, na.rm = TRUE),
    sdAUROC = sd(auroc, na.rm = TRUE),
    sdAUPRC = sd(auprc, na.rm = TRUE)
  )

# Save averaged AUC information
write.csv(group.aucs,'../results/detection_results/regression_results/tuned_model_metrics/averaged_aucs.csv',row.names = FALSE)

group.curve.axes <- group.curve.axes %>% mutate(model.id = paste0(model,'_',obs.window))
group.aucs <- group.aucs %>% 
  mutate(model.id = paste0(model,'_',obs.window),
         model.type = str_sub(model,1,4))

# Determine whether SMOTE is optimal for each model
grouped.opt.r.df <- opt.r.df %>%
  mutate(model.name = str_sub(model,1,4),
         smote.ind = str_detect(model,'smote'),
         model.id = paste0(model,'_',obs.window)) %>%
  group_by(model.name,obs.window) %>%
  summarise(SMOTE.ind = smote.ind[which.max(maxAUROC)], 
            model.ident = model.id[which.max(maxAUROC)],
            d = maxAUROCdim[which.max(maxAUROC)],
            maxAUROC = max(maxAUROC)) %>%
  mutate(formatted = sprintf('%0.02f',maxAUROC))

# Table 4 statistics
table.4.stats <- group.aucs %>%
  group_by(GCSm,obs.window,model.type) %>%
  summarise(maxModel.id.AUROC = model.id[which.max(meanAUROC)],
            maxMeanAUROC = max(meanAUROC),
            maxSdAUROC = sdAUROC[which.max(meanAUROC)],
            maxModel.id.AUPRC = model.id[which.max(meanAUPRC)],
            maxMeanAUPRC = max(meanAUPRC),
            maxSdAUPRC = sdAUPRC[which.max(meanAUPRC)]) %>%
  mutate(formatted.1 = sprintf('%.02f (%.02f)',maxMeanAUROC,maxSdAUROC),
         formatted.2 = sprintf('%.02f (%.02f)',maxMeanAUPRC,maxSdAUPRC)) %>%
  dplyr::select(GCSm, obs.window, model.type, formatted.1, formatted.2) %>%
  pivot_wider(id_cols = c(GCSm, obs.window),names_from = model.type, values_from = c(formatted.1, formatted.2))

table.4.macro.averaged <- group.aucs %>%
  group_by(GCSm,obs.window,model.type) %>%
  summarise(maxModel.id.AUROC = model.id[which.max(meanAUROC)],
            maxMeanAUROC = max(meanAUROC),
            maxSdAUROC = sdAUROC[which.max(meanAUROC)],
            maxModel.id.AUPRC = model.id[which.max(meanAUPRC)],
            maxMeanAUPRC = max(meanAUPRC),
            maxSdAUPRC = sdAUPRC[which.max(meanAUPRC)]) %>%
  group_by(obs.window,model.type) %>%
  summarise(macroMeanAUROC = mean(maxMeanAUROC),
            macroSdAUROC = sqrt( (1/n())*sum(maxSdAUROC^2)  + ((n()+1)/(n()))*(sd(maxMeanAUROC)^2) ) ,
            macroMeanAUPRC = mean(maxMeanAUPRC),
            macroSdAUPRC = sqrt( (1/n())*sum(maxSdAUPRC^2)  + ((n()+1)/(n()))*(sd(maxMeanAUPRC)^2) ))%>%
  mutate(formatted.1 = sprintf('%.02f (%.02f)',macroMeanAUROC,macroSdAUROC),
         formatted.2 = sprintf('%.02f (%.02f)',macroMeanAUPRC,macroSdAUPRC)) %>%
  dplyr::select(obs.window, model.type, formatted.1, formatted.2) %>%
  pivot_wider(id_cols = c(obs.window),names_from = model.type, values_from = c(formatted.1, formatted.2))

# Determine optimal models for each observation window per each GCSm discrimination
opt.group.aucs <- group.aucs %>%
  group_by(obs.window, GCSm) %>%
  summarise(
    maxAUROC = max(meanAUROC),
    model.id = model.id[which.max(meanAUROC)],
    sdAUROC = sdAUROC[which.max(meanAUROC)]
  ) %>%
  mutate(formatted = sprintf('%.02f (%.02f)',maxAUROC,sdAUROC))

# Determine optimal models for each observation window per each GCSm discrimination
opt.group.aucprs <- group.aucs %>%
  group_by(obs.window, GCSm) %>%
  summarise(
    maxAUPRC = max(meanAUPRC),
    model.id = model.id[which.max(meanAUPRC)],
    sdAUPRC = sdAUPRC[which.max(meanAUPRC)]
  ) %>%
  mutate(formatted = sprintf('%.02f (%.02f)',maxAUPRC,sdAUPRC))

# Filter group.curve.axes for plotting ROCs
filtered.group.curve.axes <- as.data.frame(matrix(ncol = ncol(group.curve.axes), nrow = 0))
for (i in 1:nrow(opt.group.aucs)){
  print(paste("Curve no.",i,"out of",nrow(opt.group.aucs),"filtered."))
  temp.df <- group.curve.axes %>% 
    filter(model.id == opt.group.aucs$model.id[i] & GCSm == opt.group.aucs$GCSm[i])
  filtered.group.curve.axes <- rbind(filtered.group.curve.axes,temp.df)
}

# Save filtered, averaged axes for plotting ROCs
write.csv(filtered.group.curve.axes,'../results/detection_results/regression_results/tuned_model_metrics/plot_averaged_curve_axes.csv',row.names = FALSE)

# Filter group.curve.axes for plotting PRCs
filtered.group.prc.axes <- as.data.frame(matrix(ncol = ncol(group.curve.axes), nrow = 0))
for (i in 1:nrow(opt.group.aucprs)){
  print(paste("Curve no.",i,"out of",nrow(opt.group.aucprs),"filtered."))
  temp.df <- group.curve.axes %>% 
    filter(model.id == opt.group.aucprs$model.id[i] & GCSm == opt.group.aucprs$GCSm[i])
  filtered.group.prc.axes <- rbind(filtered.group.prc.axes,temp.df)
}

# Save filtered, averaged axes for plotting PRCs
write.csv(filtered.group.prc.axes,'../results/detection_results/regression_results/tuned_model_metrics/plot_averaged_prc_axes.csv',row.names = FALSE)