fix_test_cutpoint <- function(preds_set,self_cps_set){
  outputList <- list()
  for (i in 1:length(preds_set)){
    curr_pred_space <- preds_set[[i]]
    pred_space_cps <- list()
    n_folds <- length(curr_pred_space)
    for (j in 1:length(curr_pred_space)){
      curr_fold <- curr_pred_space[[j]]
      modelNames <- colnames(curr_fold)
      tru_levels <- levels(curr_fold[[1]])
      fold_cps <- list()
      for(m in 2:ncol(curr_fold)){
        curr_cp <- cutpointr(curr_fold[[m]],
                             curr_fold[[1]],
                             direction = ">=",
                             pos_class = tru_levels[2],
                             neg_class = tru_levels[1],
                             method = oc_manual,
                             cutpoint = self_cps_set[[i]][[j]][[m-1]]$optimal_cutpoint,
                             break_ties = min)
        fold_cps[[m-1]] = curr_cp
      }
      names(fold_cps) <- modelNames[2:length(modelNames)]
      pred_space_cps[[j]] <- fold_cps
    }
    names(pred_space_cps) <- paste0(rep("Fold",n_folds),1:n_folds)
    outputList[[i]] <- pred_space_cps
  }
  names(outputList) <- names(preds_set)
  return(outputList)
}