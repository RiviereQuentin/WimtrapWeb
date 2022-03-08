#' plotROC 
#'
#' @description A function that plots the ROC curve obtained with a TFBS model
#'
#' @return ROC plot
#'
#' @noRd
plotROC <- function(TFBSmodel){
	plot(pROC::roc(TFBSmodel$ts_label, TFBSmodel$xgbpred), xlim = c(1, 0), ylim = c(0, 1))
    lines(pROC::roc(TFBSmodel$ts_label, TFBSmodel$test_matchLogPval), col = "red")
    legend(x = "bottom", legend = c("Model", "Pattern-Matching"), fill = c("black", "red"))
    mtext(text = paste0("AUC: ", round(pROC::auc(TFBSmodel$ts_label, TFBSmodel$xgbpred), 2)), side = 2)
}