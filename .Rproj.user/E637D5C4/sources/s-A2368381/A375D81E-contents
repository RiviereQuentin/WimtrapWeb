#' printConfusionMatrix 
#'
#' @description A function that prints the confusion matrix obtained with a TFBS model
#'
#' @return Confusion matrix
#'
#' @noRd
printConfusionMatrix <- function(TFBSmodel){
	print(caret::confusionMatrix(as.factor(TFBSmodel$xgbpred), as.factor(TFBSmodel$ts_label)))
}