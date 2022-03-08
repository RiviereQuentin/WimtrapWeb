#' buildTFBSmodel 
#'
#' @description A modified version of the builTFBSmodel function of Wimtrap to allow the output of the ROC, 
#' feature importance plot and the confusion matrix
#'
#' @return A list of values including the TFBS model and values that are necessary for plotting
#'
#' @noRd
buildTFBSmodel <- function(TFBSdata,
                           ChIPpeaks,
                           ChIPpeaks_length = 400,
                           TFs_validation = NULL,
                           model_assessment = TRUE,
                           xgb_modeling = TRUE){
  ChIP_regions <- Wimtrap:::listChIPRegions(ChIPpeaks, NULL, ChIPpeaks_length)
  DataSet <- data.frame()
  if (length(names(ChIPpeaks))==0 &length(TFBSdata) == 1) {names(ChIPpeaks) == names(TFBSdata)}
  for (trainingTF in names(ChIPpeaks)){
     considered <- data.table::fread(TFBSdata[trainingTF],
                                     stringsAsFactors = TRUE)
     considered$TF <- trainingTF
     considered <- GenomicRanges::makeGRangesFromDataFrame(considered,
                                                           keep.extra.columns = TRUE)
     validated_TFBS <- GenomicRanges::findOverlaps(considered, ChIP_regions[[trainingTF]], select = "all")
     considered <- as.data.frame(considered)
     considered$ChIP.peak <- 0
     considered$ChIP.peak[validated_TFBS@from[!duplicated(validated_TFBS@from)]] <- 1
     NbTrueBs <- nrow(considered[considered$ChIP.peak == 1,])
     DataSet <- rbind(DataSet,
                      considered[considered$ChIP.peak == 1,],
                      considered[sample(which(considered$ChIP.peak == 0), NbTrueBs),])
     rm(considered)
  }
  DataSet <- data.table::as.data.table(DataSet)
  TFBSdata <- DataSet[,-seq(1,5), with = FALSE]
  rm(DataSet)
  #Split the dataset into a training and a validation datasets
  if(is.null(TFs_validation)){
    trainind <- sample(seq(1,nrow(TFBSdata)), as.integer(nrow(TFBSdata)*0.8))
    testind <- seq(1,nrow(TFBSdata))[!(seq(1,nrow(TFBSdata)) %in% trainind)]
  } else {
    trainind <- which(TFBSdata$TF %in% TFs_validation)
    testind <- which(!(TFBSdata$TF) %in% TFs_validation)
  }
  TFBSdata.training <- TFBSdata[trainind,]
  TFBSdata.validation <- TFBSdata[testind,]
  #Pre-processing of the training dataset
  ##Remove columns that do not have to be taken into account
  if (length(grep(pattern = "matchScore", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "matchScore", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "TF", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "TF", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "ClosestTSS", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "ClosestTSS", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "ClosestTTS", colnames(TFBSdata.training))) > 0){
    torm <- grep(pattern = "ClosestTTS", colnames(TFBSdata.training))
    TFBSdata.training <- TFBSdata.training[, colnames(TFBSdata.training)[torm] := NULL]
  } else {}
  ## Remove the infinite p-values associated to P-M,
  ## that occurs when the PWM is not flexible (i.e. is a consensus)
  TFBSdata.training$matchLogPval[which(is.infinite(TFBSdata.training$matchLogPval))] <- max(TFBSdata$matchLogPval)
  ##Create dummy variables
  dummy <- stats::model.matrix(~.+0, data = TFBSdata.training[,-c("ChIP.peak"),with=F])
  TFBSdata.training <- cbind(dummy, TFBSdata.training[,"ChIP.peak"])
  ##Remove highly correlated features
  descrCor <- stats::cor(TFBSdata.training, use = 'complete')
  highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .95)
  filteredDescr <- TFBSdata.training[,colnames(TFBSdata.training)[highlyCorDescr] := NULL]
  NAs <- is.na(as.data.frame(filteredDescr))
  NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
  train <- filteredDescr[!NAs,]
  #Pre-processing of the validation dataset
  ##Remove columns that do not have to be taken into account
  if (length(grep(pattern = "TF", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TF", colnames(TFBSdata.validation))
     TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
   } else {}
  if (length(grep(pattern = "TSS", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TSS", colnames(TFBSdata.validation))
    TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
  } else {}
  if (length(grep(pattern = "TTS", colnames(TFBSdata.validation))) > 0){
    torm <- grep(pattern = "TTS", colnames(TFBSdata.validation))
    TFBSdata.validation <- TFBSdata.validation[, colnames(TFBSdata.validation)[torm] := NULL]
  } else {}
  ## Remove the infinite p-values associated to P-M,
  ## that occurs when the PWM is not flexible (i.e. is a consensus)
  TFBSdata.validation$matchLogPval[which(is.infinite(TFBSdata.validation$matchLogPval))] <- max(TFBSdata$matchLogPval)
  #Create dummy variables
  dummy <- stats::model.matrix(~.+0, data = TFBSdata.validation[,-c("ChIP.peak"),with=F])
  TFBSdata.validation <- cbind(dummy, TFBSdata.validation[,"ChIP.peak"])
  ##Remove highly correlated features
  descrCor <- stats::cor(TFBSdata.validation, use = 'complete')
  highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .95)
  filteredDescr <- TFBSdata.validation[,colnames(TFBSdata.validation)[highlyCorDescr] := NULL]
  NAs <- is.na(as.data.frame(filteredDescr))
  NAs <- apply(NAs, 1, function(x) {if (length(which(x==TRUE)) > 0 ) {return(TRUE)} else {return(FALSE)} })
  test <- filteredDescr[!NAs,]

  train <- train[,which(colnames(train) %in% colnames(test)), with = FALSE]
  test <- test[,which(colnames(test) %in% colnames(train)), with = FALSE]
  #Build a model by extreme gradient boosting
  labels <- train$ChIP.peak
  ts_label <- test$ChIP.peak

  new_tr <- stats::model.matrix(~.+0,data = train[,-c("ChIP.peak"),with=F])
  new_ts <- stats::model.matrix(~.+0,data = test[,-c("ChIP.peak"),with=F])

  if (xgb_modeling == TRUE) {
    rm(train)
    rm(filteredDescr)
    rm(descrCor)
    rm(highlyCorDescr)
    rm(TFBSdata.validation)
    rm(TFBSdata.training)
    rm(NAs)
    dtrain <- xgboost::xgb.DMatrix(data = new_tr,label = labels)
    dtest <- xgboost::xgb.DMatrix(data = new_ts,label=ts_label)

    params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.3, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1)
    xgbcv <- xgboost::xgb.cv( params = params, data = dtrain, nrounds = 100, nfold = 5, showsd = T, stratified = T, print.every.n = 10, early.stop.round = 100, maximize = F, metrics = "auc")
    
    xgb1 <- xgboost::xgb.train( params = params, data = dtrain, watchlist = list(train = dtrain, eval = dtest), nrounds = 100, eval_metric = "auc")

    xgbpred <- stats::predict(xgb1,dtest)

    if (model_assessment){
      xgbpred <- ifelse(xgbpred > 0.5,1,0)
      mat <- xgboost::xgb.importance(model = xgb1)
      nb_max <- ifelse(nrow(mat) < 35, nrow(mat), 35)
    } else {
    }
    return(list(model = xgb1, importance_matrix = mat[1:nb_max], ts_label = ts_label,
    	xgbpred = xgbpred, test_matchLogPval = test$matchLogPval))
  } else {
    return(list(training.dataset = train, validation.dataset = test))
  }
}

