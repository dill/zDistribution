### Check m1.2 model performance with AUC
### November 2025


library(pROC)


# check model fit
## Area under ROC (AUC) requires a test set of data
# m1.2preds <- predict(m1.2, m1.2_predictions, type = "response")
# 
# roc_object <- roc(detect_data$Detected, m1.2preds)
# auc(roc_object)