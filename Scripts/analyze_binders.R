args = commandArgs( trailingOnly=T )

binders = dir(args[1], pattern="*.recint.score.csv", full.names=T)
nobinders = dir(args[2], pattern="*.recint.score.csv", full.names=T)

load_normalized <- function(file) {
  x = read.csv(file, header = F, row.names=1)
  x$V2[x$V4 == 1]
}

load_nonnormalized <- function(file) {
  x = read.csv(file, header = F, row.names=1)
  x$V3[x$V4 == 1]
}

binders.norm <- t(sapply(binders, load_normalized))
nobinders.norm <- t(sapply(nobinders, load_normalized))

rownames(binders.norm) = gsub(".*-", "", gsub(".recint.score.csv$", "", rownames(binders.norm)))
rownames(nobinders.norm) = gsub(".*-", "", gsub(".recint.score.csv$", "", rownames(nobinders.norm)))

print(sprintf("Loaded %d binders with %d interface residues.", nrow(binders.norm), ncol(binders.norm)))
print(sprintf("Loaded %d no_binders with %d interface residues.", nrow(nobinders.norm), ncol(nobinders.norm)))

D = rbind(binders.norm, nobinders.norm)
Y = factor(rep(c("bind","nobind"), c(nrow(binders.norm), nrow(nobinders.norm) ) ))

Dnonnorm = rbind(
  t(sapply(binders, load_nonnormalized)),
  t(sapply(nobinders, load_nonnormalized))
)

## loading z-scores
zbinders   = read.csv(sprintf("%s/top_zscore.txt", args[1]), sep="\t", as.is=TRUE, row.names="Ligand")
znobinders = read.csv(sprintf("%s/top_zscore.txt", args[2]), sep="\t", as.is=TRUE, row.names="Ligand")

## order the z-scores according to binders.norm
zbinders   = zbinders[ rownames(binders.norm), , drop=FALSE ]
znobinders = znobinders[ rownames(nobinders.norm), , drop=FALSE ]
Z = rbind.data.frame(zbinders, znobinders)

if (any(is.na(zbinders)) || any(is.na(znobinders))) {
  stop("zbinders or znobinders have missing values.")
}

## extraTrees
library(extraTrees)

et = extraTrees(D, Y, numRandomCuts = 2, nodesize = 1, mtry=10, subsetSizes = 30)
yhat = predict(et, D, probability = T)

## crossvalidation
library(optunity)
library(e1071)
folds = c(
          generate_folds(num_instances = nrow(binders.norm),   num_folds = 5), 
          generate_folds(num_instances = nrow(nobinders.norm), num_folds = 5))

for (i in 1:max(folds)) {
  itrain = folds != i
  itest  = folds == i
  
  et = extraTrees(D[itrain,], Y[itrain], numRandomCuts = 2, nodesize = 1, mtry=10)
  yhat = predict(et, D[itest,], probability = T)
  
  m2 <- svm(x = D[itrain,], y=Y[itrain], kernel = "radial", probability = T)
  svm.yhat <- attributes(predict(m2, D[itest,], decision.values = T))$decision.values
  
  ## z-score predictor
  zscore.yhat <- Z$Score[itest]
  
  print(data.frame(
    et.yhat  = yhat, 
    svm.yhat = unname(svm.yhat),
    zscore.yhat = zscore.yhat,
    ytrue    = Y[itest]))
  print(sprintf("ROC(et):     %1.2f", auc_roc(ytrue = 2-as.numeric(Y[itest]), yhat[,"bind"]) ))
  print(sprintf("ROC(svm):    %1.2f", auc_roc(ytrue = 2-as.numeric(Y[itest]), svm.yhat) ))
  print(sprintf("ROC(zscore): %1.2f", auc_roc(ytrue = 2-as.numeric(Y[itest]), -zscore.yhat) ))
}

## svm
m <- svm(x = D, y=as.numeric(Y), kernel = "radial", probability = T, cost = 10, gamma=0.01)
svm.yhat <- predict(m, D, decision.values = T)
svm.yhat

m2 <- svm(x = Dnonnorm, y=as.numeric(Y), kernel = "radial", probability = T)
svm.yhat <- predict(m, Dnonnorm, decision.values = T)
svm.yhat
