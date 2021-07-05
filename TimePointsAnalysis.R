library(expm)
library(glmnet)
library(Rcpp)
library(pROC)
library(stringr)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
require(doParallel)
library(RcppArmadillo)
library(Rcpp)
library(devtools)
library(usethis)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(readr)
library(ensembl)
library(edgeR)
library(rlist)
library(matrixStats)

cv.glmnetMy <- function(x,y,alpha)
{
  weights <- y
  weights[y==1] <- 1/sum(y==1)/2
  weights[y==0] <- 1/sum(y==0)/2
  fitcv <- list()
  fitcv$fit <- glmnet(x=x, y=y, family="binomial", weights=weights, intercept = TRUE, 
                      standardize = TRUE, alpha=alpha,standardize.response = FALSE)
  fitcv$fit.preval <- matrix(0, nrow=length(y),ncol=length(fitcv$fit$lambda))
  for (i in 1:length(y))
  {
    weightsI <- y[-i]
    weightsI[y[-i]==1] <- 1/sum(y[-i]==1)/2
    weightsI[y[-i]==0] <- 1/sum(y[-i]==0)/2
    fit <- glmnet(x=x[-i,,drop=FALSE], y=y[-i], family="binomial", weights=weightsI, lambda=fitcv$fit$lambda, 
                  intercept = TRUE, standardize = TRUE, alpha=alpha,standardize.response = FALSE)
    fitcv$fit.preval[i,] <- predict(fit, newx=x[i,,drop=FALSE],type="response",s = fitcv$fit$lambda,alpha=alpha)
  }
  fitcv$auc <- sapply(1:ncol(fitcv$fit.preval),function(i){auc(roc(y,fitcv$fit.preval[,i],plot=FALSE,direction="<",quiet=TRUE))})
  fitcv$loglikelihood <- sapply(1:ncol(fitcv$fit.preval),function(i){
    exp(sum(weights*(y*log(fitcv$fit.preval[,i])+(1-y)*log(1-fitcv$fit.preval[,i]))))
    })
  fitcv$lambda <- fitcv$fit$lambda
  fitcv$lambda.min <- fitcv$lambda[which.max(fitcv$loglikelihood)]
  return(fitcv)
}
system("export OPENBLAS_NUM_THREADS=1")
system("export GOTO_NUM_THREADS=1")
system("export OMP_NUM_THREADS=1")
setwd("/home/m.sheinman/Development/precision-CaseControl/src/models/TimePoints")

############################################## import clinical data ###########################################################################################

SamplesInfo <- as.data.frame(read_tsv("/DATA/share/dcis_recurrence/2021-05-31-sample-annotation/sample_annotation.tsv", col_names = TRUE))
SamplesInfo <- SamplesInfo[SamplesInfo$material_type=="RNA",]
SamplesInfo <- SamplesInfo[SamplesInfo$event_type=="PRI",]
SamplesInfo <- SamplesInfo[SamplesInfo$experimental_technique=="RNA Sequencing",]
SamplesInfo <- SamplesInfo[which(SamplesInfo$radiotherapy=="No"),]
SamplesInfo <- SamplesInfo[which(SamplesInfo$laterality=="Ipsilateral"),]
SamplesInfo <- SamplesInfo[which(SamplesInfo$case_pathology=="Invasive Carcinoma" | is.na(SamplesInfo$case_pathology)),]
SamplesInfo <- SamplesInfo[!is.na(SamplesInfo$inv_time),]


SamplesInfo$experimental_technique_column <- 0
SamplesInfo <- unique(SamplesInfo)
rownames(SamplesInfo) <- SamplesInfo$sample_id


SamplesInfo[,"status"] <- ifelse(SamplesInfo[,"case_control"]=="Case",1,ifelse(SamplesInfo[,"case_control"]=="Control",0,-1))
SamplesInfo[,"time"] <- as.numeric(SamplesInfo[,"inv_time"])/365
# SamplesInfo[,"time"] <- as.numeric(SamplesInfo[,"follow_up"])/365

############################################## import expression levels ###########################################################################################

xKCL1 <- as.matrix(read.table("/DATA/share/dcis_recurrence/2021-04-14-rnaseq-kcl/gene-expression-kcl.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/DATA/share/dcis_recurrence/2021-04-14-rnaseq-kcl/excluded_samples_kcl.tsv"))$precision_sample
xKCL1 <- xKCL1[,!(colnames(xKCL1) %in% Excluded)]; 
xKCL1 <- xKCL1[,(colnames(xKCL1) %in% rownames(SamplesInfo))]; 
SamplesInfoKCL1 <- SamplesInfo[colnames(xKCL1),]
SamplesInfoKCL1$Set <- "KCL1"
xKCL1 <- xKCL1[which(rowMeans(xKCL1)>1),]
xx <- DGEList(xKCL1)
xx <- calcNormFactors(xx,method="TMM")
xx <- cpm(xx, log=TRUE)
xKCL1 <- xKCL1[which(rowSds(xx)>0.75),]

xNKI1 <- as.matrix(read.table("/DATA/share/dcis_recurrence/2019-07-01-rnaseq-counts/gene-expression-nki.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/home/m.sheinman/Development/precision-CaseControl/data/processed/excluded_samples_nki.tsv"))
xNKI1 <- xNKI1[,!(colnames(xNKI1) %in% Excluded$SampleID)]; 
xNKI1 <- xNKI1[,(colnames(xNKI1) %in% rownames(SamplesInfo))]; 
SamplesInfoNKI1 <- SamplesInfo[colnames(xNKI1),]
SamplesInfoNKI1$Set <- "NKI1"
xNKI1 <- xNKI1[which(rowMeans(xNKI1)>1),]
xx <- DGEList(xNKI1)
xx <- calcNormFactors(xx,method="TMM")
xx <- cpm(xx, log=TRUE)
xNKI1 <- xNKI1[which(rowSds(xx)>0.75),]

xNKI2 <- as.matrix(read.table("/DATA/share/dcis_recurrence/2020-06-11-rnaseq-counts-nki2/gene-expression-nki2.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/DATA/share/dcis_recurrence/2020-06-11-rnaseq-counts-nki2/excluded_samples_nki2.tsv"))
xNKI2 <- xNKI2[,!(colnames(xNKI2) %in% Excluded$precision_sample)]
xNKI2 <- xNKI2[,(colnames(xNKI2) %in% rownames(SamplesInfo))]; 
SamplesInfoNKI2 <- SamplesInfo[colnames(xNKI2),]
SamplesInfoNKI2$Set <- "NKI2"
xNKI2 <- xNKI2[which(rowMeans(xNKI2)>1),]
xx <- DGEList(xNKI2)
xx <- calcNormFactors(xx,method="TMM")
xx <- cpm(xx, log=TRUE)
xNKI2 <- xNKI2[which(rowSds(xx)>0.75),]

xKCL2 <- as.matrix(read.table("/DATA/share/dcis_recurrence/2020-04-23-rnaseq-counts-kcl2/gene-expression-kcl2.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/DATA/share/dcis_recurrence/2020-04-23-rnaseq-counts-kcl2/excluded_samples_kcl2.tsv"))
xKCL2 <- xKCL2[,!(colnames(xKCL2) %in% Excluded$precision_sample)]
xKCL2 <- xKCL2[,(colnames(xKCL2) %in% rownames(SamplesInfo))]; 
SamplesInfoKCL2 <- SamplesInfo[colnames(xKCL2),]
SamplesInfoKCL2$Set <- "KCL2"
xKCL2 <- xKCL2[which(rowMeans(xKCL2)>1),]
xx <- DGEList(xKCL2)
xx <- calcNormFactors(xx,method="TMM")
xx <- cpm(xx, log=TRUE)
xKCL2 <- xKCL2[which(rowSds(xx)>0.75),]

############################################## combine datasets ############################################## 
Genes <- intersect(rownames(xKCL1),rownames(xNKI1))
Genes <- intersect(Genes,rownames(xNKI2))
Genes <- intersect(Genes,rownames(xKCL2))
x <- cbind(xKCL1[Genes,],xNKI1[Genes,],xKCL2[Genes,],xNKI2[Genes,])
SamplesInfo <- rbind(SamplesInfoKCL1,SamplesInfoNKI1,SamplesInfoKCL2,SamplesInfoNKI2)

############################################## plots of follow-up times and time-to-recurrence ############################################## 
pdf(paste0("./plots/TimesHist.pdf"))
p <- ggboxplot(SamplesInfo, x = "Set", y = "time",
               color = "case_control",add="jitter",add.params = list(size = 1))
# facet(p, facet.by = "time")
p
dev.off()

############################################## change genes names ############################################## 
GenesEnsembl <- sapply(rownames(x),function(g){strsplit(g,split="[.]")[[1]][1]})
GenesNames <- ensembldb::select(EnsDb.Hsapiens.v86, key=GenesEnsembl,columns=c("SYMBOL"),keytype="GENEID")
rownames(GenesNames) <- GenesNames[,1]
GenesNames <- GenesNames[GenesEnsembl,2]
Ind <- which(!is.na(GenesNames))
x <- x[Ind,]
rownames(x) <- GenesNames[Ind]

############################################## COMBAT ############################################## 
library(sva)
print(dim(x))
x0 <- ComBat_seq(x, batch=SamplesInfo$Set, group=NULL)

############################################## voom ############################################## 
# pdf(paste0("./plots/x.pdf"))
# x0 <- t(voom(x0, model.matrix(~ 0 + SamplesInfo$Set), plot = TRUE)$E)
# dev.off()

x0 <- DGEList(x0)
x0 <- calcNormFactors(x0,method="TMM")
x0 <- t(cpm(x0, log=TRUE))

# for (i in 1:ncol(x0))
# {
#   x0[,i] <- (x0[,i] - mean(x0[,i]))/sd(x0[,i])
# }



############################################# Cox analysis ##############################################
library(doParallel)
library(survival)
library(survminer)
registerDoParallel(40)
Ins <- c("KCL1","NKI1","KCL2","NKI2")
p <- NULL
k <- 1
for (i in 1:length(Ins))
{
  print(i)
  I <- Ins[i]
  IndI <- which(SamplesInfo$Set==Ins[i])
  fitcv <- cv.glmnet(as.matrix(x0[IndI,]), as.matrix(SamplesInfo[IndI,c("status","time")]), family = "cox",
                     alpha = 0.5, 
                     nfolds = 10, 
                     parallel = TRUE,
                     type.measure = "C", relax=FALSE, keep=TRUE,grouped=FALSE)
  pdf(paste0("./plots/",I,"_cox_cv.pdf"))
  plot(fitcv)
  dev.off()
  
  for (j in 1:length(Ins))
  {
    J <- Ins[j]
    IndJ <- which(SamplesInfo$Set==Ins[j])
    if (i!=j)
    {
      predsJ <- predict(fitcv, newx=x0[IndJ,], type="response",s = fitcv$lambda.min,alpha=alpha)[,1]
    }
    else
    {
      predsJ <- fitcv$fit.preval[,which(fitcv$lambda==fitcv$lambda.min)]
    }
    dat <- data.frame(SamplesInfo[IndJ,], HR=(predsJ))
    dat$y <- ifelse(dat$HR < median(dat$HR),"low risk","high risk")
    fitKM <- survfit(Surv(time, status) ~ y,data = dat)
    p[[k]] <- ggsurvplot(fitKM, conf.int = FALSE, surv.median.line = c('hv'), data = dat, 
                         pval = TRUE, pval.method = TRUE,
                         risk.table = FALSE, title=paste0(I," vs. ",J),
                         test.for.trend = length(table(dat$y))>2)
    k <- k+1
  }
}
pdf(paste0("./plots/CoxTable.pdf"))
arrange_ggsurvplots(p, print = TRUE,ncol = 2, nrow = 2, risk.table.height = 0.4)
dev.off()

############################################## Independent time points analysis ############################################## 
Ins <- c("KCL1","NKI1","KCL2","NKI2")
tV <- seq(round(min(SamplesInfo[,"time"]),0)+0.5, round(max(SamplesInfo[,"time"]),0)-0.5,0.25)
library(doParallel)
registerDoParallel(40)
alpha <- 0.5

Results <- foreach (it = 1:length(tV), .inorder = TRUE) %dopar%{
  ResultsIJ <- list(list())
  y <- ifelse(SamplesInfo$time>tV[it],0,ifelse(SamplesInfo$status==1,1,-1))
  for (i in (1:length(Ins)))
  {
    print(i)
    IndI <- which((SamplesInfo$Set==Ins[i]) & (y %in% c(0,1)))
    if (sum(y[IndI]==1)<4 | sum(y[IndI]==0)<4){
      ResultsIJ[[Ins[i]]]$fit_exists <- FALSE
    } else{
      ResultsIJ[[Ins[i]]]$fit_exists <- TRUE
      fitcv <- cv.glmnetMy(x=x0[IndI,,drop=FALSE], y=y[IndI], alpha=alpha)
      pdf(paste0("./plots/cvs/cv_",Ins[i],"_t",tV[it],".pdf"))
      # plot(fitcv, main=paste0("Set ",Ins[i]," Time point = ",tV[it]))
      plot(log10(fitcv$lambda), fitcv$auc,main=paste0("Set ",Ins[i]," Time point = ",tV[it]))
      plot(log10(fitcv$lambda), fitcv$loglikelihood,main=paste0("Set ",Ins[i]," Time point = ",tV[it]))
      dev.off()
      pdf(paste0("./plots/cvs/box_",Ins[i],"_t",tV[it],".pdf"))
      p <- ggboxplot(data.frame(status=y[IndI],preds=fitcv$fit.preval[,which.max(fitcv$loglikelihood)]), x = "status", y = "preds",
                     color = "status",add="jitter",add.params = list(size = 1)
      ) +  stat_compare_means(method = "wilcox.test",method.args = list(alternative = "greater"),vjust=1.2) #+ theme(text = element_text(size = 10))
      print(p)
      dev.off()
      ResultsIJ[[Ins[i]]]$preds <- SamplesInfo
      ResultsIJ[[Ins[i]]]$preds[,paste0("status_T_",it)] <- y
      ResultsIJ[[Ins[i]]]$preds[IndI,paste0(Ins[i],"_pred_T_",it)] <- fitcv$fit.preval[,which.max(fitcv$auc)]
      ResultsIJ[[Ins[i]]]$preds[-IndI,paste0(Ins[i],"_pred_T_",it)] <- predict(fitcv$fit, newx=x0[-IndI,], relax = FALSE,
                                                                               type="response",
                                                                               s = fitcv$lambda.min,alpha=alpha)[,1]
    }
  }
  ResultsIJ
}

ResultsMerge <- SamplesInfo
for (it in 1:length(tV)){
  for (i in (1:length(Ins))){
    if (Results[[it]][[Ins[i]]]$fit_exists){
      ResultsMerge <- merge(ResultsMerge,Results[[it]][[Ins[i]]]$preds)
    }
  }
}
write.table(ResultsMerge,
            file = paste0("./plots/PredsTable.csv"),
            append=FALSE,row.names=FALSE,col.names=TRUE,sep = "\t",quote=FALSE)

############################################## assessments of results ############################################## 
Results <- as.data.frame(read_tsv(paste0("./plots/PredsTable.csv"), col_names = TRUE))
for (i in (1:length(Ins))){
  for (j in (1:length(Ins))){
    TimePoints <- c()
    aucTimePoints <- c()
    pValTimePoints <- c()
    dataTimePoints <- data.frame()
    for (it in 1:length(tV)){
      status <- Results[Results$Set==Ins[j],paste0("status_T_",it)]
      preds <- Results[Results$Set==Ins[j],paste0(Ins[i],"_pred_T_",it)]
      if (!is.null(preds)) {dataTimePoints <- rbind(dataTimePoints,data.frame(status=status,preds=preds,TimePoint=tV[it]))}
      if (!is.null(preds) & sum(status==0)>3 & sum(status==1)>3 ){
        pROC_obj <- roc(status[status %in% c(0,1)],preds[status %in% c(0,1)],plot=FALSE,direction="<",quiet=TRUE)
        aucTimePoints <- c(aucTimePoints,auc(pROC_obj))
        pValTimePoints <- c(pValTimePoints,wilcox.test(preds[status %in% c(1)],preds[status %in% c(0)], alternative = "greater", paired = FALSE)$p.value)
        TimePoints <- c(TimePoints, tV[it])
      }
    }
    pdf(paste0("./plots/AUC_",Ins[i],"_",Ins[j],".pdf"))
    plot(TimePoints,TimePoints/TimePoints*0.5, type = "l",ylim=c(0,1))
    points(TimePoints[pValTimePoints<0.05],aucTimePoints[pValTimePoints<0.05], pch = 19,col="black")
    points(TimePoints[pValTimePoints>=0.05],aucTimePoints[pValTimePoints>=0.05], pch = 19,col="grey")
    dev.off()
    pdf(paste0("./plots/Box_",Ins[i],"_",Ins[j],".pdf"),width=15,height=15)
    p <- ggboxplot(dataTimePoints, x = "status", y = "preds",
                   color = "status",add="jitter",add.params = list(size = 1)
                   # ,ylim = c(0, 1)
    ) +  stat_compare_means(method = "wilcox.test",comparisons = list( c("0", "1"),c("-1", "1"),c("-1", "0") ),method.args = list(alternative = "less"),vjust=1.2) #+ theme(text = element_text(size = 10))
    p <- facet(p, facet.by = "TimePoint")
    print(p)
    dev.off()
  }
}



























