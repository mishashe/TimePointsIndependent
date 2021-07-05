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


SamplesInfo[,"status"] <- ifelse(SamplesInfo[,"case_control"]=="Case",1,ifelse(SamplesInfo[,"case_control"]=="Control",1,-1))
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


x <- as.matrix(read.table("/DATA/share/dcis_recurrence/2020-06-11-rnaseq-counts-nki2/gene-expression-nki2.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/DATA/share/dcis_recurrence/2020-06-11-rnaseq-counts-nki2/excluded_samples_nki2.tsv"))
xNKI2 <- x[,!(colnames(x) %in% Excluded$precision_sample)]
x <- as.matrix(read.table("/DATA/share/dcis_recurrence/2020-04-23-rnaseq-counts-kcl2/gene-expression-kcl2.tsv", header = TRUE, row.names=1, sep="\t"))
Excluded <- data.frame(read_tsv("/DATA/share/dcis_recurrence/2020-04-23-rnaseq-counts-kcl2/excluded_samples_kcl2.tsv"))
xKCL2 <- x[,!(colnames(x) %in% Excluded$precision_sample)]
Genes <- intersect(rownames(xKCL2),rownames(xNKI2))
x2 <- cbind(xNKI2[Genes,],xKCL2[Genes,]);
x2 <- x2[,(colnames(x2) %in% rownames(SamplesInfo))]; 
SamplesInfo2 <- SamplesInfo[colnames(x2),]
SamplesInfo2$Set <- "Set2"
x2 <- x2[which(rowMeans(x2)>1),]
xx <- DGEList(x2)
xx <- calcNormFactors(xx,method="TMM")
xx <- cpm(xx, log=TRUE)
x2 <- x2[which(rowSds(xx)>0.75),]

############################################## combine datasets ############################################## 
Genes <- intersect(rownames(xKCL1),rownames(xNKI1))
Genes <- intersect(Genes,rownames(x2))
x <- cbind(xKCL1[Genes,],xNKI1[Genes,],x2[Genes,])
SamplesInfo <- rbind(SamplesInfoKCL1,SamplesInfoNKI1,SamplesInfo2)

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
Ins <- c("KCL1","NKI1","Set2")
p <- NULL
k <- 1
for (i in 1:3)
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
  
  for (j in 1:3)
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

Ins <- c("KCL1","NKI1","Set2")
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




































xNKI1 <- DGEList(xNKI1)
xNKI1 <- calcNormFactors(xNKI1,method="TMM")
cutoff <- 5
drop <- which(apply(cpm(xNKI1), 1, mean) < cutoff)
xNKI1 <- xNKI1[-drop,] 
dim(xNKI1)
Out <- as.matrix(data.frame(time=FU,status=(case_control=="Case")+0))
mm <- model.matrix(~0 + case_control)
pdf(paste0(OutDir,"plots/NKI1.pdf"))
y <- voom(xNKI1, mm, plot = TRUE)
dev.off()

FU_tot <- c(FU_tot,FU)
case_control_tot <- c(case_control_tot,case_control)


############################################## Set 2 ###########################################################################################



samples <- colnames(x2)
patients <- sapply(samples,function(s){SamplesInfo$site_accession[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
FU <- sapply(patients,function(p){FUtable$fu_months[which(FUtable$site_accession==p)][1]})
tissue_pathology <- sapply(samples,function(s){SamplesInfo$tissue_pathology[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
event_type <- sapply(samples,function(s){SamplesInfo$event_type[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_laterality <-sapply(samples,function(s){SamplesInfo$case_laterality[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
radiotherapy <- sapply(samples,function(s){SamplesInfo$radiotherapy[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_control <- sapply(samples,function(s){SamplesInfo$case_control[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_pathology <- sapply(samples,function(s){SamplesInfo$case_pathology[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
radiotherapy[which(is.na(radiotherapy))] <- 999
Ind <- which(radiotherapy %in% c(0) & !is.na(FU)  & event_type=="PRI" & case_control %in% c("Case","Control") & (case_control %in% c("Control") | case_pathology=="Invasive Carcinoma"))
x2 <- x2[,Ind]
samples <- colnames(x2)
patients <- sapply(samples,function(s){SamplesInfo$site_accession[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
FU <- sapply(patients,function(p){FUtable$fu_months[which(FUtable$site_accession==p)][1]})
Age <- sapply(patients,function(p){FUtable$age_diagnose[which(FUtable$site_accession==p)][1]})
tissue_pathology <- sapply(samples,function(s){SamplesInfo$tissue_pathology[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
event_type <- sapply(samples,function(s){SamplesInfo$event_type[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_laterality <-sapply(samples,function(s){SamplesInfo$case_laterality[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
radiotherapy <- sapply(samples,function(s){SamplesInfo$radiotherapy[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_control <- sapply(samples,function(s){SamplesInfo$case_control[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})
case_pathology <- sapply(samples,function(s){SamplesInfo$case_pathology[which(SamplesInfo$sample_id==s & SamplesInfo$experimental_technique=="RNA Sequencing" & SamplesInfo$material_type=="RNA")][1]})


GenesEnsembl <- sapply(rownames(x2),function(g){strsplit(g,split="[.]")[[1]][1]})
GenesNames <- ensembldb::select(EnsDb.Hsapiens.v86, key=GenesEnsembl,columns=c("SYMBOL"),keytype="GENEID")
rownames(GenesNames) <- GenesNames[,1]
GenesNames <- GenesNames[GenesEnsembl,2]
Ind <- which(!is.na(GenesNames))
x2 <- x2[Ind,]
rownames(x2) <- GenesNames[Ind]

x2 <- DGEList(x2)
x2 <- calcNormFactors(x2,method="TMM")
cutoff <- 5
drop <- which(apply(cpm(x2), 1, mean) < cutoff)
x2 <- x2[-drop,] 
dim(x2)
Out <- as.matrix(data.frame(time=FU,status=(case_control=="Case")+0))
mm <- model.matrix(~0 + case_control)
pdf(paste0(OutDir,"plots/Set2.pdf"))
y <- voom(x2, mm, plot = TRUE)
dev.off()

FU_tot <- c(FU_tot,FU)
case_control_tot <- c(case_control_tot,case_control)





library(sva)
Genes <- intersect(rownames(xKCL1),rownames(xNKI1))
Genes <- intersect(Genes,rownames(x2))

Institute <- c(rep("KCL1",ncol(xKCL1)),rep("NKI1",ncol(xNKI1)),rep("Set2",ncol(x2)))
# x0 <- as.matrix(cbind(xKCL1[Genes,],xNKI1[Genes,]))
x0 <- ComBat_seq(as.matrix(cbind(xKCL1[Genes,],xNKI1[Genes,],x2[Genes,])), batch=Institute, group=NULL)

FU <- FU_tot
case_control <- case_control_tot

x0 <- DGEList(x0)
x0 <- calcNormFactors(x0,method="TMM")
x0 <- t(cpm(x0,log=TRUE))

FollowUp <- FU
y0 <- (case_control=="Case") + 0

registerDoParallel(cores = 20)
detach("package:TimePointsPenalized", unload=TRUE)
remove.packages("TimePointsPenalized")

library(rlist)
library(Rcpp)
library(devtools)
library(glmnet)
library(roxygen2)
setwd("~/Documents/Development/TimePointsPenalized")
usethis::use_rcpp()
Rcpp::compileAttributes(pkgdir = ".", verbose = TRUE)
usethis::use_github_actions()
roxygen2::roxygenise(load_code = "source")
pkgbuild::compile_dll()
devtools::document()
system("git add --all .")
system("git commit -m 'debug' ")
system("git push")


library(rlist)
library(Rcpp)
library(devtools)
library(glmnet)
detach("package:TimePointsPenalized", unload=TRUE)
remove.packages("TimePointsPenalized")
install_github("mishashe/TimePointsPenalized", force=TRUE)
library(TimePointsPenalized)
nSamples <- 100
nGenes <- 300
tV <- seq(2,8,1)*12
x0 <- matrix(rnorm(nSamples*nGenes),nrow=nSamples,ncol=nGenes)
beta0 <- rep(0,nGenes)
beta0[1:5] <- 1
beta0[6:10] <- -1
y0 <- round(1/(1+exp(-x0 %*% beta0)))
beta <- rep(0,nGenes*length(tV))
tV <- seq(4,7,1)*12
FollowUp <- 2*12 + (8-2)*12*runif(nSamples)
lam1V <- 10^seq(-1.0,-4.5,-0.025)
gamma <- 10
rownames(x0) <- paste0("S",1:nSamples)
colnames(x0) <- paste0("G",1:nGenes)
fits <- fitTimePointsPenalized(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0))
fits0 <- fitTimePointsNonPenalized(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0))




###############################################################################################################################
RcppArmadillo.package.skeleton(name = "TRY", list = character(), 
                               environment = .GlobalEnv, path = "/home/misha/Documents/Development/", force = TRUE, 
                               code_files = character(), example_code = TRUE)

require("devtools")



setwd("~/Documents/Development/")
system("R CMD build TimePointsPenalized")
system("R CMD ldd /home/misha/Documents/Development/TimePointsPenalized.Rcheck/00LOCK-TimePointsPenalized/00new/TimePointsPenalized/libs/TimePointsPenalized.so")
# install.packages("TimePointsPenalized_1.0.tar.gz", repos = NULL)
system("R CMD check TimePointsPenalized")
system("R CMD INSTALL TimePointsPenalized_1.0.tar.gz")
# git add .
# git commit -a
# git push

detach("package:TimePointsPenalized", unload=TRUE)
install_github("mishashe/TimePointsPenalized")
library(TimePointsPenalized)
tV <- seq(2,8,1)*12
lam1V <- 10^seq(-1.0,-2.5,-0.025)
gamma <- 10
fitTimePointsPenalized(y0, x0, FollowUp, lam1V, gamma, tV, standardize=TRUE, Clinical0=data.frame(case_control0=y0), cores=1)









for (i in 1:ncol(x0))
{
  x0[,i] <- (x0[,i] - mean(x0[,i]))/sd(x0[,i])
}














x0 <- cbind(matrix(1,ncol=1,nrow=nrow(x0)),x0)







colnames(x0)[1] <- "Constant"

Institute0 <- Institute


Institute <- c()
tV <- seq(2,8,1)*12
y <- c()
samples <- c()
times <- c()
Set <- c()
for (it in 1:length(tV))
{
  t <- tV[it]
  case_controlT <- ifelse(FU>t,"Control",ifelse(case_control=="Case","Case","Exclude"))
  y <- c(y,case_controlT)
  samples <- c(samples,rownames(x0))
  times <- c(times,rep(t,nrow(x0)))
  Institute <- c(Institute,Institute0)
}

# x <- kronecker(diag(length(tV)),x0)
y <- matrix(y,nrow=nrow(x0)*length(tV),ncol=1)
y <- ifelse(y=="Case",1,ifelse(y=="Control",0,-1))
# colnames(x) <- 1:ncol(x)
# rownames(x) <- 1:nrow(x)
rownames(y) <- 1:(nrow(x0)*length(tV))
samplesT <- 1:(nrow(x0)*length(tV))
GenesT <- rep("",ncol(x0)*length(tV))
for (it in 1:length(tV))
{
  # colnames(x)[(1+(it-1)*ncol(x0)):(it*ncol(x0))] <- paste0(colnames(x0),"_t_",it)
  # rownames(x)[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
  rownames(y)[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
  GenesT[(1+(it-1)*ncol(x0)):(it*ncol(x0))] <- paste0(colnames(x0),"_t_",it)
  samplesT[(1+(it-1)*nrow(x0)):(it*nrow(x0))] <- paste0(rownames(x0),"_t_",it)
}

Ind <- which(y %in% c(0,1))
# x <- x[Ind,]
y <- y[Ind,,drop=TRUE]
samples <- samples[Ind]
samplesT <- samplesT[Ind]
times <- times[Ind]
Institute <- Institute[Ind]
names(times) <- samplesT
names(samples) <- samplesT
names(Institute) <- samplesT

fold <- rep(0,length(samples))
names(fold) <- samples
samples_U <- unique(samples)
for (is in 1:length(samples))
{
  fold[is] <- which(samples[is]==samples_U)
}

w <- times*0
for (Ins in unique(Institute))
{
  IndI <- which(Institute==Ins)
  for (it in 1:length(tV[IndI]))
  {
    Ind <- which(times[IndI]==tV[it] & y[IndI] %in% c(0,1))
    w[IndI][Ind][which(y[IndI][Ind]==0)] <- length(Ind)/(0.1+sum(y[Ind]==0))
    w[IndI][Ind][which(y[IndI][Ind]==1)] <- length(Ind)/(0.1+sum(y[Ind]==1))
    w[IndI][Ind] <- w[IndI][Ind]/sum(w[IndI][Ind])/length(tV)/length(Ind)*length(w[IndI])
  }
}
################################################### TRY ######################################################################







































cvlMy <- function(par,y, x0, lam1V)
{
  fits <- foreach(f = c(0,unique(fold)), .inorder=TRUE) %dopar%
    {
      if (f==0) {Ind <- which(Institute=="KCL1")
      }    else {Ind <- which(fold!=f & Institute=="KCL1")}
      beta <- rep(0,ncol(x0)*length(tV))
      betaM <- matrix(0,nrow=length(beta),ncol=length(lam1V))
      rownames(betaM) <- rep(colnames(x0),length(tV))
      IndFor0 <- c()
      IndTFor0 <- c()
      for (it in 1:length(tV))
      {
        IndT <- which(times[Ind]==tV[it])
        w[Ind][IndT][which(y[Ind][IndT]==0)] <- 1/(0.001+sum(y[Ind]==0))
        w[Ind][IndT][which(y[Ind][IndT]==1)] <- 1/(0.001+sum(y[Ind]==1))
        IndFor0 <- c(IndFor0,which(rownames(x0) %in% samples[Ind][IndT]))
        IndTFor0 <- c(IndTFor0,which(rownames(x0) %in% samples[Ind][IndT])*0+it)
      }
      for (ilam1 in 1:length(lam1V))
      {
        lam1 <- lam1V[ilam1]
        lam2 <- par*lam1
        beta <- FitRound(x0, y[Ind], tV, lam1, lam2, beta, w[Ind],IndFor0,IndTFor0);
        betaM[,ilam1] <- beta
      }
      fit <- list(beta=betaM)
      rownames(fit$beta) <- GenesT
      return(fit)
    }
  #gc(verbose = FALSE, reset = TRUE, full = TRUE)
  
  predictions <- foreach(f = unique(fold), .combine=rbind, .inorder=TRUE) %dopar%
    {
      Ind <- which(fold==f & Institute=="KCL1")
      M <- matrix(0,nrow=length(Ind),ncol=length(lam1V),byrow=TRUE)
      for (it in 1:length(tV))
      {
        IndT <- which(times[Ind]==tV[it])
        IndFor0 <- which(rownames(x0) %in% samples[Ind][IndT])
        M[IndT,] <- M[IndT,,drop=FALSE] + x0[IndFor0,,drop=FALSE] %*% fits[[f+1]]$beta[1:ncol(x0) + (it-1)*ncol(x0),,drop=FALSE]
      }
      return(data.frame(pred=1/(1+exp(-M )),status=y[Ind],
                        sample=samplesT[Ind], time=times[samplesT[Ind]], Institute=Institute[samplesT[Ind]]))
    }
  Ind <- which(Institute!="KCL1")
  M <- matrix(0,nrow=length(Ind),ncol=length(lam1V),byrow=TRUE)
  for (it in 1:length(tV))
  {
    IndT <- which(times[Ind]==tV[it])
    IndFor0 <- which(rownames(x0) %in% samples[Ind][IndT])
    M[IndT,] <- M[IndT,,drop=FALSE] + x0[IndFor0,,drop=FALSE] %*% fits[[1]]$beta[1:ncol(x0) + (it-1)*ncol(x0),,drop=FALSE]
  }
  predictions <- rbind(predictions,data.frame(pred=1/(1+exp(-M )),status=y[Ind],
                                              sample=samplesT[Ind], time=times[samplesT[Ind]], Institute=Institute[samplesT[Ind]]))
  AUCmean <- lam1V*0
  pValLog10mean <- lam1V*0
  nonZeroNumGenes <- lam1V*0
  AUC <- tV*0
  pValLog10 <- tV*0
  for(i in 1:length(lam1V))
  {
    for (it in 1:length(tV))
    {
      t <- tV[it]
      Statuses <- predictions$status[which(predictions$time==t & predictions$Institute=="KCL1")]
      Preds <- predictions[,i][which(predictions$time==t& predictions$Institute=="KCL1")]
      pROC_obj <- roc(Statuses,Preds,plot=FALSE,direction="<",quiet=TRUE)
      AUC[it] <- auc(pROC_obj)
      pValLog10[it] <- ifelse(var(Preds)==0,1,-log10(t.test(Preds[Statuses==1],Preds[Statuses==0], alternative = "greater", paired = FALSE)$p.value))
    }
    # FigureMerit[i] <- median(predictions[,i][which(predictions$status==1)]) - median(predictions[,i][which(predictions$status==0)])
    AUCmean[i] <- mean(AUC)
    pValLog10mean[i] <- mean(pValLog10)
    nonZeroNumGenes[i] <- sum(fits[[1]]$beta[,i]!=0)/length(tV)
  }
  pdf(paste0(OutDir,"plots/",par,"_cv_lambda.pdf"),height=10)
  par(mfrow = c(3, 1))
  plot(log10(lam1V),AUCmean)
  plot(log10(lam1V),pValLog10mean)
  plot(log10(lam1V),nonZeroNumGenes)
  lines(log10(lam1V),log10(lam1V)*0)
  dev.off()
  FigureMerit <- AUCmean
  
  iMax <- which.max(FigureMerit)
  FigureMeritMax <- FigureMerit[iMax]
  colnames(predictions)[colnames(predictions)=="pred"] <- "nonpred"
  colnames(predictions)[iMax] <- "pred"
  lam1Max <- lam1V[iMax]  
  Res <- data.frame(pred=predictions$pred,time=times[predictions$sample],sample=predictions$sample,status=predictions$status, nGenes = nonZeroNumGenes[iMax], Institute=predictions$Institute)
  if (nonZeroNumGenes[iMax]>1) {MakeFigures(tV=tV,Res=Res,fit=fits[[1]],iMax=iMax,par=par)}
  print(paste0("par = ",par))
  print(paste0("lam2 = ",lam1Max*par))
  print(paste0("lam1Max = ",lam1Max))
  print(FigureMeritMax)
  gc(TRUE)
  return(list(AUCmean[iMax],nonZeroNumGenes[iMax],lam1V[iMax]))
}














library(expm)
library(glmnet)
library(Rcpp)
library(pROC)
library(stringr)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
Rcpp::sourceCpp("/home/m.sheinman/Development/precision-CaseControl/src/models/TimePointsMy/Functions.cpp")
source("/home/m.sheinman/Development/precision-CaseControl/src/models/TimePointsMy/Functions.R")
registerDoMC(cores=50)






















lam2 <- 0.1
beta <- rep(0,ncol(x))
p <- 1/(1+exp(-x %*% beta))
betaM <- matrix(0,nrow=length(beta),ncol=0)
for (lam1 in 10^seq(log10(0.1),log10(0.002),-0.01))
{
  beta <- FitRound(x, y, tV, lam1, lam2, beta);
  betaM <- cbind(betaM,beta)
  nNonZero <- sum(beta!=0)
  print(lam1)
  print(nNonZero)
}




################################# CROSS VALIDATION ###############################################################################
library(penalized)
library(glmnet)
registerDoMC(cores=40)
lamV <- 10^seq(log10(0.1),log10(0.002),-0.01)
Res <- foreach (s = unique(samples), .inorder=TRUE, .combine=rbind) %dopar%
  {
    print(s)
    preds <- c()
    Ind <- which(!(s == samples))
    xI <- x[Ind,]
    yI <- y[Ind,1,drop=FALSE]
    glmnetGenes <- Calc_glmnetGenes(xI,yI,tV,nGenes=20)
    chra <- c()
    for (g in 1:length(glmnetGenes)){chra <- c(chra,rep(g,length(tV)))}
    xG <- x[,sapply(colnames(x),function(str){strsplit(str,"_")[[1]][1]}) %in% glmnetGenes]
    xG <- xG[,sort(colnames(xG))]
    startbeta <- rep(0,ncol(xG))
    
    xx <- xG[Ind,,drop=FALSE]
    yy <- yI[,1,drop=TRUE]
    
    
    beta <- chra*0
    for (lam1 in lamV)
    {
      for (lam2 in lamV)
      {
        lam1prev <- lam1
        lam2prev <- lam2
        fit <- penalized(yy, xx,lambda1 = lam1,lambda2 = lam2, fusedl=chra, model = c("cox", "logistic", "linear", "poisson")[2],
                         startbeta=beta, epsilon = 1e-5,maxiter=10000, standardize = FALSE, trace = TRUE)
        beta <- fit@penalized
        IndNew <- which((s == samples) ) 
        print(lam1)
        print(lam2)
        print(sum(beta!=0))
        preds <- cbind(preds,predict(fit, xG[IndNew,,drop=FALSE]))
      }
    }
    data.frame(pred=preds,time=times[IndNew],sample=samples[IndNew],status=y[IndNew],Set=Set[IndNew])
  }


library(pROC)
iMax <- 0 
FigureOfMerit <- -Inf
for (i in 1:(length(lamV))^2)
{
  for (it in 1:(length(tV)-1))
  {
    t <- tV[it]
    pROC_obj <- roc(Res[which(Res$time==t),]$status,Res[which(Res$time==t),][,i],plot=FALSE,direction="<")
    AUC[it] <- auc(pROC_obj)
  }
  if (min(AUC)>FigureOfMerit)
  {
    iMax <- i
    FigureOfMerit <- min(AUC)
    print(i)
    print(FigureOfMerit)
  }
}
print(iMax)

glmnetGenes <- Calc_glmnetGenes(x,y,tV,nGenes=20)
chra <- c()
for (g in 1:length(glmnetGenes)){chra <- c(chra,rep(g,length(tV)))}
xG <- x[,sapply(colnames(x),function(str){strsplit(str,"_")[[1]][1]}) %in% glmnetGenes]
xG <- xG[,sort(colnames(xG))]
startbeta <- rep(0,ncol(xG))
fit <- penalized(y, xG,lambda1 = lamV[iMax],lambda2 = lamV[iMax], fusedl=chra, model = c("cox", "logistic", "linear", "poisson")[2],
                 startbeta=startbeta, epsilon = 1e-5,maxiter=10000, standardize = FALSE, trace = TRUE)
beta <- fit@penalized
write.table(beta[beta!=0],
            file = "/home/m.sheinman/Development/precision-CaseControl/src/models/Pathways/plots/TimePoints/grplasso/genes.tsv",
            append=FALSE,row.names=TRUE,col.names=NA,sep = "\t",quote=FALSE)


colnames(Res)[1:(length(lamV)^2)] <- 1:(length(lamV)^2)
colnames(Res)[iMax] <- "pred"
AUC <- tV*0
library(pROC)
pdf("/home/m.sheinman/Development/precision-CaseControl/src/models/Pathways/plots/TimePoints/grplasso/ROC_cv.pdf")
par(mfrow=c(5,4))
for (it in 1:length(tV))
{
  t <- tV[it]
  pROC_obj <- roc(Res[which(Res$time==t),]$status,Res[which(Res$time==t),]$pred,smoothed = FALSE,
                  ci=FALSE, ci.alpha=0.95, stratified=FALSE,
                  plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                  print.auc=TRUE, show.thres=TRUE,cex.lab=1.0, cex.axis=1.0, cex.main=1.0, cex.sub=1.0,direction="<")
  # sens.ci <- ci.se(pROC_obj)
  # plot(sens.ci, type="shape", col="lightblue")
  # plot(sens.ci, type="bars")
  title(t/12)
  AUC[it] <- auc(pROC_obj)
}
dev.off()
pdf("/home/m.sheinman/Development/precision-CaseControl/src/models/Pathways/plots/TimePoints/grplasso/AUC_cv.pdf")
plot(tV/12,AUC)
dev.off()
library(ggpubr)
pdf("/home/m.sheinman/Development/precision-CaseControl/src/models/Pathways/plots/TimePoints/grplasso/box_cv.pdf")
p <- ggboxplot(Res, x = "status", y = "pred",shape="Set",
               color = "Set",add="jitter",add.params = list(size = 1),ylim = c(0, 1)) +  stat_compare_means(method = "wilcox.test")
facet(p, facet.by = "time")
dev.off()






















