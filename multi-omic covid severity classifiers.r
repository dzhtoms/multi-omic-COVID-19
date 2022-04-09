library(caret)
library(tidyverse) 
library(readr)
library(mice) 
library(viridis)
library(doParallel)
library(MLmetrics)

df.annotation = read.csv('annotation_casecontrol.csv')
df.annotation$samples = paste0('sample.',df.annotation$participant)
rownames(df.annotation) = df.annotation$samples

#metabolomics data
dat0 = read.csv('sample.uniqmetabolite.manual.annotation.csv')
rownames(dat0) = dat0$Description
dat.sample = dat0[,-c(1:19)]
dat.sample[dat.sample==0]=NA
df.met = log2(dat.sample)
#proteomics data
dat.tempt = read.csv('proteinmatrix.csv');dim(dat.tempt)
rownames(dat.tempt)=dat.tempt$X
df.prot = dat.tempt[,-1];dim(df.prot)
protein.run.df = read.csv('integratation/run_batch20210427.csv'); dim(protein.run.df)
#annotation[,1]<-gsub("-",".",annotation[,1])
rownames(protein.run.df)<-protein.run.df[,1]
protein.run.df$samples = paste0('sample.',protein.run.df$participant0); head(protein.run.df)
colnames(df.prot)=protein.run.df[colnames(df.prot),'samples'];df.prot[1:3,1:3];dim(df.prot)

min.df.prot = min(apply(df.prot,1,function(x){min(na.omit(x))})); min.df.prot
min.df.met = min(apply(df.met,1,function(x){min(na.omit(x))}));min.df.met
#imputation
df.prot.imput = as.matrix(Hmisc::impute(df.prot,(min.df.prot-1)))
df.met.imput = as.matrix(Hmisc::impute(df.met,(min.df.met-1)))
#
df.matrix0 = rbind(df.prot.imput,df.met.imput); dim(df.matrix0)
df.matrix1 = t(df.matrix0);dim(df.matrix1);df.matrix2 =as.data.frame(df.matrix1)
df.matrix2$group = annot_filter$group2
df.matrix2$age = df.annotation[rownames(df.matrix2),'age']
df.matrix2$sex = df.annotation[rownames(df.matrix2),'sex']
df.matrix3 = df.matrix2[,c(1500:1502,1:1499)]
df.matrix4 = fastDummies::dummy_cols(df.matrix3, select_columns ='sex')
df.matrix5 = df.matrix4[,-3]
df.matrix5$group = factor(df.matrix5$group,levels = c('H','M','S'))

# split the dataset into training and validation=======================================================================
set.seed(42)
index <- createDataPartition(df.matrix5$group, p = 0.7, list = FALSE)
train_data <- df.matrix5[index,];dim(train_data)
test_data  <- df.matrix5[-index,];dim(test_data)
#Handling class imbalance with weighted methods
weight_groupclass = ifelse(train_data$group=='H',sum(table(train_data$group))/(3*table(train_data$group)[1]),
                           ifelse(train_data$group=='M',sum(table(train_data$group))/(3*table(train_data$group)[2]),
                                  sum(table(train_data$group))/(3*table(train_data$group)[3])))
#normalization
mean = apply(train_data[,-c(1,1501:1502)],2,function(x) mean(x)); length(mean)
sd = apply(train_data[,-c(1,1501:1502)],2,function(x) sd(x)); length(sd)
scaler = data.frame(mean = mean, sd = sd); rownames(scaler) = names(train_data)[-c(1,1501:1502)]
#normalization function
norm_fun = function(data,scalerr){
  data.raw = data; dataset.norm = data; 
  scal = scalerr[names(data.raw),]
  nrow_dataset = nrow(data.raw)
  for(i in 1:nrow_dataset){
    dataset.norm[i,]=(data.raw[i,]-scal[,'mean'])/scal[,'sd']
  }
  return(dataset.norm)
}
traindata_norm = norm_fun(train_data[,-c(1,1501:1502)],scaler)
trainData2 = cbind(train_data[,c(1,1501:1502)],traindata_norm)

testdata_norm = norm_fun(test_data[,-c(1,1501:1502)],scaler)
testdata = cbind(test_data[,c(1,1501:1502)],testdata_norm)
#Feature Selection===============================================================================
#Correlation--------------------------------------------------------------------
library(corrplot)
# calculate correlation matrix
corMatMy <- cor(trainData2[, -1])
highlyCor <- colnames(trainData2[, -1])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
train_data_cor <- trainData2[, which(!colnames(trainData2) %in% highlyCor)]

#Boruta-feature selection (i.e. pick important variables) using Boruta---
library(Boruta)
set.seed(123)
results_boruta <- Boruta(x = trainData2[, -1], y = as.factor(trainData2$group),doTrace = 2)
final.boruta <- TentativeRoughFix(results_boruta)
finalvars = getSelectedAttributes(final.boruta, withTentative = F)
train_data_boruta <- trainData2[, c(1, which(colnames(trainData2) %in% finalvars))]

#Recursive Feature Elimination (RFE)---------------------------------------------
set.seed(7)
library(randomForest)
results_rfe <- rfe(x = trainData2[, -1], 
                   y = as.factor(trainData2$group), 
                   sizes = 2^(2:10), 
                   rfeControl = rfeControl(functions = rfFuncs, method = "repeatedcv", number = 10))

# chosen features
predictors(results_rfe)
train_data_rfe <- trainData2[, c(1, which(colnames(trainData2) %in% predictors(results_rfe)))]

#Genetic Algorithm (GA)-----------------------------------
model_ga <- gafs(x = trainData2[, -1], 
                 y = as.factor(trainData2$group),
                 iters = 20, 
                 popSize = 20, 
                 gafsControl = gafsControl(functions = rfGA, 
                                           method = "cv",    
                                           genParallel = TRUE, 
                                           allowParallel = TRUE))
train_data_ga <- trainData2[, c(1, which(colnames(trainData2) %in% model_ga$ga$final))]

#feature selection using MXM----------------------
library(MXM)
feature = NULL;
for (i in 1:nrow(trainData2)) {
  mmpco <- MXM::MMPC( target  = as.factor(trainData2[-i,]$group),            
                      dataset = as.matrix(trainData2[-i, -1]),            
                      max_k = 3,          
                      threshold = 0.05,                                         
                      test = 'testIndFisher',   
                      ini = NULL,                                                
                      hash =  TRUE,      
                      hashObject = NULL,                                        
                      ncores = 3,         
                      backward = TRUE)  
  feature_temp = colnames(trainData2[,-1])[mmpco@selectedVarsOrder]
  feature = c(feature,feature_temp)
}
var1 = as.data.frame(table(feature))
var2 = subset(var1,Freq >= nrow(trainData2)*0.2)
train_data_mmpc <- trainData2[, c(1, which(colnames(trainData2) %in% var2$feature))]

#Machine Learning packages======================================================================================================================
# configure multicore
cl <- makeCluster(detectCores())
registerDoParallel(cl)

#Define the training control-----------
fitControl <- caret::trainControl(
  method = "cv", 
  number = 5,  
  returnResamp = 'final',
  verboseIter = FALSE, 
  allowParallel = TRUE, 
  savePredictions = 'all',  
  classProbs = T,
  summaryFunction=multiClassSummary
)

#model comparation function============================
multiple_model_comparation = function(dataset,methodlist,fitcontrols,modelnames_tail) {
  modellist = list(NULL)
  for(i in 1:length(methodlist)){
    tryCatch({
      print(methodlist[i])
      model_train.tempt = caret::train(group ~ .,data = dataset,method = methodlist[i],trControl = fitcontrols,weights = weight_groupclass)
      modellist[[i]]=model_train.tempt
      names(modellist)[i]=paste0(methodlist[i],'.',modelnames_tail)
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(modellist)
}
# Define the model methodlist-------------------------------------------
method.list = c('rf','kknn','C5.0Tree','pls')
#model with Cor selected features====================
set.seed(42)
modellist.traindat_cor = multiple_model_comparation(train_data_cor,method.list,fitControl,'traindata_cor')

#model with Recursive Feature Elimination selected features========
set.seed(43)
modellist.traindat_rfe = multiple_model_comparation(train_data_rfe,method.list,fitControl,'traindata_rfe')

#model with Genetic Algorithm selected features===================================================
set.seed(44)
modellist.traindat_ga = multiple_model_comparation(train_data_ga,method.list,fitControl,'traindata_ga')

#model with boruta selected features===================================================
set.seed(45)
modellist.traindat_boruta = multiple_model_comparation(train_data_boruta,method.list,fitControl,'traindata_boruta')

#model with mmpc selected features===================================================
set.seed(46)
modellist.traindat_mmpc = multiple_model_comparation(train_data_mmpc,method.list,fitControl,'traindata_mmpc')

#model with all features===================================================
set.seed(47)
modellist.traindat_all = multiple_model_comparation(trainData2,method.list,fitControl,'traindata_all')

#model comparisons=====
model_list_traindata = c(modellist.traindat_all,modellist.traindat_cor,#
               modellist.traindat_boruta,modellist.traindat_mmpc,
               modellist.traindat_rfe,#
               modellist.traindat_ga)
model_list_traindata = model_list_traindata[-which(sapply(model_list_traindata, is.null))]
resample_models_traindat <- resamples(model_list_traindata)
resample_results = resample_models_traindat

bwplot(resample_results , metric = c('prAUC','logLoss','Accuracy'),layout = c(3,1),xlim=c(0,1.2),
       ylab ='Multi-omic models comparisons in train_data')

# similar to value or rank on the leaderboard
#selected models--************************************************************----------------------------------------
RF_CA = model_list_traindata[[match('rf.traindata_cor',names(model_list_traindata))]]
plot(varImp(RF_CA),top = 10,main = 'Top 10 features of multiomic RF_CA')#4*5
#model features
feature_imp_fun = function(modellists,selectedmodels) {
  model.list = modellists
  selected_models = selectedmodels
  features.df = data.frame()
  for(i in 1:length(selected_models)){
    tryCatch({
      df_vap = varImp(model.list[[match(selected_models[i],names(model.list))]])$importance; 
      df_vap$model = selected_models[i]
      df_vap$featurelabel = str_split_fixed(selected_models[i],"_",2)[2]
      df_vap$features =gsub("[`]","",rownames(df_vap))
      df_vap$algorithm = str_split_fixed(selected_models[i],"[.]",2)[1]
      df_vap$rank = rank(-df_vap$Overal)
      features.df = rbind(features.df,df_vap[,c('features','model','algorithm','featurelabel','Overall','rank')])
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(features.df)
}
feature_df = feature_imp_fun(model_list_traindata,'rf.traindata_cor'); head(feature_df)
feature_df.top10 = subset(feature_df,rank<11);dim(feature_df.top10); dim(feature_df.top10)
#==========================
final_model_fun = function(datast,featuredf,fitcontrols) {
  modellist = list(NULL)
  dataset = datast
  df.feature = featuredf
  feature.model = unique(df.feature$model)
  for(i in 1:length(feature.model)){
    tryCatch({
      features = c(c('group'),df.feature[df.feature$model==feature.model[i],'features'])
      meth = unique(df.feature[df.feature$model==feature.model[i],'algorithm'])
      model_train.tempt = caret::train(group ~ .,data = dataset[,features],method = meth,trControl = fitcontrols,preProcess = NULL,tuneLength=50,weights = weight_groupclass)
      bestparameter = model_train.tempt$bestTune
      model_train = caret::train(group ~ .,data = dataset[,features],method = meth, trControl = fitcontrols,tuneGrid = expand.grid(bestparameter),weights = weight_groupclass)
      modellist[[i]]=model_train
      names(modellist)[i]=feature.model[i]
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(modellist)
}

set.seed(145)
modellist.final = final_model_fun(trainData2,feature_df.top10,fitControl)

#
confution_df_fun = function(modellists,featurelist,datast) {
  model.list = modellists
  feature.list = featurelist
  dataset = datast
  con_matrix = data.frame()
  for(i in names(model.list)){
    tryCatch({
      features = c(c('group'),feature.list[[match(i,names(feature.list))]])
      model = model.list[[match(i,names(model.list))]]
      val.datast = dataset[,features]
      predictions <- predict(model, val.datast)
      con.matrix <- as.data.frame(table(predictions, val.datast$group))
      con.matrix$goodbad = ifelse(con.matrix$predictions == con.matrix$Var2, "good", "bad")
      con.matrix$prop = con.matrix$Freq/sum(con.matrix$Freq)
      con.matrix$model = i;con.matrix$featureslabel =i
      con_matrix = rbind(con_matrix,con.matrix)
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(con_matrix)
}
confusion_matrix.batch = confution_df_fun(modellist.final,feature_list,testdata)
confusion_matrix.batch$model = as.factor(confusion_matrix.batch$model)
#confution_plots
ggplot(data = confusion_matrix.batch, mapping = aes(x = predictions,y = Var2, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "blue", bad = "red"))+
  theme_bw()+facet_wrap(confusion_matrix.batch$model)+
  xlim(rev(levels(confusion_matrix.batch$Var2)))
#followup cohort annotation-----
annot_fu<-protein.run.df[protein.run.df$group=="M.fu"|protein.run.df$group=="S.fu"|protein.run.df$group=="F.fu",];dim(annot_fu)
annot_fu$group2 = ifelse(annot_fu$group=='M.fu','M',ifelse(annot_fu$group=='S.fu','S','F'))
annot_fu.val = subset(annot_fu,group2!='F')
annot_fu.val$sample = paste0('sample.',annot_fu.val$participant)
#------
metabolomics.followup = read.csv('concatenated.norm.uniq.followup.sample.uniqmetabolite.manual.annotation.csv');
rownames(metabolomics.followup)=metabolomics.followup[,2]
metabolomics.followup = log2(metabolomics.followup[,-c(1:2)]); metabolomics.followup[1:3,1:3]
metabo.fu = as.data.frame(t(metabolomics.followup[,annot_fu.val$sample])); dim(metabo.fu)
min.met = min(apply(metabo.fu,2,function(x){min(na.omit(x))})); min.met
metabo.fu.imput = as.matrix(Hmisc::impute(metabo.fu,(min.met-1)))
metabo.fu.imput = as.data.frame(metabo.fu.imput)
#-------
proteomics.followup = read.csv('dat.sample.followup.csv');
rownames(proteomics.followup)=proteomics.followup[,1]
proteomics.followup = proteomics.followup[,-1]; proteomics.followup[1:3,1:3]
colnames(proteomics.followup)=protein.run.df[colnames(proteomics.followup),'samples']; proteomics.followup[1:3,1:3];dim(proteomics.followup)
proteomic.fu = as.data.frame(t(proteomics.followup[,annot_fu.val$samples])); dim(proteomic.fu)
proteomic.fu[proteomic.fu==0]=NA; table(proteomic.fu==0)
min.prot = min(apply(proteomic.fu,2,function(x){min(na.omit(x))}))
proteomic.fu.imput = as.matrix(Hmisc::impute(proteomic.fu,(min.prot-1)))
proteomic.fu.imput = as.data.frame(proteomic.fu.imput)
#pheno df
fu.participant = data.frame(annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'samples'],annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'participant0'],str_split_fixed(annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'participant0'], "_", 2))
names(fu.participant) = c('filename','sample','paticipant','time'); rownames(fu.participant)=fu.participant$filename
fu.pheno = data.frame(participant = c('329', '333', '326', '328', '330', '331', '332'),age=c(56, 43, 57, 82, 83, 45, 61),sex = c('M', 'M', 'F', 'M', 'F', 'F', 'M'))
fu.pheno.dummy = fastDummies::dummy_cols(fu.pheno, select_columns ='sex'); rownames(fu.pheno.dummy) = fu.pheno.dummy$participant
fu.participant.pheno = cbind(fu.participant[,-1],fu.pheno.dummy[fu.participant$paticipant,2:5])
#integrate omic matrix-------------------------------
df.val.raw = cbind(metabo.fu.imput,proteomic.fu.imput,fu.participant.pheno[,4:7]); table(is.na(df.val.raw))
df.val.raw$group = annot_fu.val[,'group2']; table(rownames(df.val.raw)==annot_fu.val$sample)
df.val.raw$group = factor(df.val.raw$group,levels = c('H','M','S'))
df.val_norm = norm_fun(df.val.raw[,-c(1590:1593)],scaler)
df.val = cbind(df.val.raw[,c(1590:1593)],df.val_norm)
#-----------------------------------------------------------------------------
predict_prob_fun = function(modellists,featurelist,datast) {
  model.list = modellists
  feature.list = featurelist
  dataset = datast
  df.predictprob = data.frame()
  for(i in names(model.list)){
    tryCatch({
      features = c(c('group'),feature.list[[match(i,names(feature.list))]])
      model = model.list[[match(i,names(model.list))]]
      val.datast = dataset[,features]
      prediction.prob <- predict(model, val.datast,type = 'prob')
      df.rt = data.frame(prediction = prediction.prob[,2:3],reference = val.datast$group)
      samples = annot_fu.val[match(rownames(df.rt),annot_fu.val$sample),'participant0']
      df.rt$sample = samples
      df.rt$model = i
      df.predictprob = rbind(df.predictprob,df.rt)
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(df.predictprob)
}
#
df.prediction = predict_prob_fun(modellist.final,feature_list,df.val)
df.prediction$model = as.factor(df.prediction$model)
#scores barplots
df.bar = df.prediction; dim(df.bar)
df.bar$participant = paste0('P',substr( df.bar$sample,start=1,stop=3))
df.bar$participant = factor(df.bar$participant,levels = c('P333','P332','P330','P326','P328','P329','P331'))
df.bar$time = paste0('T',substr( df.bar$sample,start=5,stop=5))
df.bar$scores = df.bar$prediction.S-df.bar$prediction.M
ggplot(df.bar, aes(fill=time, y=scores, x=participant)) + theme_bw()+
  ggtitle('RF classifiers different multiomic feature panels')+theme(plot.title = element_text(hjust = 0.5))+
  geom_bar(position="dodge", stat="identity",width = 0.5)+ylab('prediction.prob.S-prediction.prob.M')+
  scale_fill_viridis(discrete = T,option = 'C')+
  theme_bw()+facet_wrap(df.bar$model)
#########function to return pr or roc curve lines and mean micro-averaged auc-----------------------------
pr_roc_auc_fun = function(inputmodellist,auctype_pr.or.roc){
  library(ROCR) # for ROC curves
  plotmodellist = inputmodellist
  auctype = auctype_pr.or.roc
  if(auctype_pr.or.roc=='pr'){
    plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
         ylab="Precision",
         xlab="Recall",
         bty='n')
  }else{
    plot(x=NA, y=NA, xlim=c(0,1), ylim=c(0,1),
         ylab="True positive rate",
         xlab="False positive rate",
         bty='n')
  }
  
  lines.list = list()
  AUPR_mean_micro <- vector("list", length(plotmodellist)) # store AUCs
  micro.average.auc = vector("list", length(plotmodellist))
  for (i in seq_along(plotmodellist)){
    plotmodel = plotmodellist[[i]]
    cvfold = names(table(plotmodel$pred$Resample))
    micro_AUCPR_cvfold = vector("list",length(cvfold))
    score.initeration = list()
    labels.ininteration = list()
    for (j in seq_along(cvfold)){
      df.pred.cvfold = plotmodel$pred[plotmodel$pred$Resample==cvfold[j],]
      response = df.pred.cvfold$obs
      scorelist = vector("list",length(levels(response)))
      labellist = vector("list",length(levels(response)))
      auc.class.infold = list()
      for (k in seq_along(levels(response))) {
        cur.class <- levels(response)[k] #positive class
        score = df.pred.cvfold[,cur.class]# posterior for  positive class
        labels = ifelse(response==cur.class,1,0)
        pred <- prediction(score, labels)
        if(auctype=='pr'){perf = performance(pred, "prec", "rec")}else{perf = performance(pred, "tpr", "fpr")}
        roc.x <- unlist(perf@x.values)#recall
        roc.y <- unlist(perf@y.values)#precision
        lines(roc.y ~ roc.x, lwd = 1)
        #test
        line = cbind(roc.x,roc.y)
        lines.list[[i*j*k]]=line#####
        #store score and label for class k
        scorelist[[k]] = score
        labellist[[k]] = labels
        if(auctype=='pr'){auc.class = MLmetrics::PRAUC(y_pred = score, y_true = labels)} else
        {auc.class = MLmetrics::AUC(y_pred = score, y_true = labels)}
        #auc.fold <- performance(pred.fold, "aucpr")
        auc.class.infold[[k]] <- auc.class
      }
      # get micro-average AUPR in fold j in i th iteration
      scores.fold = unlist(scorelist)
      labels.fold = unlist(labellist)
      if(auctype=='pr'){auc.fold = MLmetrics::PRAUC(y_pred = scores.fold, y_true = labels.fold)} else
      {auc.fold = MLmetrics::AUC(y_pred = scores.fold, y_true = labels.fold)}
      #micro.averaged_AUCPR <- auc.fold
      micro_AUCPR_cvfold[[j]] = unlist(auc.class.infold)
      score.initeration[[j]] = scores.fold
      labels.ininteration[[j]] = labels.fold
    }
    scores.iteration = unlist(score.initeration)
    labels.iteration = unlist(labels.ininteration)
    if(auctype=='pr'){auc.iter = MLmetrics::PRAUC(y_pred = scores.iteration, y_true = labels.iteration)} else
    {auc.iter = MLmetrics::AUC(y_pred = scores.iteration, y_true = labels.iteration)}
    micro.average.auc[[i]] = auc.iter
    AUPR_mean_micro[[i]] = unlist(micro_AUCPR_cvfold)
  }
  
  #print mean AUPR and 95% confidence interval
  mean_AUPR = mean(unlist(AUPR_mean_micro))
  sd = sd(unlist(AUPR_mean_micro))
  n = length(unlist(AUPR_mean_micro))
  margin = qt(0.975,df=n-1)*sd/sqrt(n)
  low = mean_AUPR-margin
  high = mean_AUPR+margin
  lines(x=c(0,1), c(0,1))
  if(auctype_pr.or.roc=='pr'){
    legend("bottomleft", paste0("Mean AUPR = ",round(mean_AUPR,4),"\n","95% CI = ",round(low,4)," - ",round(high,4)))
  }else{
    legend("bottomleft", paste0("Mean AUROC = ",round(mean_AUPR,4),"\n","95% CI = ",round(low,4)," - ",round(high,4)))
  }
  
  return(list(metrics = cbind(mean.AUC = mean_AUPR,CI95low = low, CI95high = high),lines = lines.list,
              auc.list = AUPR_mean_micro,micro.auc.list = micro.average.auc))
}
#rasterized densityplot function-------------
plot_tune_roc.or.pr = function(iteration.model.list,plotypte.roc.or.pr){
  iteration_model_list = iteration.model.list
  plottype = plotypte.roc.or.pr
  pt = pr_roc_auc_fun(iteration_model_list,plottype)
  lines.obj = pt[[2]]
  lines.obj.re_null = lines.obj[-which(sapply(lines.obj, is.null))]
  lines.obj.re_null.reNA = sapply(lines.obj.re_null, function(x) na.omit(x));lines.obj.re_null.reNA
  library(raster)
  prlines = spLines(lines.obj.re_null.reNA); plot(prlines)
  r <- raster(ncols=200, nrows=200)
  extent(r) = c(0,1.1,0,1.1)
  r <- rasterize(prlines, r,fun='count')
  fun_color_range <- colorRampPalette(c("blue", "red"))   
  my_colors <- fun_color_range(100) 
  if(plotypte.roc.or.pr=='pr'){
    plot(r, col=my_colors,
         ylab="Precision", xlab="Recall",bty='n',lwd = 0.51)
    legend("bottomleft", paste0("Mean AUPR = ",round(pt$metrics[1],4),"\n","95% CI = ",round(pt$metrics[2],4)," - ",round(pt$metrics[3],4)))
  }else{
    plot(r, col=my_colors,
         ylab="True positive rate", xlab="False positive rate",bty='n',lwd = 0.51)
    legend("bottomleft", paste0("Mean AUROC = ",round(pt$metrics[1],4),"\n","95% CI = ",round(pt$metrics[2],4)," - ",round(pt$metrics[3],4)))
  }
  
}
#validation by 5-fold cv and iterations--------------
fit.Control <- caret::trainControl(
  method = "cv", 
  number = 5, 
  returnResamp = 'all',
  verboseIter = FALSE, 
  allowParallel = TRUE, 
  savePredictions = 'all',       
  classProbs = T,                 
  summaryFunction=multiClassSummary
)
#iteraction----------------
iteration_model_list = list()
for (i in 1:100){
  iteration_model_list[[i]] = caret::train(group ~., data = trainData2[,c(c('group'),feature_list[[1]])],
    trControl = fit.Control,
    tuneLength = 20,
    weights = weight_groupclass,
    method = "rf",
    verbose = TRUE
  )
}
#-
length(iteration_model_list)
#compute metrics---------------
zz = pr_roc_auc_fun(iteration_model_list,'pr'); hist(unlist(zz$auc.list)); median(unlist(zz$auc.list))
ww = pr_roc_auc_fun(iteration_model_list,'roc'); hist(unlist(ww$auc.list)); median(unlist(ww$auc.list))
#plot get the PR-roc lines objects from defined function
pr.tuneplot = plot_tune_roc.or.pr(iteration_model_list,'pr')
roc.tuneplot = plot_tune_roc.or.pr(iteration_model_list,'roc')
#----------------------------------------------------------------------------------------------------
sessionInfo()

