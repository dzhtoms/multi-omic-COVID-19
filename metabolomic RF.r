library(caret)
library(tidyverse) 
library(readr)
library(mice) 
library(viridis)
library(doParallel)
library(MLmetrics)

source('multi-omic covid severity classifiers.r')
#-----------------------
trainData.m = trainData2[,c(1:4,632:1502)]; dim(trainData.m)
corMatMy <- cor(trainData.m[, -1])
highlyCor <- colnames(trainData.m[, -1])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
train_data_cor.m <- trainData.m[, which(!colnames(trainData.m) %in% highlyCor)]
dim(train_data_cor.m)
#metabolite RF model
m.rf.ca = caret::train(group ~ .,data = train_data_cor.m,method = 'rf',trControl = fit.Control,preProcess = NULL,weights = weight_groupclass)
varimp.m = varImp(m.rf.ca)$importance; varimp.m$rank = rank(-varimp.m$Overall); varimp.m$name = gsub("[`]","",rownames(varimp.m));head(varimp.m)
m.features = varimp.m[varimp.m$rank<11,'name']
retrained.m.rf.tune = caret::train(group ~ .,data = trainData.m[,c(c('group'),m.features)],method = 'rf',trControl = fit.Control,preProcess = NULL,tuneLength=20,weights = weight_groupclass)
retrained.m.rf = caret::train(group ~ .,data = trainData.m[,c(c('group'),m.features)],method = 'rf',trControl = fit.Control,preProcess = NULL,tuneGrid = expand.grid(mtry = 3),weights = weight_groupclass)
confusion_matrix.batch = confution_df_fun(retrained.m.rf,m.features,testdata)
confusion_matrix.batch$model = as.factor(confusion_matrix.batch$model)
#confution_plots
ggplot(data = confusion_matrix.batch, mapping = aes(x = predictions,y = Var2, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "blue", bad = "red"))+
  theme_bw()+facet_wrap(confusion_matrix.batch$model)+
  xlim(rev(levels(confusion_matrix.batch$Var2)))#4*8
#disease clock
metabolomics.followup = read.csv('concatenated.norm.uniq.followup.sample.uniqmetabolite.manual.annotation.csv');
rownames(metabolomics.followup)=metabolomics.followup[,2]
metabolomics.followup = log2(metabolomics.followup[,-c(1:2)]); metabolomics.followup[1:3,1:3]
metabo.fu = as.data.frame(t(metabolomics.followup[,annot_fu.val$sample])); dim(metabo.fu)
min.met = min(apply(metabo.fu,2,function(x){min(na.omit(x))})); min.met
metabo.fu.imput = as.matrix(Hmisc::impute(metabo.fu,(min.met-1)))
metabo.fu.imput = as.data.frame(metabo.fu.imput)
#pheno df
fu.participant = data.frame(annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'samples'],
                            annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'participant0'],
                            str_split_fixed(annot_fu.val[match(rownames(metabo.fu),annot_fu.val$sample),'participant0'], "_", 2))
names(fu.participant) = c('filename','sample','paticipant','time'); rownames(fu.participant)=fu.participant$filename
fu.pheno = data.frame(participant = c('329', '333', '326', '328', '330', '331', '332'),age=c(56, 43, 57, 82, 83, 45, 61),sex = c('M', 'M', 'F', 'M', 'F', 'F', 'M'))
fu.pheno.dummy = fastDummies::dummy_cols(fu.pheno, select_columns ='sex'); rownames(fu.pheno.dummy) = fu.pheno.dummy$participant
fu.participant.pheno = cbind(fu.participant[,-1],fu.pheno.dummy[fu.participant$paticipant,c(2,4,5)])
#integrate omic matrix-------------------------------
df.val.raw = cbind(metabo.fu.imput,fu.participant.pheno[,4:6]); table(is.na(df.val.raw))
df.val.raw$group = annot_fu.val[,'group2']; table(rownames(df.val.raw)==annot_fu.val$sample)
df.val.raw$group = factor(df.val.raw$group,levels = c('H','M','S'))
df.val_norm = norm_fun(df.val.raw[,-c(873:875)],scaler)
df.val = cbind(df.val.raw[,c(875,873:874)],df.val_norm); df.val[1:3,1:6]
#
df.prediction = predict_prob_fun(retrained.m.rf,m.features,df.val)
df.prediction$model = as.factor(df.prediction$model)
#scores barplots
df.bar = df.prediction; dim(df.bar)
df.bar$participant = paste0('P',substr( df.bar$sample,start=1,stop=3))
df.bar$participant = factor(df.bar$participant,levels = c('P333','P332','P330','P326','P328','P329','P331'))
df.bar$time = paste0('T',substr( df.bar$sample,start=5,stop=5))
df.bar$scores = df.bar$prediction.S-df.bar$prediction.M
ggplot(df.bar, aes(fill=time, y=scores, x=participant)) + theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_bar(position="dodge", stat="identity",width = 0.5)+ylab('prediction.prob.S-prediction.prob.M')+
  scale_fill_viridis(discrete = T,option = 'C')+
  theme_bw()+facet_wrap(df.bar$model)
#retrained model comparsons------------
retrained.modellist = list(multiomic.rRF = modellist.final$rf.traindata_cor,proteomic.rRF = retrained.p.rf,metabolomic_rRF = retrained.m.rf)
resample_models_revise <- resamples(retrained.modellist)
bwplot(resample_models_revise , metric = c('Accuracy','logLoss','Mean_F1'),layout = c(3,1),xlim=c(0.85,1.05),
       ylab ='Top10_feature based retrained RF in train_data')

#----------------------------------------------------------------------------------------------------
sessionInfo()

