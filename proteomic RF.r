library(caret)
library(tidyverse) 
library(readr) 
library(mice)
library(viridis)
library(doParallel)
library(MLmetrics)

source('multi-omic covid severity classifiers.r')
df.protein0 = df.prot.imput
annot_filter.p = subset(df.annotation,group!='B');dim(annot_filter.p)
df.protein1 = t(df.protein0);dim(df.protein1);df.protein1[1:3,1:3]
df.protein2 =as.data.frame(df.protein1)
df.protein2$group = df.annotation[rownames(df.protein2),'group2']
df.protein2$age = df.annotation[rownames(df.protein2),'age']
df.protein2$sex = df.annotation[rownames(df.protein2),'sex']; dim(df.protein2)
df.protein3 = df.protein2[,c(629:631,1:628)]; df.protein3[1:3,1:6]
str(df.protein3[,1:6])
df.protein4 = fastDummies::dummy_cols(df.protein3, select_columns ='sex');dim(df.protein4); df.protein4[1:3,c(1:6,631:633)]
df.protein5 = df.protein4[,-3];df.protein5[1:3,c(1:6,631:632)];dim(df.protein5)
df.protein5$group = factor(df.protein5$group,levels = c('H','M','S'))
set.seed(42)
index <- createDataPartition(df.protein5$group, p = 0.7, list = FALSE)
train_data.p <- df.protein5[index, ];dim(train_data.p)
test_data.p  <- df.protein5[-index,];dim(test_data.p)
#Handling class imbalance with weighted methods
weight_groupclass.p = ifelse(train_data.p$group=='H',sum(table(train_data.p$group))/(3*table(train_data.p$group)[1]),
                           ifelse(train_data.p$group=='M',sum(table(train_data.p$group))/(3*table(train_data.p$group)[2]),
                                  sum(table(train_data.p$group))/(3*table(train_data.p$group)[3])))
mean = apply(train_data.p[,-c(1,630:631)],2,function(x) mean(x)); length(mean)
sd = apply(train_data.p[,-c(1,630:631)],2,function(x) sd(x)); length(sd)
scaler.p = data.frame(mean = mean, sd = sd); rownames(scaler.p) = names(train_data.p)[-c(1,630:631)]
#normalization
traindata_norm.p = norm_fun(train_data.p[,-c(1,630:631)],scaler.p)
trainData.p = cbind(train_data.p[,c(1,630:631)],traindata_norm.p)

testdata_norm.p = norm_fun(test_data.p[,-c(1,630:631)],scaler.p)
testdata.p = cbind(test_data.p[,c(1,630:631)],testdata_norm.p)
#cor filter features
corMatMy <- cor(trainData.p[, -1])
highlyCor <- colnames(trainData.p[, -1])[findCorrelation(corMatMy, cutoff = 0.7, verbose = TRUE)]
train_data_cor.p <- trainData.p[, which(!colnames(trainData.p) %in% highlyCor)]
dim(train_data_cor.p)#120*466
#protein RF model
p.rf.ca = caret::train(group ~ .,data = train_data_cor.p,method = 'rf',trControl = fit.Control,preProcess = NULL,weights = weight_groupclass.p)
varimp.p = varImp(p.rf.ca)$importance; varimp.p$rank = rank(-varimp.p$Overall); varimp.p$name = gsub("[`]","",rownames(varimp.p));head(varimp.p)
p.features = varimp.p[varimp.p$rank<11,'name']
retrained.p.rf.tune = caret::train(group ~ .,data = trainData.p[,c(c('group'),p.features)],method = 'rf',trControl = fit.Control,preProcess = NULL,tuneLength=20,weights = weight_groupclass.p)
retrained.p.rf = caret::train(group ~ .,data = trainData.p[,c(c('group'),p.features)],method = 'rf',trControl = fit.Control,preProcess = NULL,tuneGrid = expand.grid(mtry = 2),weights = weight_groupclass.p)

confusion_matrix.batch = confution_df_fun(retrained.p.rf,p.features,testdata)
confusion_matrix.batch$model = as.factor(confusion_matrix.batch$model)
#confution_plots
ggplot(data = confusion_matrix.batch, mapping = aes(x = predictions,y = Var2, fill = goodbad, alpha = prop)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1, fontface  = "bold", alpha = 1) +
  scale_fill_manual(values = c(good = "blue", bad = "red"))+
  theme_bw()+facet_wrap(confusion_matrix.batch$model)+
  xlim(rev(levels(confusion_matrix.batch$Var2)))

#disease clock
#followup cohort annotation-----
annot_fu<-protein.run.df[protein.run.df$group=="M.fu"|protein.run.df$group=="S.fu"|protein.run.df$group=="F.fu",];dim(annot_fu)
annot_fu$group2 = ifelse(annot_fu$group=='M.fu','M',ifelse(annot_fu$group=='S.fu','S','F'))
annot_fu.val = subset(annot_fu,group2!='F')
annot_fu.val$sample = paste0('sample.',annot_fu.val$participant)
#-------
proteomics.followup = read.csv('dat.sample.followup.csv');
rownames(proteomics.followup)=proteomics.followup[,1]
proteomics.followup = proteomics.followup[,-1]; proteomics.followup[1:3,1:3]
colnames(proteomics.followup)=protein.run.df[colnames(proteomics.followup),'samples']; proteomics.followup[1:3,1:3];dim(proteomics.followup)
proteomic.fu = as.data.frame(t(proteomics.followup[,annot_fu.val$samples])); dim(proteomic.fu)
proteomic.fu[proteomic.fu==0]=NA; table(proteomic.fu==0)
min.prot = min(apply(proteomic.fu,2,function(x){min(na.omit(x))})); min.prot
proteomic.fu.imput = as.matrix(Hmisc::impute(proteomic.fu,(min.prot-1)))
proteomic.fu.imput = as.data.frame(proteomic.fu.imput)
#pheno df
fu.participant = data.frame(rownames(proteomic.fu),annot_fu.val[match(rownames(proteomic.fu),annot_fu.val$samples),'participant0'],str_split_fixed(annot_fu.val[match(rownames(proteomic.fu),annot_fu.val$samples),'participant0'], "_", 2))
names(fu.participant) = c('filename','sample','paticipant','time'); rownames(fu.participant)=fu.participant$filename
fu.pheno = data.frame(participant = c('329', '333', '326', '328', '330', '331', '332'),age=c(56, 43, 57, 82, 83, 45, 61),sex = c('M', 'M', 'F', 'M', 'F', 'F', 'M'))
fu.pheno.dummy = fastDummies::dummy_cols(fu.pheno, select_columns ='sex'); rownames(fu.pheno.dummy) = fu.pheno.dummy$participant
fu.participant.pheno = cbind(fu.participant[,-1],fu.pheno.dummy[fu.participant$paticipant,c(2,4:5)])
#proteommic matrix-------------------------------
df.val.raw = cbind(proteomic.fu.imput,fu.participant.pheno[,4:6]); table(is.na(df.val.raw))
df.val.raw$group = annot_fu.val[,'group2']; table(rownames(df.val.raw)==annot_fu.val$samples)
df.val.raw$group = factor(df.val.raw$group,levels = c('H','M','S'))
df.val_norm = norm_fun(df.val.raw[,-c(719:721)],scaler)
df.val = cbind(df.val.raw[,c(721,719,720)],df.val_norm)

df.prediction = predict_prob_fun(retrained.p.rf,p.features,df.val)
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
#----------------------------------------------------------------------------------------------------
sessionInfo()

