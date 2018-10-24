
library(pROC) 
library(factoextra)
library(FactoMineR)
library(cluster)
library(dendextend)
library(ggsci)
library(dplyr)
library(class)
library(BKPC)
library(kernlab)
library(ggplot2)
library(purrr)
library(tidyr)
library(broom)
library(sigr)
library(WVPlots)
library(vtreat)
library(caret)
library(dbscan)
library(fpc)
library(ecodist)
library(ISLR)

#Preparing the data
data <-read.table("wholedata.txt",header = TRUE, sep = ",")
#Selecting the continous predictors
datas <- data[4:16]
#Rownames
row.names(datas) <- data$sj
#Standarizing the dataset
data.scale<-scale(datas)


#PCA ANALYSIS by Princomp
pca.data2<-princomp(data.scale)
#PCA Analysis by PCA
pca.data<-PCA(data.scale)

#Elbow plot
fviz_screeplot(pca.data)
#####################################################################################
#Distance matrix
dist<-distance(data.scale,method = "mahalanobis")
#Clustering object
data.clust<-hclust(dist,method = "complete")

#Hierarchycal clustering
dendh <- as.dendrogram(data.clust) %>% set("branches_k_color", value = 2:5, k = 3) %>%
  set("labels",data$sj)%>%  plot(main = "Cancer clustering by Hierarchical methods")
#####################################################################################
#KMeans
data.km3<-kmeans(data.scale,3,iter.max = 50)
fviz_cluster(data.km3,data.scale)
#####################################################################################
#Silhouettes plot

avg_sil <- function(k) {
  km.res <- kmeans(data.scale, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster,dist)
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "o", pch = 15, frame = FALSE, 
     col= 10,
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")

# Prediction

scores<-pca.data2$scores[,1:8]

#sample.data<-sample(1:nrow(abs(scores)),floor(0.70*nrow(abs(scores))),replace=F)
sample.data<-sample(1:nrow(data.scale),floor(0.70*nrow(data.scale)),replace=F)
train<-data[sample.data,]%>% select(-sj,-G,-Emotion)
test<-data[-sample.data,]%>% select(-sj,-G,-Emotion)
labels<-data[sample.data,]$Emotion
label.test<-data[-sample.data,]$Emotion
#Maximun 60% for k =11

accuracy <- function(mc){
  ac<-diag(mc)/sum(mc)
  return(sum(ac))
}

knnc=knn(train,test,labels,k=11)
mc=table(as.matrix(label.test),as.matrix(knnc))
accuracy(mc)



#No linear PCA
#New dataset
data.fw<-datas%>%select(scrm,scrstd,scrdr,hrrmssd)
data.bw<-datas%>%select(scrstd,scrdr,hrm,hrmtime,hrrmssd)
kpcData<-kPCA(data.scale)
kpcData.fw <- kPCA(scale(data.fw))
kpcData.bw <- kPCA(scale(data.bw))

plot(kpcData.fw$Es/sum(kpcData.fw$Es), pch=6,col="blue",main="Principal Components by
     using the extrated features from Stepwise Regression")
points(kpcData.bw$Es/sum(kpcData.bw$Es), pch=17, col="red")
points(kpcData$Es/sum(kpcData$Es), pch=12, col="green")

#Computing eigen values for each KPCA 

e.fw<-kpcData.fw$Es/sum(kpcData.fw$Es)
eig.fw<-sum(e.fw[1:3])
e.bw<-kpcData.bw$Es/sum(kpcData.bw$Es)
eig.bw<-sum(e.bw[1:3])
e.full<-kpcData$Es/sum(kpcData$Es)
eig.full<-sum(e.full[1:3])

legend("right", c("Backward", "Forward","Whole variables"),
       col = c("red", "blue","green"), lty = c(1, 2,3), ncol = 1,
       cex = 0.7)
#KPCA compontents of BW method
kpcs.bw<-kpcData.bw$KPCs[,1:3]

#Splitting the data set into 70% train and the rest test
sample<-sample(1:nrow(kpcs.bw),floor(0.70*nrow(kpcs.bw)),replace=F)
train<-data.scale[sample,]
test<-data.scale[-sample,]
label.test<-data[-sample,]$Emotion
lab<-data[sample,]
labels=lab$Emotion

# 67.85% con k=20
  knn.kpca<-knn(train,test,labels,k=20)
  mc.kpca<-table(as.matrix(label.test),as.matrix(knn.kpca))
  accuracy(mc.kpca)
  
#####################################################################################
#Cross-Validation
set.seed(100)
data.fs<-data %>% select(-sj,-G,hrmtime)
indxTrain <- createDataPartition(y = data.fs$Emotion,p = 0.70,list = FALSE)
training <- data.fs[indxTrain,]
testing  <- data.fs[-indxTrain,]
control<-trainControl(method="repeatedcv",repeats = 10,savePredictions = TRUE)


#knn
model<-train(Emotion~.,data=training,method="knn",preProcess=c("center","scale"),
             trControl=control, tuneLength=50)

#Testing 
test_pred <- predict(model, newdata = testing)
cm.caret<-confusionMatrix(test_pred,testing$Emotion)
cm.caret

#multinomial regression
mreg<-train(Emotion~.,data=training,method="multinom",preProcess=c("center","scale"),
            family=binomial(link="logit"), trControl=control,tuneLenght=50 )

pred.multinom<- predict(mreg, newdata = testing)
multinom.caret<-confusionMatrix(pred.multinom,testing$Emotion)

#knn 81% Bno
#svmpoly 74%
#svmRadial 62%
#svmRadialCost 70.37%
#naive_bayes 81,48% Bno
#stepLDA 77.78% Bno
svmr<-train(Emotion~.,data=training,method="adaboost",preProcess=c("center","scale"),
            trControl=control)

pred.svm<- predict(svmr, newdata = testing)
svm.caret<-confusionMatrix(pred.svm,testing$Emotion)

#rfe
rfe<-train(Emotion~.,data=training,method="rf",preProcess=c("center","scale"),
            trControl=control)
pred.rf<-predict(rfe, newdata = testing)
rf.caret<-confusionMatrix(pred.rf,testing$Emotion)


steplda<-train(Emotion~.,data=training,method="stepLDA",preProcess=c("center","scale"),
            trControl=control)
pred.steplda<-predict(steplda, newdata = testing)
adab.caret<-confusionMatrix(pred.steplda,testing$Emotion)

rmda<-train(Emotion~.,data=training,method="gaussprLinear",preProcess=c("center","scale"),
               trControl=control)
pred.rmda<-predict(rmda, newdata = testing)
rmda.caret<-confusionMatrix(pred.rmda,testing$Emotion)




 #Plot K-Neighboors vs Accuracy
info<-model$results
ggplot(info,aes(x=info$k,y=info$Accuracy))+geom_line()+geom_point()

###Caret KNN with Stepwise Backward model
Em<-data$Emotion
data.bw.knn<-data.bw%>%mutate(Em)
indxTrain.kpca <- createDataPartition(y = data.bw.knn$Em,p = 0.70,list = FALSE)
training.kpca <- data.bw.knn[indxTrain.kpca,]
testing.kpca  <- data.bw.knn[-indxTrain.kpca,]
control.kpca<-trainControl(method="repeatedcv",repeats = 3)
model.kpca<-train(Em~.,data=training.kpca,method="knn",
             trControl=control.kpca, tuneLength=50)

test_pred.kpca <- predict(model.kpca, newdata = testing.kpca)
confusionMatrix(test_pred.kpca,testing$Em)


#####################################################################################
#Forward selection
FitStart<-lm(scrm~1,data=data[-1:-3])
tmodel<-lm(scrm~.,data=data[-1:-3])
#Backward selection
bm<-step(tmodel,direction="backward") #Step:  AIC=-144.49
fm<-step(FitStart,direction="forward",scope=formula(tmodel)) #AIC=-145.59

#####################################################################################
#Logistic regression

#Discrimining Amuse and Sad excluding Neutral
no.neutral.df<-data

no.neutral.df$dummy.amuse<-ifelse(no.neutral.df[,3]=="A",1,0)

amuse.model<-no.neutral.df%>%glm(formula=dummy.amuse~scrm+scrstd+scrdr+scraonv+scrpnv+
                   hrm+hrstd+hrdr+hrmode+hrssdn+hrrmssd,family="binomial")

no.neutral.df$amuse.prob<-predict(amuse.model,newdata=NULL,type="response")

no.neutral.df$amuse.pred<-ifelse(no.neutral.df$amuse.prob>mean(no.neutral.df$dummy.amuse),1,0)

mean(no.neutral.df$dummy.amuse==no.neutral.df$amuse.pred)

roc.amuse<-roc(no.neutral.df$dummy.amuse,no.neutral.df$amuse.pred)

auc.amuse<-auc(roc.amuse)

##############################################################################

#Discrimining by Sad-Neutral excluding Amuse
no.amuse.df<-data

no.amuse.df$dummy.sad<-ifelse(no.neutral.df[,3]=="S",1,0)

sad.model<-no.amuse.df%>%glm(formula=dummy.sad~scrm+scrstd+scrdr+scraonv+scrpnv+
                                   hrm+hrstd+hrdr+hrmode+hrssdn+hrrmssd,family="binomial")

no.amuse.df$sad.prob<-predict(sad.model,newdata=NULL,type="response")
  
no.amuse.df$sad.pred<-ifelse(no.amuse.df$sad.prob>mean(no.amuse.df$dummy.sad),1,0)

mean(no.amuse.df$dummy.sad==no.amuse.df$sad.pred)

roc.sad<-roc(no.amuse.df$dummy.sad,no.amuse.df$sad.pred)

auc.sad<-auc(roc.sad)


##############################################################################

#Discrimining by Amusement-Neutral excluding sad
no.sad.df<-data

no.sad.df$dummy.neutral<-ifelse(no.sad.df[,3]=="N",1,0)

neutral.model<-no.sad.df%>%glm(formula=dummy.neutral~scrm+scrstd+scrdr+scraonv+scrpnv+
                               hrm+hrstd+hrdr+hrmode+hrssdn+hrrmssd,family="binomial")

no.sad.df$neutral.prob<-predict(neutral.model,newdata=NULL,type="response")

no.sad.df$neutral.pred<-ifelse(no.sad.df$neutral.prob>mean(no.sad.df$dummy.neutral),1,0)

mean(no.sad.df$dummy.neutral==no.sad.df$neutral.pred)

roc.neutral<-roc(no.sad.df$dummy.neutral,no.sad.df$neutral.pred)

auc.neutral<-auc(roc.neutral)

####################################################################

plot(roc.amuse,col="red",main="ROC for the whole dataset")
par(new=TRUE)
plot(roc.neutral,col="blue",type="l",lty=2)
par(new=TRUE)
plot(roc.sad,col="green",type="l",lty=2,lwd=3)

legend("right", c("Amusement", "Neutral","Sad"),
       col = c("red", "blue","green"), lty = c(1, 2,3), ncol = 1,
       cex = 0.7)
#####################################################################


#####################################################################
#Density clustering

data.dc<-as.matrix(kpcs.bw)
dist.dc<-full(distance(data.dc,"manhattan"))
res.fpc <- fpc::dbscan(dist.dc, eps = 1.5, MinPts = 3,method = "dist")
fviz_cluster(res.fpc, dist.dc, geom = "point")
