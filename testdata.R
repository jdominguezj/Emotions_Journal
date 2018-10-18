library(factoextra)
library(FactoMineR)
library(cluster)
library(dendextend)
library(ggsci)
library(dplyr)
library(class)
library(care)
library(BKPC)
library(kernlab)
library(ggplot2)

set.seed(1000)

data <-read.table("wholedata.txt",header = TRUE, sep = ",")
datas <- data[9:21]
row.names(datas) <- data$sj
data.scale<-scale(datas)


#PCA ANALYSIS
pca.data2<-princomp(data.scale)
pca.data<-PCA(data.scale)
#fviz_pca_biplot(pca.data, col.var = "contrib")+
#  scale_color_gradient2(low="green",mid="blue", high="red", midpoint=5)+theme_minimal()

fviz_pca_biplot(pca.data,
                col.ind = data$Emotion, palette = "aas", legend.title="Emotions",
                col.var="black")

fviz_screeplot(pca.data)

dist<-dist(data.scale,method = "euclidean")
data.clust<-hclust(dist,method = "complete")

dendh <- as.dendrogram(data.clust) %>% set("branches_k_color", value = 2:5, k = 3) %>%
  set("labels",data$sj)%>%  plot(main = "Cancer clustering by Hierarchical methods")

data.km3<-kmeans(data.scale,3,iter.max = 50)
fviz_cluster(data.km3,data.scale)



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


# Statistical info

dplyr::sample_n(datas,66)

v1 <- datas$scrm
v2 <- datas$scrstd
v3 <- data$scrdr
v4 <- data$scravd
v5 <- data$scraonv
v6 <- data$scrpnv
v7 <- data$hrm
v8 <- data$hrstd
# MANOVA test
res.man <- manova(cbind(scrm, scrstd,scrdr,scravd,scraonv,scrpnv,hrm,hrstd) ~ Emotion, data = data)
summary(res.man)


# Prediction

scores<-pca.data2$scores[,1:8]
pc.comp1<-scores[,1]
pc.comp2<-scores[,2]

#sample.data<-sample(1:nrow(abs(scores)),floor(0.70*nrow(abs(scores))),replace=F)
sample.data<-sample(1:nrow(data.scale),floor(0.70*nrow(data.scale)),replace=F)
train<-data[sample.data,]%>% select(-sj,-G,-Emotion,-M,-F,-N,-S,-A)
test<-data[-sample.data,]%>% select(-sj,-G,-Emotion,-M,-F,-N,-S,-A)

X.train=train[,-1]
y.train=train$Emotion

X.test=test[,-1]
y.test=test$Emotion

lab<-data[sample.data,]
tlab<-data[-sample.data,]
label.test<-tlab$Emotion
labels=lab$Emotion

#Maximun 60% for k =11

knnc=knn(train,test,labels,k=15)
mc=table(as.matrix(label.test),as.matrix(knnc))
accuracy(mc)


accuracy <- function(mc){
  ac<-diag(mc)/sum(mc)
  return(sum(ac))
}


#No linear ppca


trk <- as.matrix(data.scale[sample.data ,])
tek <- as.matrix(data.scale[-sample.data ,])

gk<-gaussKern(trk)
Ktrain <- gk$K
image(Ktrain)

gk2 <- gaussKern(trk, tek, gk$theta) 
Ktest <- gk2$K
image(Ktest)

kfunc <- laplacedot(sigma = 1)
Ktrain2 <- kernelMatrix(kfunc, trk)
image(Ktrain2)

# make testing set kernel using kernelMatrix {kernlab}

Ktest2 <- kernelMatrix(kfunc, tek, trk)

kpcData <- kPCA(trk)
kpcKern <- kPCA(Ktrain)
kpcKern2 <- kPCA(Ktrain2)


pairs(kpcData$KPCs[ , 1 : 3], col = data.scale[-sample.data,])

KPCA.new<-kpca(data.scale, kernel = "rbfdot", kpar = list(sigma = 0.1),
               features = 0, th = 1e-4)

plot(rotated(KPCA.new),col=as.integer(data[-sample.data,3]),
     xlab="1st Principal Component",ylab="2nd Principal Component")