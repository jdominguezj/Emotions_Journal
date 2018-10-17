library(factoextra)
library(FactoMineR)
library(cluster)
library(dendextend)
library(ggsci)

set.seed(1000)

data <-read.table("wholedata.txt",header = TRUE, sep = ",")
datas <- data[4:21]
row.names(datas) <- data$sj
data.scale<-scale(datas)
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

pam.model<-pam(data.scale,k=3)
