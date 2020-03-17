##Load seurat data set, dataMerge_all is a seurat object including all types of cells
load("~/scRNA.RData")

##Subset seurate data and isolate progenitor cells
data.subsetE10E14 <- SubsetData(object = dataMerge_all, subset.raw=T,
                                ident.use = c("E_AP1","E_AP2","L_AP","Mut_L_AP",
                                              "L_AP_mix","E_AP"))
##Extract proginitor cells for analysis
cell_list<-list()
oidentity<-unique(data.subsetE10E14@meta.data$orig.ident)
orig.ident<-data.subsetE10E14@meta.data$orig.ident
orig.ident<-as.data.frame(orig.ident)
colnames(orig.ident)<-"ident"
rownames(orig.ident)<-data.subsetE10E14@raw.data@Dimnames[[2]]
for (i in 1:10)
{
  cell_list[[i]]<-
    WhichCells(data.subsetE10E14,subset.name="orig.ident",accept.value=oidentity[i])
}
data_list<-list()
for(i in 1:10)
{
  data_list[[i]]<-data.use[genes.use, cell_list[[i]], drop = F]
  
}

##Input birthdate weight using published data (Telley et al., 2019)
Temporal_weight<-read.table(file="Temporal_weight_science.txt",sep="\t",head=T,row.names=1)
##"RP23-379C24.2" "Leprel1"       "RP23-14P23.9"  "Mir99ahg"      "Yam1"    were removed, no these genes in my data

##Use temporal fate gene and their data
genes.use_weight<-rownames(Temporal_weight)
genes.use_weight<-intersect(genes.use_weight,data.use@Dimnames[[1]])

####Calculate birthday score for scData
data_list_weight<-list()
for(i in 1:8)
{
  data_list_weight[[i]]<-data.use[genes.use_weight, cell_list[[i]], drop = F]
}
Brithday_score1<-NULL
for(i in 1:8)
{
  Brithday_score<-t(Temporal_weight)%*%as.matrix(data_list_weight[[i]])
  Brithday_score<-t(Brithday_score)
  Brithday_score1<-rbind(Brithday_score1,Brithday_score)
}
Brithday_score1<-Brithday_score1*(-1)
cellid_weight<-orig_ident$orig.ident[match(rownames(Brithday_score1),rownames(orig_ident))]
Brithday_score1<-cbind(Brithday_score1,cellid_weight)
colnames(Brithday_score1)<-c("Birthday_Score","Orig_ID")

##Output of birthdate score
write.table(x=Brithday_score1,"191118_Brithday_score.txt",sep="\t",row.names=T)
Brithday_score1<-read.table(file="Brithday_score.txt")

##Plot of the result
Brithday_score1$Orig_ID=factor(Brithday_score1$Orig_ID,levels=c("E10Con","E10Het","E10Dko",
                                                                "E12het","E12mut","E14Con","E14Het","E14Dko","10X_E13Con"))
beeswarm <- beeswarm(Brithday_score1$Birthday_Score ~ Brithday_score1$Orig_ID,cex=0.03,
                     data = Brithday_score1, method = 'swarm',ylim=c(-100, 150),col=brewer.pal(10, "Set3"),
                     xlab="Orig_ID",
                     ylab="Birthday_Score")
boxplot(Brithday_score1$Birthday_Score ~ Brithday_score1$Orig_ID,boxwex=0.5,col=brewer.pal(10, "Set3"),
        data = Brithday_score1, add = T,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5))


##Input differentiation weight using published data (Telley et al., 2019)
Dif_weight<-read.table(file="differentiation_weight_science.txt",sep="\t",head=T,row.names=1)
genes.use_Difweight<-rownames(Dif_weight)
genes.use_Difweight<-intersect(genes.use_Difweight,data.use@Dimnames[[1]])

####Calculate differentiation score for scData
data_list_weight_Dif<-list()
for(i in 1:10)
{
  data_list_weight_Dif[[i]]<-data.use[genes.use_Difweight, cell_list[[i]], drop = F]
}
Dif_score1<-NULL
for(i in 1:10)
{
  Dif_score<-t(Dif_weight)%*%as.matrix(data_list_weight_Dif[[i]])
  Dif_score<-t(Dif_score)
  Dif_score1<-rbind(Dif_score1,Dif_score)
}
Dif_score1<-as.data.frame(Dif_score1)
Dif_score1$ident<-orig.ident$ident[match(rownames(Dif_score1),rownames(orig.ident))]
colnames(Dif_score1)<-c("Dif_Score","Orig_ID")

##Output of differentiation score
write.table(x=Dif_score1,"191118_Dif_score.txt",sep="\t",row.names=T)
Dif_score1<-read.table(file="Dif_score.txt")

##Plot of the result
Dif_score1$Orig_ID=factor(Dif_score1$Orig_ID,levels=c("E10Con","E10Het","E10Dko",
                                                      "E12het","E12mut","E14Con","E14Het","E14Dko","10X_E13Con"))
beeswarm <- beeswarm(Dif_score1$Dif_Score ~ Dif_score1$Orig_ID,cex=0.03,
                     data = Dif_score1, method = 'swarm',ylim=c(-120, 200),col=brewer.pal(10, "Set3"),
                     xlab="Orig_ID",
                     ylab="Differentiation_Score")
boxplot(Dif_score1$Dif_Score ~ Dif_score1$Orig_ID,boxwex=0.5,col=brewer.pal(10, "Set3"),
        data = Dif_score1, add = T,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5))
