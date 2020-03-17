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
write.table(x=Brithday_score1,"Brithday_score.txt",sep="\t",row.names=T)
Brithday_score1<-read.table(file="Brithday_score.txt")

##Plot of the result
Brithday_score1$Orig_ID=factor(Brithday_score1$Orig_ID,levels=c("E10Con","E10Het","E10Dko",
                                                                "E12het","E12mut","E14Con","E14Het","E14Dko"))
beeswarm <- beeswarm(Brithday_score1$Birthday_Score ~ Brithday_score1$Orig_ID,cex=0.03,
                     data = Brithday_score1, method = 'swarm',ylim=c(-150, 100),col=brewer.pal(8, "Set1"),
                     xlab="Orig_ID",
                     ylab="Birthday_Score")
boxplot(Brithday_score1$Birthday_Score ~ Brithday_score1$Orig_ID,boxwex=0.5,col=brewer.pal(8, "Set1"),
        data = Brithday_score1, add = T,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5))

##Input differentiation weight using published data (Telley et al., 2019)
Dif_weight<-read.table(file="differentiation_weight_science.txt",sep="\t",head=T,row.names=1)
genes.use_Difweight<-rownames(Dif_weight)
genes.use_Difweight<-intersect(genes.use_Difweight,data.use@Dimnames[[1]])

####Calculate differentiation score for scData
data_list_weight_Dif<-list()
for(i in 1:8)
{
  data_list_weight_Dif[[i]]<-data.use[genes.use_Difweight, cell_list[[i]], drop = F]
}
Dif_score1<-NULL
for(i in 1:8)
{
  Dif_score<-t(Dif_weight)%*%as.matrix(data_list_weight_Dif[[i]])
  Dif_score<-t(Dif_score)
  Dif_score1<-rbind(Dif_score1,Dif_score)
}
cellid_weight<-orig_ident$orig.ident[match(rownames(Dif_score1),rownames(orig_ident))]
Dif_score1<-cbind(Dif_score1,cellid_weight)
colnames(Dif_score1)<-c("Dif_Score","Orig_ID")

##Output of differentiation score
write.table(x=Dif_score1,"Dif_score.txt",sep="\t",row.names=T)
Dif_score1<-read.table(file="Dif_score.txt")
Dif_score1$Orig_ID=factor(Dif_score1$Orig_ID,levels=c("E10Con","E10Het","E10Dko",
                                                      "E12het","E12mut","E14Con","E14Het","E14Dko"))

##Plot of the result
beeswarm <- beeswarm(Dif_score1$Dif_Score ~ Dif_score1$Orig_ID,cex=0.03,
                     data = Dif_score1, method = 'swarm',ylim=c(-120, 200),col=brewer.pal(8, "Set1"),
                     xlab="Orig_ID",
                     ylab="Differentiation_Score")
boxplot(Dif_score1$Dif_Score ~ Dif_score1$Orig_ID,boxwex=0.5,col=brewer.pal(8, "Set1"),
        data = Dif_score1, add = T,at=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5))