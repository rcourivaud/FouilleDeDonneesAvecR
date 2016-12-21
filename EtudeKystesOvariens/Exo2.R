#===================================================
#         premier exemple ; ACP-r
#===================================================
rm(list = ls()) 
library(FactoMineR)
library(corrplot)
library(ggplot2)

directory.down<-'//'
file.down<-'PB_II_ACP.txt' ;  

source("ACP_TOOLS.R")  

#-------------------------- formatage du data frame ---------------------------------
df<-read.table(file.down,header = TRUE)
id.names<-df[,1]

#On supprime les points les plus espacés
df <- df[-c(116, 33, 93, 67, 79, 35),]

#--------------------------------- Options -----------------------------------------------
file.nameEXL <-'PCA_ACTIVITE.xlsx'
stat.graph <- 0 # visualisation des graphiques stat desc.

#----------------------------- Pre requis avant ACP --------------------------------------

# 1.-> statistiques g?n?rales
va.num<-which(sapply(df,is.numeric))
va.cat<-which(sapply(df,is.factor))

#-> Box plot full model
df.num<- df[,va.num] ;
df.cat<- df[,va.cat] ;
nb.ind <- dim(df.num)[1]

if (stat.graph) 
{ 
  df.num.scale <-apply(df.num,2,scale)
  PROC_BOXPLOTALL(df.num.scale, p = c(1,1), main.name = 'donn?es standardis?es')
}

#-> Statistiques univari?es
stat.sum <-apply(df.num,2,summary) ; 
stat.sd  <-apply(df.num,2,sd) 
pca.stat <-rbind(stat.sum,stat.sd)

#-> Statistiques cat?gorielles: implementation dans une liste
stat.cat.list<-list() ; iter <-1
for (i in 1:length(va.num))
{
  stat.cat.list[[iter]]<-aggregate(df.num[,i],list(df.cat),summary)
  names(stat.cat.list)[iter]<-names(df.num)[i]
  iter<-iter+1
}
#-> Statistiques cat autre forlumation
# aggregate(df.num[,i], list(df.cat),mean) }

#-> Matrice des corr?lations  
pca.cor<-cor(df.num)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png("Corrélation.png")
corrplot(pca.cor, method="shade", shade.col=NA, tl.col="black", tl.srt=45,col=col(200), addCoef.col="black", addcolorlabel="no", order="AOE")
#-> Graphique en pairs
dev.off()
pairs(df.num)

#=====================================================================================
#                                      ACP
#=====================================================================================
pca<-PCA(df.num, graph = FALSE)

#-------------------- NB AXES ET INERTIE -----------------------
#-> inertie
inertie<-matrix(c(seq(1,length(va.num),1),pca$eig[,3]),ncol = 2)
colnames(inertie)<-c('axe','% inertie cumul?e')
png('AxeKaiser.png')
plot(inertie[,2]~inertie[,1],type = 'b',xlab='axe',ylab='% inertie cumul?e')
dev.off()

#-> indice de Kaiser
axe.kaiser<-which(pca$eig[,1] >= 1)

#-> intervalle de confiance d'Anderson
Icm<-pca$eig[,1]*exp(-1.96*sqrt(2/(nb.ind-1)))  
Icp<-pca$eig[,1]*exp(+1.96*sqrt(2/(nb.ind-1)))  
axe.anderson<-as.matrix(cbind(Icm,pca$eig[,1],Icp),ncol = 3)

#-> boostrap total sur les axes
B = 2000 ; alpha = 0.1 ; nb.axe <-4
lam.star <-matrix(0,ncol = nb.axe,nrow = B)  
for (i in 1 : B)
{ boot.lam<- sample(size<-seq(1,nb.ind), replace = TRUE)
df.num.star <- df.num[boot.lam,]
pca.star<-PCA(df.num.star,graph = FALSE)
for (j in 1 : nb.axe)
{ lam.star[i,j]<-pca.star$eig[j,1]}
}  

lam.star.mean <-mean(lam.star[,1]) ; 
lam.star.sd <- sd(lam.star[,1])
qt<-quantile(lam.star[,1], c(alpha/2, 1 - alpha/2)) ;  
ICm<-qt[[1]] ;ICp<-qt[[2]]  
# histogramme
hist(lam.star[,1],nclass = 50,cex.main = 0.8,freq = FALSE, cex.lab = 0.7,proba=TRUE, main = paste("f2 boostrap : nb = ", B,sep = "" ))
s<-seq(min(lam.star[,1]),max(lam.star[,1]),le=50)
# distribution normale et densit?
prov<- dnorm(s,lam.star.mean,lam.star.sd) 
lines(prov~s,col = 'red')  
lines(density(lam.star[,1]),col = 'blue',lty = 2)
# limite des intervalles de confiance et moyenne + m?diane
abline(v=mean(lam.star[,1]),col = 'red')
abline(v=median(lam.star[,1]),col = 'blue')
abline(v=ICm,col = 'red',lty = 2)
abline(v=ICp,col = 'red',lty = 2)       
# graphique des densit?s des axes s?lectionn?s 
plot(density(lam.star[,1]),col  = 'blue',lty = 1,type ='l',xlim=c(0,10), ylim =c(0,max(lam.star)),
     main = 'densit? axes',cex.main = 0.6,xlab = '',ylab = '')
text(x=mean(lam.star[,1]),y = 1,label = paste('Axe ',1,sep = ''),cex =0.7,col = 'red')
for (i in 2: nb.axe)
{ 
  lines(density(lam.star[,i]),col = 'blue' ,lty = 1)
  text(x=mean(lam.star[,i]),y = 1,label = paste('Axe ',i,sep = ''),cex =0.7,col = 'red')
}

#-------------------- VARIABLES ET QUALITE DE REPRESENTATION -----------------------
pca.var<- pca$var$coord ; colnames(pca.var)<-paste('axes ',seq(1,5,1),sep ='')
pca.var.qlt<-pca$var$cos2[,c(1,2)]
pca.var.qlt<-cbind(pca.var.qlt,(apply(pca.var.qlt,1,sum))) ; colnames(pca.var.qlt)[3]<-'Sum qtl'

#-------------------- INDIVIDUS ET CONTRIBUTION RELATIVE -----------------------
pca.ind     <- pca$ind$coord ; colnames(pca.ind)<-paste('axes ',seq(1,5,1),sep ='')
pca.ind.ctr <- pca$ind$contrib[,c(1,2)]

#-------------------- GRAPHIQUE PCA  ----------------------- -------------------------
plot.PCA(pca,axes = c(1,2),choix = 'var',cex.main = 0.8)
plot.PCA(pca,axes = c(1,2),choix = 'ind',cex.main = 0.8)

#-------------------- GRAPHIQUES PCA AVEC VARIABLES CATEGORIELLEE ET ELLIPSOIDE DE CONFIANCE -----------------------

  df.acp <-cbind(df.num,df.cat);
  names(df.acp)[(length(va.num)+1)]<-'CLA';
  res.pca<-PCA(df.acp,ncp = 2,quali.sup = (length(va.num)+1),graph = FALSE)
  concat.data   <-cbind.data.frame(df.acp[,11],res.pca$ind$coor)
  ellipse.coord <- coord.ellipse(concat.data,bary = TRUE)
  plot.PCA(res.pca,axes = c(1,2),choix = 'ind',cex.main = 0.8,title = names(df.acp)[11],
           col.quali = 'red',col.ind = 'blue',ellipse = ellipse.coord)
  
#=====================================================================================
#                            EXPORTATION EXCEL
#=====================================================================================

export2excel<-0
if(export2excel)
{
  library(excel.link)
  xls=xl.get.excel() ; 
  # test l'existance du fichier
  if(file.exists(paste(directory.down,file.nameEXL,sep = '')))
  {xl.workbook.open(paste(directory.down,file.nameEXL,sep = ''))        
  }else{xl.workbook.add() ;  wb<-xl.workbooks()}        
  # ouverture d'exc et addition d'un classeur
  sheets<-xl.sheets() ; 
  #-> Teste l'existence de la feuille
  if(length(unique(is.element(sheets,'PCA_STAT')))== 1){xl.sheet.add("PCA_STAT")}
  # addition d'une feuille
  xl.sheet.activate("PCA_STAT")
  #--> REFERENCES
  idl<- 3 ; rng=xls[["Activesheet"]]$Cells(idl,3)               ; nxt=xl.write("PCA",rng)
  idl<- idl + 1 ; rng=xls[["Activesheet"]]$Cells(idl,3)         ; nxt=xl.write('STATISTIQUES DESCRIPTIVES' ,rng)
  idl<- idl + 1 ; rng=xls[["Activesheet"]]$Cells(idl,3)         ; nxt=xl.write('statistiques globales' ,rng)
  idl<- idl + 1 ; rng=xls[["Activesheet"]]$Cells(idl,3)         ; nxt=xl.write(pca.stat,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 9 ; rng=xls[["Activesheet"]]$Cells(idl,3)         ; nxt=xl.write('Statistiques cat?gorielles' ,rng)
  for (i in 1 : length(stat.cat.list))
  {   temp <- as.matrix(stat.cat.list[[i]]) ; nl<-dim(temp)[1] 
  colnames(temp)[1] <- names(stat.cat.list)[i]
  idl<- idl + nl + 2 ; rng=xls[["Activesheet"]]$Cells(idl,3)    ; nxt=xl.write(temp,rng,row.names=TRUE,col.names=TRUE)        
  }
  
  if(length(unique(is.element(sheets,'PCA')))== 1){xl.sheet.add("PCA")}
  # addition d'une feuille
  xl.sheet.activate("PCA")
  idl <- 3       ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('Matrice des corr?lations',rng)    
  idl<- idl + 1 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.cor,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.stat,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('variables coordonn?es',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.var,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(' qualit? de r?pr?sentaion des variables',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.var.qlt,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('individus corrdonn?es',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.ind,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 16 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('Qualit? de repr?sentation individus',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.ind.qlt,rng,row.names=TRUE,col.names=TRUE)
  
  
  # xl.workbook.save(paste(directory.down,file.nameEXL,sep = ''))
}


#-------------------- calcul des distances  -----------------------

dist_eucli <- function(x1, x2, y1, y2){
  return(sqrt((x1-x2)^2 +(y1-y2)^2))
}

concat.data.trait <- concat.data
colnames(concat.data.trait) <- c("CLA", "Dim1", "Dim2")
aggdata <-aggregate(concat.data.trait[,2:3], by=list(cat = concat.data.trait$CLA),FUN=mean, na.rm=TRUE)
print(aggdata)


dist.KBEnd <- dist_eucli(concat.data.trait[,2], aggdata[1,2], concat.data.trait[,3], aggdata[1,3])
dist.KBSer <- dist_eucli(concat.data.trait[,2], aggdata[2,2], concat.data.trait[,3], aggdata[2,3])
dist.KBSM <- dist_eucli(concat.data.trait[,2], aggdata[3,2], concat.data.trait[,3], aggdata[3,3])
dist.KF <- dist_eucli(concat.data.trait[,2], aggdata[4,2], concat.data.trait[,3], aggdata[4,3])
dist.KOM <- dist_eucli(concat.data.trait[,2], aggdata[5,2], concat.data.trait[,3], aggdata[5,3])

min_dist_CAT <- function(a, b, c, d, e){
  if(a==min(a,b,c,d,e)){
    return('KBEnd')
  }
  else if(b==min(a,b,c,d,e)){
    return('KBSer')
  }
  else if(c==min(a,b,c,d,e)){
    return('KBSM')
  }
  else if(d==min(a,b,c,d,e)){
    return('KF')
  }
  else if(e==min(a,b,c,d,e)){
    return('KOM')
  }
}

result <- c()

for(i in 1:length(dist.KOM)){
  result <- c(result, min_dist_CAT(dist.KBEnd[i], dist.KBSer[i], dist.KBSM[i], dist.KF[i], dist.KOM[i]))
}
concat.data.trait.clean <- cbind(concat.data.trait, result)
png("plottest1.png")
plot(concat.data.trait.clean$Dim1, concat.data.trait.clean$Dim2, col=concat.data.trait.clean$result, type='p')
points(aggdata$Dim1, aggdata$Dim2, col= aggdata$result, cex=3, bg=c('blue', 'green'))
text(aggdata$Dim1+1, aggdata$Dim2,aggdata$cat, cex=1.5)

aggdata['result']<- aggdata$cat
all.data <- merge(concat.data.trait.clean, aggdata, by='result')
for(i in 1: length(concat.data.trait.clean$CLA)){
  x.val <- all.data[i,3]
  y.val <- all.data[i,4]
  x.cat <- all.data[i,6]
  y.cat <- all.data[i,7]
  segments(x.val,y.val , x1 = x.cat, y1 = y.cat, col = all.data[i,1], lty = par("lty"), lwd = par("lwd"))
}


diff_calc <- !as.logical(as.numeric(concat.data.trait.clean$CLA) -  as.numeric(concat.data.trait.clean$result))

summary(diff_calc)

concat.data.trait.diff <- cbind(concat.data.trait.clean, diff_calc)

#Affiche les 
points(concat.data.trait.diff[concat.data.trait.diff$diff_calc,]$Dim1, concat.data.trait.diff[concat.data.trait.diff$diff_calc,]$Dim2, col=concat.data.trait.diff[concat.data.trait.diff$diff_calc,]$result, type='p', pch=19, cex=2)
dev.off()


df$num <- 1

ggplot(df, aes(x=df$cla, y=df$num)) +
  geom_bar(stat="identity", fill="lightblue")
ggsave("NbClasse.png")

ggplot(df, aes(y=df$SURVIE, x=df$K_SCORE, colour=df$TT, shape=df$NEO)) + 
  geom_point(size=4) +
  scale_colour_brewer(palette="Set1")
ggsave("KscoreSUrvie.png")



