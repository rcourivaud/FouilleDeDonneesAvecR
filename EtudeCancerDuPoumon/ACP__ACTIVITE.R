#===================================================
#         premier exemple ; ACP-r
#===================================================
rm(list = ls()) 
library(FactoMineR)
f  <-"F" ; 
  directory<-paste(f,'://ENSEIGNEMENTS_2015//STAT_E3//SCRIPT_R',sep = '') 
  directory.down<-paste(f,'://ENSEIGNEMENTS_2015//STAT_E3//DATA//',sep = '') 
  file.down<-'EX_PCA_ACTIVITE.txt' ;  

  source("F://ENSEIGNEMENTS_2015//STAT_E3//SCRIPT_R//ACP_TOOLS.R")  

  #-------------------------- formatage du data frame ---------------------------------
  df<-read.table(paste(directory.down,file.down,sep =''),header = TRUE)
  id.names<-df[,1]
  df <- df[,2:15] ; rownames(df)<-id.names
  
  # Nommage des variables catégorielles
  df$SEX  <-factor(df$SEX, levels = c(1,2)    ,labels=c("Homme","Femme"))  
  df$ACT  <-factor(df$ACT, levels = c(1,2,9)  ,labels=c("Actif","Non.Actif","Non.précise"))  
  df$CIV  <-factor(df$CIV, levels = c(9,2,1)  ,labels=c("Non.précise","Marie","Celibataire"))  
  df$PAY  <-factor(df$PAY, levels = c(1,2,3,4),labels=c("USA", "Ouest","Est","Youg")) 
  
  #--------------------------------- Options -----------------------------------------------
  file.nameEXL <-'PCA_ACTIVITE.xlsx'
  stat.graph <- 0 # visualisation des graphiques stat desc.
  
  #----------------------------- Pre requis avant ACP --------------------------------------
  
  # 1.-> statistiques générales
  va.num<-which(sapply(df,is.numeric))
  va.cat<-which(sapply(df,is.factor))
  
  #-> Box plot full model
  df.num<- df[,va.num] ;
  df.cat<- df[,va.cat] ;
  nb.ind <- dim(df.num)[1]
  
  if (stat.graph) 
    { 
     # 2.-> Box plot catégoriel
     PROC_BOXPLOT(df[,1] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[1] )
     PROC_BOXPLOT(df[,2] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[2] )
     PROC_BOXPLOT(df[,3] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[3] )
     PROC_BOXPLOT(df[,4] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[4] )
     PROC_BOXPLOT(df[,5] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[5] )
     PROC_BOXPLOT(df[,6] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[6] )
     PROC_BOXPLOT(df[,7] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[7] )
     PROC_BOXPLOT(df[,8] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[8] )
     PROC_BOXPLOT(df[,9] , Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[9] )
     PROC_BOXPLOT(df[,10], Datacat = df[,va.cat], p = c(2,2), main.name = names(df)[10])
     
     # box plot alls variables
     PROC_BOXPLOTALL(df.num, p = c(1,1))
     df.num.scale <-apply(df.num,2,scale) ; PROC_BOXPLOTALL(df.num.scale, p = c(1,1), main.name = 'données standardisées')
     }

  #-> Statistiques univariées
  stat.sum <-apply(df.num,2,summary) ; 
  stat.sd  <-apply(df.num,2,sd) 
  pca.stat <-rbind(stat.sum,stat.sd)
  
  #-> Statistiques catégorielles: implementation dans une liste
  stat.cat.list<-list() ; iter <-1
      for (i in 1:length(va.num))
        { for (j in 1:length(va.cat))
            {
                stat.cat.list[[iter]]<-aggregate(df.num[,i],list(df.cat[,j]),summary); names(stat.cat.list)[iter]<-names(df.num)[i]
                iter<-iter+1
            }
      }
  #-> Statistiques cat autre forlumation
  # aggregate(df.num[,i], list(df.cat[,1]),mean) }
  
  #-> Matrice des corrélations  
  pca.cor<-cor(df.num)
  
  #-> Graphique en pairs
  pairs(df.num)
 
#=====================================================================================
#                                      ACP
#=====================================================================================
  pca<-PCA(df.num, graph = FALSE)
 
    #-------------------- NB AXES ET INERTIE -----------------------
    #-> inertie
    inertie<-matrix(c(seq(1,length(va.num),1),pca$eig[,3]),ncol = 2) ; colnames(inertie)<-c('axe','% inertie cumulée')
    plot(inertie[,2]~inertie[,1],type = 'b',xlab='axe',ylab='% inertie cumulée')
    
    #-> indice de Kaiser
    axe.kaiser<-which(pca$eig[,1] >= 1)
   
    #-> intervalle de confiance d'Anderson
    Icm<-pca$eig[,1]*exp(-1.96*sqrt(2/(nb.ind-1)))  
    Icp<-pca$eig[,1]*exp(+1.96*sqrt(2/(nb.ind-1)))  
    axe.anderson<-as.matrix(cbind(Icm,pca$eig[,1],Icp),ncol = 3)
    
    #-> boostrap total sur les axes
    B = 1000 ; alpha = 0.1 ; nb.axe <-4
    lam.star <-matrix(0,ncol = nb.axe,nrow = B)  
    for (i in 1 : B)
      { boot.lam<- sample(size<-seq(1,nb.ind), replace = TRUE)
        df.num.star <- df.num[boot.lam,]
        pca.star<-PCA(df.num.star,graph = FALSE)
        for (j in 1 : nb.axe)
          { lam.star[i,j]<-pca.star$eig[j,1]}
      }  
   
    lam.star.mean <-mean(lam.star[,1]) ; lam.star.sd <- sd(lam.star[,1])
    qt<-quantile(lam.star[,1], c(alpha/2, 1 - alpha/2)) ;  ICm<-qt[[1]] ;ICp<-qt[[2]]  
    # histogramme
    hist(lam.star[,1],nclass = 50,cex.main = 0.8,freq = FALSE, cex.lab = 0.7,proba=TRUE, main = paste("f2 boostrap : nb = ", B,sep = "" ))
    s<-seq(min(lam.star[,1]),max(lam.star[,1]),le=50)
    # distribution normale et densité
    prov<- dnorm(s,lam.star.mean,lam.star.sd) 
    lines(prov~s,col = 'red')  
    lines(density(lam.star[,1]),col = 'blue',lty = 2)
    # limite des intervalles de confiance et moyenne + médiane
    abline(v=mean(lam.star[,1]),col = 'red')
    abline(v=median(lam.star[,1]),col = 'blue')
    abline(v=ICm,col = 'red',lty = 2)
    abline(v=ICp,col = 'red',lty = 2)       
    # graphique des densités des axes sélectionnés 
    plot(density(lam.star[,1]),col  = 'blue',lty = 1,type ='l',xlim=c(0,10), ylim =c(0,max(lam.star)),
         main = 'densité axes',cex.main = 0.6,xlab = '',ylab = '')
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
     for (i in 1 : length(va.cat))
        { 
          df.acp <-cbind(df.num,df.cat[,i]);  names(df.acp)[(length(va.num)+1)]<-names(df.cat)[i];
          res.pca<-PCA(df.acp,ncp = 2,quali.sup = (length(va.num)+1),graph = FALSE)
          concat.data   <-cbind.data.frame(df.acp[,11],res.pca$ind$coor)
          ellipse.coord <- coord.ellipse(concat.data,bary = TRUE)
          plot.PCA(res.pca,axes = c(1,2),choix = 'ind',cex.main = 0.8,title = names(df.acp)[11],
                   col.quali = 'red',col.ind = 'blue',ellipse = ellipse.coord)
         
        }

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
  idl<- idl + 9 ; rng=xls[["Activesheet"]]$Cells(idl,3)         ; nxt=xl.write('Statistiques catégorielles' ,rng)
  for (i in 1 : length(stat.cat.list))
      {   temp <- as.matrix(stat.cat.list[[i]]) ; nl<-dim(temp)[1] 
          colnames(temp)[1] <- names(stat.cat.list)[i]
          idl<- idl + nl + 2 ; rng=xls[["Activesheet"]]$Cells(idl,3)    ; nxt=xl.write(temp,rng,row.names=TRUE,col.names=TRUE)        
      }
  
   if(length(unique(is.element(sheets,'PCA')))== 1){xl.sheet.add("PCA")}
  # addition d'une feuille
  xl.sheet.activate("PCA")
  idl <- 3       ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('Matrice des corrélations',rng)    
  idl<- idl + 1 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.cor,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.stat,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('variables coordonnées',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.var,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(' qualité de réprésentaion des variables',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.var.qlt,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 10 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('individus corrdonnées',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.ind,rng,row.names=TRUE,col.names=TRUE)
  idl<- idl + 16 ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write('Qualité de représentation individus',rng)
  idl<- idl + 1  ; rng=xls[["Activesheet"]]$Cells(idl,3)  ; nxt=xl.write(pca.ind.qlt,rng,row.names=TRUE,col.names=TRUE)
  
  
 # xl.workbook.save(paste(directory.down,file.nameEXL,sep = ''))
}


