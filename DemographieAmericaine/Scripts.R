# install.packages("maps")
# install.packages("FactoMineR")

library(FactoMineR)
library(maps)
library(RColorBrewer)
library(ade4)
library(corrplot)

source('ACP_TOOLS.R')

df <- read.csv('DataImi.csv', sep=';', header=TRUE)
summary(df)
head(df)

# =========================================================
#                       Cleaning Data
# =========================================================

#On met le nom des états en minuscule pour l'affichage sur une carte
df$State <- tolower(df$State)
png("HistoIm_loc.png")
hist(df$IM_Loc_01, breaks=100)
dev.off()
hist(df$IM_Loc_00, breaks=100)

df.2001 <- df[,c(1, 3, 5, 7, 9, 11, 13, 15)]

IM_Loc_01_cat <- as.factor(ifelse(df.2001$IM_Loc_01<= -6000, 'VHD',
                        ifelse(df.2001$IM_Loc_01>=-6000 & df.2001$IM_Loc_01<0 , 'HD', 
                                      ifelse(df.2001$IM_Loc_01>=0 & df.2001$IM_Loc_01<6000 , 'HA', 
                                             ifelse(df.2001$IM_Loc_01>=6000 , 'VHA', 0)))))

POP_01_cat <- as.factor(ifelse(df.2001$Pop_Tot_01< 1000000, 'LP',
                                  ifelse(df.2001$Pop_Tot_01>=1000000 & df.2001$Pop_Tot_01<3500000 , 'MP', 
                                         ifelse(df.2001$Pop_Tot_01>=3500000 & df.2001$Pop_Tot_01<8000000 , 'HP', 
                                                ifelse(df.2001$Pop_Tot_01>=8000000 , 'VHP', 0)))))

hist(df.2001$IM_Int_01, breaks=100)
IM_INT_01_cat <- as.factor(ifelse(df.2001$IM_Int_01< 10000, 'LD',
                               ifelse(df.2001$IM_Int_01>=10000 & df.2001$IM_Int_01<25000 , 'MD', 
                                      ifelse(df.2001$IM_Int_01>=25000 & df.2001$IM_Int_01<50000 , 'HD', 
                                             ifelse(df.2001$IM_Int_01>=50000 , 'VHD', 0)))))
hist(df.2001$Nb_Nai_01, breaks=100)
Nb_Nai_01_cat <- as.factor(ifelse(df.2001$Nb_Nai_01< 25000, 'LN',
                               ifelse(df.2001$Nb_Nai_01>=25000 & df.2001$Nb_Nai_01<50000 , 'MN', 
                                      ifelse(df.2001$Nb_Nai_01>=50000 & df.2001$Nb_Nai_01<110000 , 'HN', 
                                             ifelse(df.2001$Nb_Nai_01>=110000 , 'VHN', 0)))))

hist(df.2001$Nb_DC_01, breaks=100)
Nb_DC_01_cat <- as.factor(ifelse(df.2001$Nb_DC_01< 15000, 'LDC',
                               ifelse(df.2001$Nb_DC_01>=15000 & df.2001$Nb_DC_01<30000 , 'MDC', 
                                      ifelse(df.2001$Nb_DC_01>=30000 & df.2001$Nb_DC_01<60000 , 'HDC', 
                                             ifelse(df.2001$Nb_DC_01>=60000 , 'VHDC', 0)))))

hist(df.2001$Pop_inf_65_01, breaks=100)
Pop_inf_65_01_cat <- as.factor(ifelse(df.2001$Pop_inf_65_01< 1000000, 'LP',
                               ifelse(df.2001$Pop_inf_65_01>=1000000 & df.2001$Pop_inf_65_01<3500000 , 'MP', 
                                      ifelse(df.2001$Pop_inf_65_01>=3500000 & df.2001$Pop_inf_65_01<8000000 , 'HP', 
                                             ifelse(df.2001$Pop_inf_65_01>=8000000 , 'VHP', 0)))))
hist(df.2001$Pop_Sup_65_01, breaks=100)
Pop_Sup_65_01_cat <- as.factor(ifelse(df.2001$Pop_Sup_65_01< 400000, 'LP',
                                      ifelse(df.2001$Pop_Sup_65_01>=400000 & df.2001$Pop_Sup_65_01<800000 , 'MP', 
                                             ifelse(df.2001$Pop_Sup_65_01>=800000 & df.2001$Pop_Sup_65_01<1100000 , 'HP', 
                                                    ifelse(df.2001$Pop_Sup_65_01>=1100000 , 'VHP', 0)))))


df.2001['IM_LOC'] <- IM_Loc_01_cat
df.2001['POP'] <- POP_01_cat
df.2001['IM_INT'] <- IM_INT_01_cat
df.2001['NB_NAI'] <- Nb_Nai_01_cat
df.2001['NB_DC'] <- Nb_DC_01_cat
df.2001['POP_INF'] <- Pop_inf_65_01_cat
df.2001['POP_SUP'] <- Pop_Sup_65_01_cat

df.2001.AFCM <- df.2001[,c(9:15)]
  
res.mca = MCA(df.2001.AFCM)

dimdesc(res.mca)
png("EllispseDC_NAI.png")
plotellipses(res.mca,keepvar=c(5,4))
dev.off()

png("EllispseINF_SUP.png")
plotellipses(res.mca,keepvar=c(6,7))
dev.off()

png("EllispsePop_IM.png")
plotellipses(res.mca,keepvar=c(1,2,3))
dev.off()
  
colones_interet <- c(1,3, 9, 11, 13, 15,16)
df<- df[,colones_interet]

stat.graph <- 1
va.num<-which(sapply(df,is.numeric))
va.cat<-which(sapply(df,is.factor))
df.num<- df[,va.num]
df.cat<- df[,va.cat]

nb.ind <- dim(df.2001.AFCM)[1]

if (stat.graph) 
{ 
  df.num.scale <-apply(df.num,2,scale)
  PROC_BOXPLOTALL(df.num.scale, p = c(1,1), main.name = 'donn?es standardis?es')
}


# =========================================================
#  --------------- POPULATION -----------------------
# =========================================================
states.latlon <- read.csv("state_latlon.csv", header=TRUE, sep=',')

US = map("state", fill = TRUE, plot = FALSE)
Pop01 <- df[order(df$Pop_Tot_01),]
Pop00 <- df[order(df$Pop_Tot_00),]
dpt01 <- Pop01$State
dpt00 <- Pop00$State
match01 <- match.map(US, dpt01, exact=FALSE)
match00 <- match.map(US, dpt00, exact=FALSE)
blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
colors01 <- blues[match01]
colors00 <- blues[match00]
png("pop01.png")
map("state", fill=TRUE, col=colors01, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
dev.off()

png("pop00.png")
map("state", fill=TRUE, col=colors00, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
dev.off()


# =======================================================
#  --------------- Local Immigration -----------------------
# =========================================================

US = map("state", fill = TRUE, plot = FALSE)
Im.Loc.01 <- df[order(df$IM_Loc_00),]
Im.Loc.00 <- df[order(df$IM_Loc_01),]
dpt01.Im <- Im.Loc.01$State
dpt00.Im <- Im.Loc.00$State
match01.Im <- match.map(US, dpt01.Im, exact=FALSE)
match00.Im <- match.map(US, dpt00.Im, exact=FALSE)
Reds <- colorRampPalette(brewer.pal(9,"Reds"))(100)
colors01.Im <- Reds[match01.Im]
colors00.Im <- Reds[match00.Im]

png("Map_IM_Loc_01.png")
map("state", fill=TRUE, col=colors01.Im, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
min <- Im.Loc.01$IM_Loc_01[1]
max <-  Im.Loc.01$IM_Loc_01[length(Im.Loc.01$IM_Loc_01)]
legend("bottomright", legend = trunc(seq(min, max, abs(max-min)/10)), pch = 20, col = Reds[seq(0,63, 63/11)])
dev.off()

png("Map_IM_Loc_00.png")
map("state", fill=TRUE, col=colors00.Im, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
min <- Im.Loc.00$IM_Loc_00[1]
max <-  Im.Loc.00$IM_Loc_00[length(Im.Loc.00$IM_Loc_00)]
legend("bottomright", legend = trunc(seq(min, max, abs(max-min)/10)), pch = 20, col = Reds[seq(0,63, 63/11)])
dev.off()

# =======================================================
#  --------------- International Immigration -----------------------
# =========================================================

US = map("state", fill = TRUE, plot = FALSE)
Im.Int.01 <- df[order(df$IM_Int_00),]
Im.Int.00 <- df[order(df$IM_Int_01),]
dpt01.Im.Int <- Im.Int.01$State
dpt00.Im.Int <- Im.Int.00$State
match01.Im.Int <- match.map(US, dpt01.Im.Int, exact=FALSE)
match00.Im.Int <- match.map(US, dpt00.Im.Int, exact=FALSE)
Purples <- colorRampPalette(brewer.pal(9,"Purples"))(100)
colors01.Im.Int <- Purples[match01.Im.Int]
colors00.Im.Int <- Purples[match00.Im.Int]

png("Map_IM_Int_01.png")
map("state", fill=TRUE, col=colors01.Im.Int, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
min <- Im.Int.01$IM_Int_01[1]
max <-  Im.Int.01$IM_Int_01[length(Im.Int.01$IM_Int_01)]
legend("bottomright", legend = trunc(seq(min, max, abs(max-min)/10)), pch = 20, col = Purples[seq(0,63, 63/11)])
dev.off()

png("Map_IM_Int_00.png")
map("state", fill=TRUE, col=colors00.Im.Int, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
min <- Im.Int.00$IM_Int_00[1]
max <-  Im.Int.00$IM_Int_00[length(Im.Int.00$IM_Int_00)]
legend("bottomright", legend = trunc(seq(min, max, abs(max-min)/10)), pch = 20, col = Purples[seq(0,63, 63/11)])
dev.off()
# =========================================================
#  --------------- POPULATION -----------------------
# =========================================================

US = map("state", fill = TRUE, plot = FALSE)
Sup_65_01 <- df[order(df$Pop_Sup_65_01),]
Sup_65_00 <- df[order(df$Pop_sup_65_00),]
dpt01.sup <- Sup_65_01$State
dpt00.sup <- Sup_65_00$State
match01.sup <- match.map(US, dpt01.sup, exact=FALSE)
match00.sup <- match.map(US, dpt00.sup, exact=FALSE)
Greens <- colorRampPalette(brewer.pal(9,"Greens"))(100)
colors01.sup <- Greens[match01.sup]
colors00.sup <- Greens[match00.sup]

png("Map_Pop_sup_01.png")
map("state", fill=TRUE, col=colors01.sup, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
min <- Sup_65_01$Pop_Sup_65_01[1]
max <-  Sup_65_01$Pop_Sup_65_01[length(Sup_65_01$Pop_Sup_65_01)]
# legend("bottomright", legend = trunc(seq(min, max, abs(max-min)/10)), pch = 20, col = Reds[seq(0,63, 63/11)])
dev.off()


png("Map_Pop_sup_00.png")
map("state", fill=TRUE, col=colors00.sup, resolution=0)
text(states.latlon$state, x=states.latlon$longitude, y=states.latlon$latitude)
dev.off()

# =========================================================
#
# =========================================================
#-> Statistiques univari?es

stat.sum <-apply(df.2001.AFCM,2,summary) ; 
stat.sd  <-apply(df.2001.AFCM,2,sd) 
pca.stat <-rbind(stat.sum,stat.sd)


#-> Matrice des corr?lations

df01 <- df[,c(3,5,7,9,11,13,15)]
df01.cor<-cor(df01)
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

png("Corrélation.png")
corrplot(df01.cor, method="shade", shade.col=NA, tl.col="black", tl.srt=45,col=col(200), addCoef.col="black", addcolorlabel="no", order="AOE")
dev.off()

#-> Graphique en pairs
pairs(df.num)

#=====================================================================================
#                                      AFCM
#=====================================================================================

res.mca = MCA(df.2001.AFCM)


#-------------------- NB AXES ET INERTIE -----------------------
#-> inertie
inertie<-matrix(c(seq(1,21,1),res.mca$eig[,3]),ncol = 2) ; 
colnames(inertie)<-c('axe','% inertie cumul?e')
png("AxeKaiser.png")
plot(inertie[,2]~inertie[,1],type = 'b',xlab='axe',ylab='% inertie cumul?e')
dev.off()
#-> indice de Kaiser
axe.kaiser<-which(res.mca$eig[,1] >= 1)

#-> intervalle de confiance d'Anderson
Icm<-ca$eig[,1]*exp(-1.96*sqrt(2/(nb.ind-1)))  
Icp<-ca$eig[,1]*exp(+1.96*sqrt(2/(nb.ind-1)))  
axe.anderson<-as.matrix(cbind(Icm,res.mca$eig[,1],Icp),ncol = 3)

#-> boostrap total sur les axes
B = 2000 ; alpha = 0.1 ; nb.axe <-4
lam.star <-matrix(0,ncol = nb.axe,nrow = B)  
for (i in 1 : B)
{ boot.lam<- sample(size<-seq(1,nb.ind), replace = TRUE)
df.num.star <- df.2001.AFCM[boot.lam,]
mca.star<-MCA(df.num.star,graph = FALSE)
for (j in 1 : nb.axe)
{ lam.star[i,j]<-mca.star$eig[j,1]}
}  

lam.star.mean <-mean(lam.star[,1]) 
lam.star.sd <- sd(lam.star[,1])
qt<-quantile(lam.star[,1], c(alpha/2, 1 - alpha/2))  
ICm<-qt[[1]] ;
ICp<-qt[[2]]  
# histogramme
png("Bootstrap.png")
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
dev.off()
# graphique des densit?s des axes s?lectionn?s 

png("AxesDensity.png")
plot(density(lam.star[,1]),col  = 'blue',lty = 1,type ='l',xlim=c(0,1), ylim =c(0,20),
     main = 'densit? axes',cex.main = 0.6,xlab = '',ylab = '')
text(x=mean(lam.star[,1]),y = 1,label = paste('Axe ',1,sep = ''),cex =0.7,col = 'red')
for (i in 2: nb.axe)
{ 
  lines(density(lam.star[,i]),col = 'blue' ,lty = 1)
  text(x=mean(lam.star[,i]),y = 1,label = paste('Axe ',i,sep = ''),cex =0.7,col = 'red')
}
dev.off()

#-------------------- VARIABLES ET QUALITE DE REPRESENTATION -----------------------
mca.var<- res.mca$var$coord ; 
colnames(mca.var)<-paste('axes ',seq(1,5,1),sep ='')
mca.var.qlt<-res.mca$var$cos2[,c(1,2)]
mca.var.qlt<-cbind(mca.var.qlt,(apply(mca.var.qlt,1,sum))) ; 
colnames(mca.var.qlt)[3]<-'Sum qtl'

#-------------------- INDIVIDUS ET CONTRIBUTION RELATIVE -----------------------
mca.ind     <- res.mca$ind$coord ; 
colnames(mca.ind)<-paste('axes ',seq(1,5,1),sep ='')
ca.ind.ctr <- ca$row$contrib[,c(1,2)]

#-------------------- GRAPHIQUE AFC  ----------------------- -------------------------

png("MCAPlot.png")
plot(res.mca, axe=c(1,2))
dev.off()
png("VariableMCA.png")
plot.MCA(res.mca, invisible=c("ind"))
dev.off()
png("IndividusMCA.png")
plot.MCA(res.mca, invisible=c("var"))
dev.off()

df.2001.AFCM$num <- 1

ggplot(df.2001.AFCM, aes(x=df.2001.AFCM$IM_LOC, y=df.2001.AFCM$num)) +
  geom_bar(stat="identity", fill="lightblue")
ggsave("HistoPopu.png")

ggplot(df, aes(y=df$State, x=df$Pop_Tot_01)) + 
  geom_point(size=4) +
  scale_colour_brewer(palette="Set1")
ggsave("StatePopu.png")

ggplot(df, aes(x=df$IM_Loc_01, y=df$IM_Int_01)) + 
  geom_text(aes(label=df$State), vjust=-1) +
  geom_point() +
  scale_colour_brewer(palette="Set1")
ggsave("ImiLocInt.png")

ggplot(df, aes(y=df$State, x=df$IM_Loc_00)) + 
  geom_point(size=4) +
  scale_colour_brewer(palette="Set1")
ggsave("ImLOC00.png")


df$DiffPop <- df$Pop_Tot_01 - df$Pop_Tot_00
mean.diff.pop <- mean(df$DiffPop)

hist(df$DiffPop, breaks=100)
ggplot(df, aes(x=df$DiffPop)) + 
  geom_vline(xintercept = mean.diff.pop, color='red') +
  geom_bar()
ggsave("DiffPop.png")