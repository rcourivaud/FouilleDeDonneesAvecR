#-------------------------------------------------------
#   première approche du boostrap
#-------------------------------------------------------
# estimation de l'erreur standard de la médiance
rm(list= ls())

X<-c(25,24,23,22,41,40,59,58,47,56,45,35,38,45,29,34)
X.size <- length(X)
X.med <- median(X)

# nombre d'échantillonage et vecteur de sauvegarde des médianes
B<- 5000 ;theta.star = rep(0,B)
# graine pour le générateur de VA
set.seed(2)
# Echantillonage et calcul de la médiane
for (i in 1 : B)
    {   # générateur d'échantillons 
        X.sample      <- sample(c(1:X.size),replace = TRUE)
        X.boot        <- X[X.sample]
        theta.star[i] <- median(X.boot)
    }
# erreur standard
theta.se <- sqrt(1/(B-1)*sum((theta.star - mean(theta.star))^2))

# équivalent à ;
theta.se.r <-sd(theta.star)
theta.mean.r<-mean(theta.star)

# Intervalle de confiance
alpha =0.10 # risque de première espèce

# Estimation des intervalles de confiances
qt<-quantile(theta.star, c(alpha/2, 1 - alpha/2)) ;  ICm<-qt[[1]] ;ICp<-qt[[2]]

# Graphique
xlab <-'médiane'
hist(theta.star,nclass = 15,cex.main = 0.8,freq = FALSE,xlim=c(0,100),xlab = xlab , cex.lab = 0.7, main = paste("boostrap : nb = ", B,sep = "" ))
abline(v=ICm,col ='red',lty = 2);abline(v=ICp,col ='red',lty = 2);abline(v=theta.mean.r,col ='blue',lty = 1)


print(theta.se)
