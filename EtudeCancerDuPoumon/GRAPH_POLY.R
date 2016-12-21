rm(list = ls()) 

pop <-c(5.1,4.9,0.2,3.5,4.2,8.4,56.58,7.6,9.9,55.6)
nat <-c(11,13,16,18,12,12,13,12,12,14)

popc<-pop - mean(pop)
natc<-nat - mean(nat)

popsc <-scale(pop)
natsc <-scale(nat)

plot(nat~pop,cex.lab = 0.8, xlab='population',ylab ='natalite',pch = 19, col = 'blue',ylim = c(-5,max(nat)), xlim = c(-15,max(pop)))
points(mean(nat),mean(pop), col = 'red',cex =1,pch = 19)
text(mean(nat),mean(pop),'g',pos =  3,cex = 0.8)
points(natc~popc,pch = 19, col = 'gray',ylim = c(0,max(nat)))
points(mean(natc),mean(popc), col = 'black',cex =1,pch = 19)
text(mean(natc),mean(popc),'E(w)',pos =  3,cex = 0.8)
abline( h = mean(natc),v = mean(popc), lty = 2)

plot(natsc~popsc,cex.lab = 0.8, xlab='population ',ylab ='natalite',pch = 19, col = 'blue',ylim = c(-1.5,+2.5), xlim = c(-1,+2.5) )
abline( h = mean(natsc),v = mean(popsc), lty = 2)
points(mean(natc)~mean(popc),pch = 19, col = 'red')
