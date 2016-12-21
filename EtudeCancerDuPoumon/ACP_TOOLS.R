PROC_BOXPLOT<-function(Datanum, Datacat = NULL, p = c(1,1),main.name = NULL)
{ 
  if (is.null(Datacat)){procwin(7,7);boxplot(Datanum,col = "lavender",main = main.name)}
  else
  { 
    nb.graph <-dim(Datacat)[2]
    procwin(10,10,p)
    for (i in 1: nb.graph)
    { boxplot(Datanum~Datacat[,i],col = "lavender",main = main.name)}
  }
}   

#--------------------------------------------------------------------------------
PROC_BOXPLOTALL<-function(Datanum,p = c(1,1), main.name = NULL)
{  
   nb.graph <-dim(Datanum)[2]
   procwin(10,10,p)
   boxplot(Datanum,col = "lavender", main = main.name )
}

#--------------------------------------------------------------------------------
procwin <- function(we,he,p)
{   hfig<-windows(width = we, height = he, xpos = 175, ypos = 175)
par(cex = 1, cex.main = 0.8,font.main = 3, font.lab = 3 ,bty = "u",pch = 21,bg =  "lightblue",mfrow = p)
}