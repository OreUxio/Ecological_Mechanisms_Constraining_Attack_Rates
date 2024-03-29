#$Id: densograph.R 849 2007-01-19 07:45:45Z cvsrep $
library(KernSmooth)
library(lattice)
source("/home/axel/NewWeb/src/my.bkde3D.R")

bandwidth=1
range=list(c(-13.5,3),c(-9,-1))
gridsize=c(201,201)
#gridsize=c(51,51)
#gridsize=c(15,15)

ww<-function(xx,yy){
 pmax=1/sqrt(var(xx)*var(yy));
 rho<-length(xx)*pmax;
 bandwidth*rho^(-1/6)*(pmax*c(1/var(xx),1/var(yy)))^(-1/3);
}

#dxy<-matrix(c(-9,-4),ncol=2)
#estKernel<-bkde2D(dxy,bandwidth,gridsize=gridsize,range=range)

d<-scan("L1.dat")
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL1<-my.bkde2D(dxy,ww(dx,dy),gridsize=gridsize,range=range)

maxL1=expand.grid(x=estL1$x1, y=estL1$x2)[estL1$fhat==max(estL1$fhat)]
print(maxL1)
write(maxL1,"maxL1")


d<-scan("L2.dat")
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL2<-my.bkde2D(dxy,ww(dx,dy),gridsize=gridsize,range=range)

maxL2=expand.grid(x=estL2$x1, y=estL2$x2)[estL2$fhat==max(estL2$fhat)]
print(maxL2)
write(maxL2,"maxL2")


d<-scan("L3.dat")
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL3<-my.bkde2D(dxy,ww(dx,dy),gridsize=gridsize,range=range)

maxL3=expand.grid(x=estL3$x1, y=estL3$x2)[estL3$fhat==max(estL3$fhat)]
print(maxL3)
write(maxL3,"maxL3")

d<-scan("L4.dat")
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL4<-my.bkde2D(dxy,ww(dx,dy),gridsize=gridsize,range=range)

maxL4=expand.grid(x=estL4$x1, y=estL4$x2)[estL4$fhat==max(estL4$fhat)]
print(maxL4)
write(maxL4,"maxL4")



grid <- expand.grid(x=estL3$x1, y=estL3$x2,kinds=c("producers", "herbivores", "carnivores", "super carnivores") )




grid$z <- c(    
            as.vector(estL1$fhat)/max(estL1$fhat)*0.99,
            as.vector(estL2$fhat)/max(estL2$fhat)*0.99,
            as.vector(estL3$fhat)/max(estL3$fhat)*0.99,
            as.vector(estL4$fhat)/max(estL4$fhat)*0.99
        )

 pf<-function(x,y,z,at=(1:3),region=TRUE,contour=TRUE,...){
     ll<-100
     panel.levelplot(x,y,z,at=(1:ll)/ll,region=TRUE,...); 
    panel.contourplot(x,y,z,at=c(1:1)/10,contour=TRUE,region=FALSE,lty=3,...);
     panel.contourplot(x,y,z,at=c(1:5)/5,contour=TRUE,region=FALSE,...); 
     panel.abline(+0.79118,1/4,lty=2);
     panel.abline(-2.1,1/4,lty=2);
     panel.abline(-4.6,0,lty=2);
     panel.abline(-7,1);

#     panel.lines(c(-3,-3),y=c(-10,0));
#     panel.lines(c(-14,-14),y=c(-10,0));
 }


 cp<-contourplot(z~x*y|kinds, grid,
        col.regions=gray.colors(256, start = 1, end = 0.1, gamma=0.7),
        colorkey=FALSE,
        labels=FALSE,
        scales = list(x = list(alternating = 1,
			       tick.number=10,
                               tck=-1),
                      y = list(alternating = 1,
			       tick.number=5,
                               tck=-1)
                ),
        xlab=expression(paste(log[10]," body mass ",group("[",kg,"]"))),
	ylab=
	expression(paste(log[10]," biomass density ",group("[",kg/m^2,"]"))),
        layout=c(1,4),
        aspect="iso",
        as.table=TRUE,
        panel = pf,
        )


postscript(file="densograph.ps", horizontal = FALSE,
     onefile = FALSE, paper = "special", width = 4, height = 10)
plot(cp)

dev.off()
