#$Id: densograph_all.R 849 2007-01-19 07:45:45Z cvsrep $
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

d<-scan(pipe("cat L?.dat"))
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL1<-my.bkde2D(dxy,ww(dx,dy),gridsize=gridsize,range=range)

maxL1=expand.grid(x=estL1$x1, y=estL1$x2)[estL1$fhat==max(estL1$fhat)]
print(maxL1)
write(maxL1,"maxL1")



grid <- expand.grid(x=estL1$x1, y=estL1$x2,kinds=c("all species") )

grid$z <- c(    
            as.vector(estL1$fhat)/max(estL1$fhat)*0.99
        )


 pf<-function(x,y,z,at=(1:3),region=TRUE,contour=TRUE,...){
     ll<-150
     panel.levelplot(x,y,z,at=((1:ll)/ll)^1.7,region=TRUE,...); 
     panel.contourplot(x,y,z,at=c(1:1)/40,contour=TRUE,region=FALSE,lty=4,...);
     panel.contourplot(x,y,z,at=c(1:5)/5,contour=TRUE,region=FALSE,...); 
#     panel.abline(+0.79118,1/4,lty=2);
#     panel.abline(-2.1,1/4,lty=2);
#     panel.abline(-4.6,0,lty=2);
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
        layout=c(1,1),
        aspect="iso",
        as.table=TRUE,
        panel = pf,
        )


postscript(file="densograph_all.ps", horizontal = FALSE,
     onefile = FALSE, paper = "special", width = 4, height = 10)
plot(cp)

dev.off()
