#$Id: MB.R 630 2006-09-29 02:50:37Z cvsrep $
library(KernSmooth)
library(lattice)



bandwidth=0.6*c(1,0.25)
range=list(c(-15,-2),c(-10,0))
gridsize=c(201,201)
#gridsize=c(15,15)


dxy<-matrix(c(-9,-4),ncol=2)
estKernel<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("P1.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estP1<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("P2.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estP2<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("P3.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estP3<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("H2.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estH2<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("H3.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estH3<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)
d<-scan("C3.dat")
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estC3<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)

d<-c(scan("P3.dat"),scan("H3.dat"),scan("C3.dat"))
dx<-d[seq(1,length(d),3)]
dy<-d[seq(2,length(d),3)]
dxy<-matrix(c(dx,dy),ncol=2)
estAll3<-bkde2D(dxy,bandwidth,gridsize=gridsize,range.x=range)


grid <- expand.grid(x=estC3$x1, y=estC3$x2,lines=c("no animals","no carnivores", "no restrictions"),kinds=c("plants", "herbivores", "carnivores") )

grid$z <- c(
            as.vector(estP1$fhat)/max(estP1$fhat),
            as.vector(estP2$fhat)/max(estP2$fhat),
            as.vector(estP3$fhat)/max(estP3$fhat),
	    as.vector(estAll3$fhat)/max(estAll3$fhat),
	    as.vector(estH2$fhat)/max(estH2$fhat),
            as.vector(estH3$fhat)/max(estH3$fhat),
            as.vector(estKernel$fhat)/max(estKernel$fhat)/2,
            0*as.vector(estC3$fhat),
            as.vector(estC3$fhat)/max(estC3$fhat))

 pf<-function(x,y,z,at=(1:3),region=TRUE,contour=TRUE,...){
     ll<-50
     panel.levelplot(x,y,z,at=(1:ll)/ll,region=TRUE,...); 
     panel.contourplot(x,y,z,at=c(1:5)/5,contour=TRUE,region=FALSE,...); 
     panel.abline(-0.2,1/4);
     panel.abline(-2.9,1/4);
     panel.abline(-6,0);
     panel.abline(-log(3.16*10^4)/log(10),1);
     panel.lines(c(-3,-3),y=c(-10,0));
     panel.lines(c(-14,-14),y=c(-10,0));
 }

 cp<-contourplot(z~x*y|kinds+lines, grid,
        col.regions=gray.colors(256, start = 1, end = 0.01, gamma=0.7),
	colorkey=FALSE,
        labels=FALSE,
  	scales = list(x = list(alternating = 1,
  			       tck=-1),
  		      y = list(alternating = 1,
  			       tck=-1)
  		),
#  	xlab=NULL,ylab=NULL,
	layout=c(3,3),
	aspect="iso",
	as.table=TRUE,
	panel = pf,
  	)

pdf(file="test.pdf", width = 9, height = 9)
plot(cp)
dev.off()