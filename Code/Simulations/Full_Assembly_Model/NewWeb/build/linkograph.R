#$Id: linkograph.R 781 2006-11-22 02:31:04Z cvsrep $
library(KernSmooth)
library(lattice)

bandwidth=200*c(1,1)
range=list(c(-13,1),c(-13,1))
gridsize=c(201,201)
gridsize=c(51,51)
#gridsize=c(15,15)

#dxy<-matrix(c(-9,-4),ncol=2)
#estKernel<-bkde2D(dxy,bandwidth,gridsize=gridsize,range=range)

d<-scan(pipe("cat links.tab|cut -d' ' -f 1-2"))
dx<-d[seq(1,length(d),2)]
dy<-d[seq(2,length(d),2)]
dxy<-matrix(c(dx,dy),ncol=2)
estL1<-bkde2D(dxy,bandwidth/sqrt(length(dx)),gridsize=gridsize,range=range)

maxL1<-expand.grid(x=estL1$x1, y=estL1$x2)[estL1$fhat==max(estL1$fhat)]
#print(maxL1)
#write(maxL1,"maxL1")

masses<-scan(pipe("cat species.tab|cut -d' ' -f 1"))
print("p0")
#mass.density<-bkde(masses,bandwidth=bandwidth[1],gridsize=gridsize[1],range.x=range[1])
mass.density<-bkde(masses,bandwidth=bandwidth[1]/sqrt(length(dx)),range.x=c(-13,1),gridsize=gridsize[1])
print("p1")
print("p2")


grid <- expand.grid(x=estL1$x1, y=estL1$x2,kinds=c("connection probability") )

grid$z <- c(    
#            as.vector(estL1$fhat)/(mass.density$y %o% (0*mass.density$y+1))
#            as.vector(estL1$fhat)/((mass.density$y/seq(gridsize[1],1)) %o% mass.density$y)
            as.vector(estL1$fhat)/(mass.density$y %o% mass.density$y)
        )   

#grid$z<-as.vector((seq(gridsize[1],1) %o%  (0*mass.density$y+1)))

#grid$z<-log(grid$z)

#grid$z <- (grid$z-min(grid$z))/(max(grid$z)-min(grid$z))*2
grid$z <- (grid$z)/(max(grid$z))*0.99*4

 pf<-function(x,y,z,at=(1:3),region=TRUE,contour=TRUE,...){
     ll<-50
     panel.levelplot(x,y,z,at=(1:ll)/ll,region=TRUE,...); 
     panel.contourplot(x,y,z,at=c(1:5)/5,contour=TRUE,region=FALSE,...); 
     panel.abline(+0,1,lty=3);
#     panel.abline(-2.8,1/4,lty=2);
#     panel.abline(-5.3,0,lty=2);
#     panel.abline(-4,1);
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
			       tick.number=10,
                               tck=-1)
                ),
        xlab="log10 resource body mass [kg]",ylab="log10 consumer body mass [kg]",
        layout=c(1,1),
        aspect="iso",
        as.table=TRUE,
        panel = pf,
        )


pdf(file="densograph.pdf", width = 5, height = 12)

mass.density$y<-log10(mass.density$y)
plot(cp)
plot(mass.density, type="l")


dev.off()