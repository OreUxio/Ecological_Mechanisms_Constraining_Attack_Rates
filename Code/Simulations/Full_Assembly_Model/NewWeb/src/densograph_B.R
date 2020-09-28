#$Id: densograph_B.R 849 2007-01-19 07:45:45Z cvsrep $
library(KernSmooth)
library(lattice)

nwebs=2939
bandwidth=1
rangex=c(-13,3)
rangey=c(-11,-1)
range=list(rangex,rangey)
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

bw=ww(dx,dy)

print("p0")
estL1<-bkde2D(dxy,bw,gridsize=gridsize,range=range)
print("p1")
mass.density<-bkde(dx,bandwidth=bw[1],range.x=rangex,gridsize=gridsize[1])
print("p2")


grid <- expand.grid(x=estL1$x1, y=estL1$x2,kinds=c("all species, p(log B|log M)") )

pcond<-estL1$fhat/(mass.density$y %o% ((1:gridsize[2])*0+1))

specx <- c();
specy <- c();
domix <- c();
domiy <- c();
for(ix in 1:gridsize[1]){
  s<-0;
  bsum<-0;
  for(iy in 1:gridsize[2]){
    s<-s+pcond[ix,iy]
    bsum <- bsum+(length(dx)/nwebs)*estL1$fhat[ix,iy]*10^(rangey[1]+(iy+0.5)*(rangey[2]-rangey[1])/gridsize[2])*(rangex[2]-rangex[1])/gridsize[1]
  };
  species.in.class<-(length(dx)/nwebs)*(mass.density$y[ix]*(gridsize[2]/(rangey[2]-rangey[1])));
  percentile<-1/species.in.class;
  sum.target<-s*percentile;
  lower.sum<-0;
  upper.sum<-0;
  next.lower.y<-1;
  next.upper.y<-gridsize[2];
  while(lower.sum+upper.sum<sum.target  && next.lower.y<next.upper.y){
    if(lower.sum<upper.sum){
      lower.sum<-lower.sum+pcond[ix,next.lower.y];
      next.lower.y <- next.lower.y+1;
    }else{
      upper.sum<-upper.sum+pcond[ix,next.upper.y];
      next.upper.y <- next.upper.y-1;
    }
  }
#  print(c(next.lower.y,next.upper.y));
  if(next.lower.y<next.upper.y ){
    domix <- c(domix, rangex[1]+(ix+0.5)*(rangex[2]-rangex[1])/gridsize[1]);
    domiy <- c(domiy, rangey[1]+(next.upper.y+0.5)*(rangey[2]-rangey[1])/gridsize[2]);
  }
  specx <- c(specx, rangex[1]+(ix+0.5)*(rangex[2]-rangex[1])/gridsize[1]);
  specy <- c(specy,log10(bsum));
}
print(lm(specy~specx));

print("p3")
grid$z <- c(    
            as.vector(estL1$fhat)/max(estL1$fhat)/(mass.density$y %o% ((1:gridsize[2])*0+1))
        )
print("p4")

grid$z <- grid$z*(0.99/max(grid$z))

pf<-function(x,y,z,at=(1:3),region=TRUE,contour=TRUE,...){
     ll<-100
     panel.levelplot(x,y,z,at=((1:ll)/ll),region=TRUE,...); 
#     panel.contourplot(x,y,z,at=c(1:1)/30,contour=TRUE,region=FALSE,lty=4,...); 
     panel.contourplot(x,y,z,at=c(1:5)/5,contour=TRUE,region=FALSE,...); 

#     panel.abline(+0.8,1/4,lty=2);
#     panel.abline(-2.7,1/4,lty=2);
#     panel.abline(-5.4,0,lty=2);
     panel.abline(-8,1);

#     panel.xyplot(domix,domiy,type="l",lty=6);
     panel.xyplot(specx,specy,type="l",lty=1,lwd=2,col="black");
     panel.lmline(specx,specy,lty=2,lwd=2);

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
        xlab="log10 M [kg]",ylab="log10 B/A [kg/m^2]",
        layout=c(1,1),
        aspect="iso",
        as.table=TRUE,
        panel = pf,
        )


postscript(file="densograph_B.ps", horizontal = FALSE,
     onefile = FALSE, paper = "special", width = 4, height = 10)

plot(cp)

dev.off()
