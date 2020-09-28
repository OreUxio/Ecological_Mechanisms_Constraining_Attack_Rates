#ww<-function(xx,yy){
# covmat<-var(matrix(c(xx,yy),nrow = length(dx), ncol=2))
# pmax=1/sqrt(det(covmat));
# rho<-length(xx)*pmax;
# scalar.factor<-bandwidth*rho^(-1/6)*pmax^(-1/3);
# eig<-eigen(covmat);
# mat.factor=(eig$vectors) %*% (t(eig$vectors)*eig$val^(1/3));
# scalar.factor * mat.factor
#}


my.bkde2D <- function(x,bandwidth,gridsize=c(51,51),range.x,truncate=TRUE)

# Last changed: 25/08/95

{
   # Rename common variables

   n <- nrow(x)
   M <- gridsize
   h <- bandwidth
   tau <- 3.4      # For bivariate normal kernel.

   # Use same bandwidth in each
   # direction if only a single
   # bandwidth is given.

   if (length(h)==1) h <- c(h,h)

   # If range.x is not specified then
   # set it at its default value.

   if (missing(range.x))
   {
      range.x <- list(0,0)
      for (id in (1:2))
      {
         range.x[[id]] <- c(min(x[,id])-1.5*h[id],max(x[,id])+1.5*h[id])
      }
   }

   a <- c(range.x[[1]][1],range.x[[2]][1])
   b <- c(range.x[[1]][2],range.x[[2]][2])

   # Set up grid points and bin the data

   gpoints1 <- seq(a[1],b[1],length=M[1])
   gpoints2 <- seq(a[2],b[2],length=M[2])

   gcounts <- KernSmooth:::linbin2D(x,gpoints1,gpoints2)

   # Compute kernel weights

   L <- numeric()
   fac <- numeric()
   kapid <- list(0,0)
   for (id in (1:2))
   {
      L[id] <- min(floor(tau*h[id]*(M[id]-1)/(b[id]-a[id])),(M[id]-1))
      lvecid <- (0:L[id])
      fac[id] <- (b[id]-a[id])/(M[id]-1)
   }
   # Now combine weight and counts using the FFT
   # to obtain estimate

   P <- 2^(ceiling(log(M+L)/log(2))) # smallest powers of 2 >= M+L
   L1 <- L[1] ; L2 <- L[2]
   M1 <- M[1] ; M2 <- M[2]
   P1 <- P[1] ; P2 <- P[2]

   coord1=c((0:(P1/2-1)),((P1/2):(P1-1))-P1)*fac[1]
   coord2=c((0:(P2/2-1)),((P2/2):(P2-1))-P2)*fac[2]

   cat(P1,"\n")
   cat(length(coord1),"\n")
   cat(P2,"\n")
   cat(length(coord2),"\n")

   loc1=coord1%o%(coord2*0+1)
   loc2=(coord1*0+1)%o%coord2

   cat(length(x),"\n")
   ww <- 1.5*var(x)/length(x)^(1/3)
   print(sqrt(ww[1,1]))
   print(sqrt(ww[2,2]))

   iww <- solve(ww) #bandwidth quadric

   rp <- matrix(exp(-0.5*(loc1*loc1*iww[1,1]+2*loc2*loc1*iww[1,2]+loc2*loc2*iww[2,2]))*(sqrt(det(iww))/(2*pi)),P1,P2)

   sp <- matrix(0,P1,P2)
   sp[1:M1,1:M2] <- gcounts
                               # zero-padded version of "gcounts"

   rp <- fft(rp)                       # Obtain FFT's of r and s
   sp <- fft(sp)
   rp <- Re(fft(rp*sp,inverse=TRUE)/(P1*P2))[1:M1,1:M2]
                             # invert element-wise product of FFT's
                             # and truncate and normalise it

   # Ensure that rp is non-negative

   rp <- rp*matrix(as.numeric(rp>0),nrow(rp),ncol(rp))

   return(list(x1=gpoints1,x2=gpoints2,fhat=rp))
}


