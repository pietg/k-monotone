rm(list=ls())
library(Rcpp)
library(purrr)
sourceCpp("1_monotone.cpp")

## Poisson mixture

N1<-100
N2 <-1000
N3 <-10000
NumIt<-100

A = matrix(0,N1,1, byrow = FALSE)
B = matrix(0,N2,1, byrow = FALSE)
C = matrix(0,N3,1, byrow = FALSE)
values <- matrix(0,nrow= NumIt, ncol= 3)
colnames(values) <- c("n=100","n=1000","n=10,000")


for (j in 1: NumIt)
{	
	sim = 1001+j
	print(j)
	set.seed(sim)
	#Generate data:
	
	## Sample N1 random variates:
	Y <-rpois(N1,10)
	for (i in 1:N1)
	{
		A[i]<-rdunif(1,0,Y[i])
	}

  	monotone1 <- monotone_1(A)
	value1 <- monotone1$value
	
	## Sample N2 random variates:
	Y <-rpois(N2,10)
	for (i in 1:N2)
	{
		B[i]<-rdunif(1,0,Y[i])
	}

  	monotone2 <- monotone_1(B)
	value2 <- monotone2$value
	
	## Sample N3 random variates:
	Y <-rpois(N3,10)
	for (i in 1:N3)
	{
		C[i]<-rdunif(1,0,Y[i])
	}

  	monotone3 <- monotone_1(C)
	value3 <- monotone3$value
	
	values[j,] <- c(value1,value2,value3)
	
	  write(c(value1,value2,value3),file = "results_monotone1.txt",ncol =3,append = TRUE)
}

pdf("BoxPlot_1monotone.pdf")
boxplot(values, main="Boxplot F(10)", las=1) 
 
dev.off()
 
   D<-monotone3$F
   
   events <- 0:30
prob <- ppois(q = events, lambda = 10, lower.tail = TRUE)
   x1<-D[,1]
   y1<-D[,2]
   x2<-events
   y2<-prob

   plot(c(-1000,-1000),xlim=c(0,30), ylim=c(0,1), main= "", ylab="",xlab="")
   lines(c(x1,30),c(y1,1),lwd=1,type='s')
   lines(x2,y2,lwd=1,type='s', col="red")
  

