# Reproducibility-of-RSS-via-NP-Bootstrapping
valuating the statistical reproducibility of RSS methods using nonparametric predictive bootstrapping

set.seed(2)
N=1000
mu=100
sd=25
rho=0.40
dt=10
bt=10
m=6
v=m/2
v1=(m+2)/2
Y=rnorm(N,mu,sd)
e=rnorm(N,0,1)
X=rho*Y+e*sqrt(1-(rho^2))
df=data.frame(Y,X)
Yb=mean(Y)
Xb=mean(X)

DIFF.rss=c(); DIFF.mrss=c(); DIFF.erss=c(); DIFF.prss=c()
Arss=c(); Amrss=c(); Aerss=c(); Aprss=c()
MRSS=c(); MMRSS=c(); MERSS=c(); MPRSS=c()



for(l in 1:dt){ 				# Times Data changed
rss=c()
for(k in 1:m){		
#s=rnorm(m,mu,sd); rss=rbind(rss, sort(s))
s=df[sample(1:nrow(df), m, replace=TRUE),]; s0=s[order(s$Y),]; rss=rbind(rss, s0[,1])
}
mrss=mean(c(diag(rss)))										   # RSS
mmrss=mean(c(rss[1:v,v], rss[v1:m,v1]))							   # MRSS
merss=mean(c(rss[1:v,1], rss[v1:m,m]))							   # ERSS
mprss=mean(c(diag(rss[1:((m+1)/2),1:((m+1)/2)]), rss[3,4], rss[2,5], rss[1,6]))  # PRSS


####### NPIB-RSS   

mfrss=c(); mfmrss=c(); mferss=c(); mfprss=c()
for(i in 1:bt){				# Times Original RSS replicated 
frss=c()
for(j in 1:nrow(rss)){		

d1=sort(rss[1,])	
I11=rnorm(1, (d1[1]+d1[m])/2, (d1[m]-((d1[1]+d1[m])/2))/(((pnorm((d1[m]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I11<d1[1]){break}; I11=rnorm(1, (d1[1]+d1[m])/2, (d1[m]-((d1[1]+d1[m])/2))/(((pnorm((d1[m]-mu)/sd))^(-1))*(m/(m+1))))}
I12=runif(1,d1[1],d1[2]); I13=runif(1,d1[2],d1[3]); I14=runif(1,d1[3],d1[4]); I15=runif(1,d1[4],d1[5])
I16=rnorm(1, (d1[1]+d1[m])/2, (d1[m]-((d1[1]+d1[m])/2))/(((pnorm((d1[m]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I16>d1[m]){break}; I16=rnorm(1, (d1[1]+d1[m])/2, (d1[m]-((d1[1]+d1[m])/2))/(((pnorm((d1[m]-mu)/sd))^(-1))*(m/(m+1))))}
I1=c(I11,I12,I13,I14,I15,I16); f1=sample(I1,1); f1	

d2=sort(c(d1,f1)); #print(d2); length(d2)
I21=rnorm(1, (d2[1]+d2[m+1])/2, (d2[m+1]-((d2[1]+d2[m+1])/2))/(((pnorm((d2[m+1]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I21<d2[1]){break}; I21=rnorm(1, (d2[1]+d2[m+1])/2, (d2[m+1]-((d2[1]+d2[m+1])/2))/(((pnorm((d2[m+1]-mu)/sd))^(-1))*(m/(m+1))))}
I22=runif(1,d2[1],d2[2]); I23=runif(1,d2[2],d2[3]); I24=runif(1,d2[3],d2[4]); I25=runif(1,d2[4],d2[5]); I26=runif(1,d2[5],d2[6])
I27=rnorm(1, (d2[1]+d2[m+1])/2, (d2[m+1]-((d2[1]+d2[m+1])/2))/(((pnorm((d2[m+1]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I27>d2[m+1]){break}; I27=rnorm(1, (d2[1]+d2[m+1])/2, (d2[m+1]-((d2[1]+d2[m+1])/2))/(((pnorm((d2[m+1]-mu)/sd))^(-1))*(m/(m+1))))}
I2=c(I21,I22,I23,I24,I25,I26,I27); f2=sample(I2,1); f2	

d3=sort(c(d2,f2)); 
I31=rnorm(1, (d3[1]+d3[m+2])/2, (d3[m+2]-((d3[1]+d3[m+2])/2))/(((pnorm((d3[m+2]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I31<d3[1]){break}; I31=rnorm(1, (d3[1]+d3[m+2])/2, (d3[m+2]-((d3[1]+d3[m+2])/2))/(((pnorm((d3[m+2]-mu)/sd))^(-1))*(m/(m+1))))}
I32=runif(1,d3[1],d3[2]); I33=runif(1,d3[2],d3[3]); I34=runif(1,d3[3],d3[4]); I35=runif(1,d3[4],d3[5]); I36=runif(1,d3[5],d3[6])
I37=runif(1,d3[6],d3[7])
I38=rnorm(1, (d3[1]+d3[m+2])/2, (d3[m+2]-((d3[1]+d3[m+2])/2))/(((pnorm((d3[m+2]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I38>d3[m+2]){break}; I38=rnorm(1, (d3[1]+d3[m+2])/2, (d3[m+2]-((d3[1]+d3[m+2])/2))/(((pnorm((d3[m+2]-mu)/sd))^(-1))*(m/(m+1))))}
I3=c(I31,I32,I33,I34,I35,I36,I37,I38); f3=sample(I3,1); f3	

d4=sort(c(d3,f3))
I41=rnorm(1, (d4[1]+d4[m+3])/2, (d4[m+3]-((d4[1]+d4[m+3])/2))/(((pnorm((d4[m+3]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I41<d4[1]){break}; I41=rnorm(1, (d4[1]+d4[m+3])/2, (d4[m+3]-((d4[1]+d4[m+3])/2))/(((pnorm((d4[m+3]-mu)/sd))^(-1))*(m/(m+1))))}
I42=runif(1,d4[1],d4[2]); I43=runif(1,d4[2],d4[3]); I44=runif(1,d4[3],d4[4]); I45=runif(1,d4[4],d4[5])
I46=runif(1,d4[5],d4[6]); I47=runif(1,d4[6],d4[7]); I48=runif(1,d4[7],d4[8]) 
I49=rnorm(1, (d4[1]+d4[m+3])/2, (d4[m+3]-((d4[1]+d4[m+3])/2))/(((pnorm((d4[m+3]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I49>d4[m+3]){break}; I49=rnorm(1, (d4[1]+d4[m+3])/2, (d4[m+3]-((d4[1]+d4[m+3])/2))/(((pnorm((d4[m+3]-mu)/sd))^(-1))*(m/(m+1))))}
I4=c(I41,I42,I43,I44,I45,I46,I47,I48,I49); f4=sample(I4,1); f4	

d5=sort(c(d4,f4))
I51=rnorm(1, (d5[1]+d5[m+4])/2, (d5[m+4]-((d5[1]+d5[m+4])/2))/(((pnorm((d5[m+4]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I51<d5[1]){break}; I51=rnorm(1, (d5[1]+d5[m+4])/2, (d5[m+4]-((d5[1]+d5[m+4])/2))/(((pnorm((d5[m+4]-mu)/sd))^(-1))*(m/(m+1))))}
I52=runif(1,d5[1],d5[2]); I53=runif(1,d5[2],d5[3]); I54=runif(1,d5[3],d5[4]); I55=runif(1,d5[4],d5[5]);I56=runif(1,d5[5],d5[6])
I57=runif(1,d5[6],d5[7]); I58=runif(1,d5[7],d5[8]); I59=runif(1,d5[8],d5[9])
I510=rnorm(1, (d5[1]+d5[m+4])/2, (d5[m+4]-((d5[1]+d5[m+4])/2))/(((pnorm((d5[m+4]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I510>d5[m+4]){break}; I510=rnorm(1, (d5[1]+d5[m+4])/2, (d5[m+4]-((d5[1]+d5[m+4])/2))/(((pnorm((d5[m+4]-mu)/sd))^(-1))*(m/(m+1))))}
I5=c(I51,I52,I53,I54,I55,I56,I57,I58,I59,I510); f5=sample(I5,1); f5	

d6=sort(c(d5,f5))
I61=rnorm(1, (d6[1]+d6[m+5])/2, (d6[m+5]-((d6[1]+d6[m+5])/2))/(((pnorm((d6[m+5]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I61<d6[1]){break}; I61=rnorm(1, (d6[1]+d6[m+5])/2, (d6[m+5]-((d6[1]+d6[m+5])/2))/(((pnorm((d6[m+5]-mu)/sd))^(-1))*(m/(m+1))))}
I62=runif(1,d6[1],d6[2]); I63=runif(1,d6[2],d6[3]); I64=runif(1,d6[3],d6[4]); I65=runif(1,d6[4],d6[5]);I66=runif(1,d6[5],d6[6])
I67=runif(1,d6[6],d6[7]); I68=runif(1,d6[7],d6[8]); I69=runif(1,d6[8],d6[9]); I610=runif(1,d6[9],d6[10]); I611=runif(1,d6[10],d6[11])
I612=runif(1,d6[10],d6[11]); I613=rnorm(1, (d6[1]+d6[m+5])/2, (d6[m+5]-((d6[1]+d6[m+5])/2))/(((pnorm((d6[m+5]-mu)/sd))^(-1))*(m/(m+1))))
repeat{if(I613>d6[m+5]){break}; I613=rnorm(1, (d6[1]+d6[m+5])/2, (d6[m+5]-((d6[1]+d6[m+5])/2))/(((pnorm((d6[m+5]-mu)/sd))^(-1))*(m/(m+1))))}
I6=c(I61,I62,I63,I64,I65,I66,I67,I68,I69,I610,I611,I612); f6=sample(I6,1); f6

f=sort(c(f1,f2,f3,f4,f5,f6)); frss=rbind(frss,f)		
}
mfrss=c(mfrss, mean(c(diag(frss)))); mfmrss=c(mfmrss, mean(c(frss[1:v,v], frss[v1:m,v1]))) 
mferss=c(mferss, mean(c(frss[1:v,1], frss[v1:m,m])))
mfprss=c(mfprss, mean(c(diag(frss[1:((m+1)/2),1:((m+1)/2)]), frss[3,4], frss[2,5], frss[1,6])))
}
Arss=cbind(Arss, mfrss); Amrss=cbind(Amrss, mfmrss); Aerss=cbind(Aerss, mferss); Aprss=cbind(Aprss, mfprss)
MRSS=c(MRSS, mrss); MMRSS=c(MMRSS, mmrss); MERSS=c(MERSS, merss); MPRSS=c(MPRSS, mprss)

DIFF.rss=cbind(DIFF.rss, mrss-mfrss); DIFF.mrss=cbind(DIFF.mrss, mmrss-mfmrss)
DIFF.erss=cbind(DIFF.erss, merss-mferss); DIFF.prss=cbind(DIFF.prss, mprss-mfprss)
}
		 
############################### AD and MSD ##############################
   AD.rss=c(); MSD.rss=c()
   for(i in 1:ncol(Arss)){ AD.rss=c(AD.rss, mean(MRSS[i]-Arss[,i])); MSD.rss=c(MSD.rss, sum(((MRSS[i]-Arss[,i])^2)/bt)) }

   AD.mrss=c(); MSD.mrss=c()
   for(i in 1:ncol(Amrss)){AD.mrss=c(AD.mrss, mean(MMRSS[i]-Amrss[,i])); MSD.mrss=c(MSD.mrss, sum(((MMRSS[i]-Amrss[,i])^2)/bt)) }

   AD.erss=c(); MSD.erss=c()
   for(i in 1:ncol(Aerss)){AD.erss=c(AD.erss, mean(MERSS[i]-Aerss[,i])); MSD.erss=c(MSD.erss, sum(((MERSS[i]-Aerss[,i])^2)/bt))}

   AD.prss=c(); MSD.prss=c()
   for(i in 1:ncol(Aprss)){AD.prss=c(AD.prss, mean(MPRSS[i]-Aprss[,i])); MSD.prss=c(MSD.prss, sum(((MPRSS[i]-Aprss[,i])^2)/bt)) }

## RESULTS1
   DATA=1:dt
   Results1=cbind(DATA, AD.rss, MSD.mrss, AD.erss, MSD.prss, AD.rss, MSD.mrss, AD.erss, MSD.prss); Results1

# CDF FOR ABSOLUTE AD; library(lattice); library(latticeExtra)
  val=data.frame(RSS=MSD.rss, MRSS=MSD.mrss, ERSS=MSD.erss, PRSS=MSD.prss)
  ecdfplot(~RSS+MRSS+ERSS+PRSS, val, auto.key=list(space='right'), xlab="AD", lwd=2, ylim=c(0,1), ylab="density", main="CDF for AD")

# CDF FOR MSD ; library(lattice); library(latticeExtra)
  val=data.frame(RSS=MSD.rss, MRSS=MSD.mrss, ERSS=MSD.erss, PRSS=MSD.prss)
  ecdfplot(~RSS+MRSS+ERSS+PRSS, val, auto.key=list(space='right') ,xlab="MSD", lwd=2, ylab="density", ylim=c(0,1), main="CDF for MSD")

