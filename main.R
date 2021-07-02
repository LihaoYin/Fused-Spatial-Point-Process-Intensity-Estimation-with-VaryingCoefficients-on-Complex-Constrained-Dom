library(igraph)
library(Matrix)
library(ggplot2)
library(deldir)
library(MASS)
library(spatstat)
library(RandomFields)
library(fields)
library(cccd)
library(scales)
library(gridExtra)

setwd("/Users/yinlihao/Desktop/codes")
source("functions.R")

###### Gnerate Spatial Covariates Pattern #####
load('covariates.RData') # load the spatial covariates pattern
# The covariates are generated from Gaussian processes with mean 0 via the code below;
# set.seed(10); Phi = 0.3; np = 50
# z1.pattern = Covariates(Phi=Phi, np=np) # 2 spatil covariates
# z2.pattern = Covariates(Phi=Phi, np=np)

###### Clustered Coefficients Setting ######
ClustType=function(x,y,type=1){ 
  if(type==1){
    mclust=cut(x+y,breaks=c(-10000,1,10000),labels=FALSE);
  }
  if(type==2){
    mclust=cut(y-x,breaks=c(-10000,0,10000),labels=FALSE);
  }
  if(type==3){
    mclust=cut(x+y,breaks=c(-10000,2/3,4/3,10000),labels=FALSE);
  }
  
  return(mclust)
}

a1=c(-1,1)
a2=c(1,-1)
a0=c(-1,0,1)

#######################################
###### Scenario 2: Simulations on 2D windows ######
p=3; # 2 coefficients and 1 intercept
lx=ly=32 # length of each side of the observation window
window=matrix(c(0,lx,0,ly),2,2) # observation window
area=lx*ly

### Generate Data
set.seed(20)
coords.data=Sim_Planar(window)
colnames(coords.data)=c('lon','lat')
n.data=nrow(coords.data)
n.data # number of points

### Add Dummy Points
qnx=qny=50 # Total 60*60 dummy points added; We suggest the number of dummy points approximate to n.data
Data.PL=Add_Dummy_Planar(coords.data,window,qnx,qny,method='Poisson') # Add dummy points when using Poisson likelihood
Data.LRL=Add_Dummy_Planar(coords.data,window,qnx,qny,method='Logistic') # Add dummy points when using Logistic Regression likelihood

n.dummy=Data.PL$n.dummy; n.comb=Data.PL$n.comb
### Generate connection graphs
g <- Edgeset_Planar(Data.PL$coords,type='DT') 
## type='DT': Delaunay triangulation;
## type='KNN': K-nearestneighbor graph
## type='MST-D': Minimum spanning tree from delaunay triangulation graphs
## type='MST-K': Minimum spanning tree from KNN graphs
nlambda=20;lambda.list=10^seq(-1,1,length=nlambda)
### Using Poisson likelihood based SVCI
fit.PL=SVCI(Data.PL,g,lambda.list,method='Poisson',Adaptive=TRUE) 
Plot_Coef_Planar(beta.hat=fit.PL,beta.true=Data.PL$beta[(1+n.data):n.comb,],Data.PL$coords,window,qnx)

### Using Logistic regression likelihood based SVCI
g <- Edgeset_Planar(Data.LRL$coords,type='DT') 
n.dummy=Data.LRL$n.dummy; n.comb=Data.LRL$n.comb
fit.LRL=SVCI(Data.LRL,g,lambda.list,method='Logistic',Adaptive=TRUE)
Plot_Coef_Planar(beta.hat=fit.LRL,beta.true=Data.LRL$beta[(1+n.data):n.comb,],Data.LRL$coords,window,qnx)



###################################################
###################################################
###### Scenario2 : Simulations on Chicago Linear Netwrok ######
domain=chicago$domain # Chicago linear network
window=domain$window
vtx=domain$vertices
plot(domain)

scale_rate = 0.06 # adjust the size the linear network to control the number of spatial points
window=as.owin(scale_rate*c(window$xrange,window$yrange))
vtx$window=window;vtx$x=scale_rate*vtx$x;vtx$y=scale_rate*vtx$y
domain=linnet(vertices=vtx,edges=cbind(domain$from,domain$to))
g.chicago=graph_from_edgelist(cbind(domain$from,domain$to))


### Generate Data
set.seed(10)
data=Sim_Linnet(domain)
n.data=nrow(data$Data) # number of points
n.data

### Add Dummy Points and Construct connection graphs
Data=Add_Dummy_Linnet(data,domain,n.dummy=3000)
g <- graph_from_edgelist(cbind(Data$from,Data$to))
H=Get_H(g)

### Using Poisson likelihood based SVCI
nlambda=20;lambda.list=10^seq(-1,1,length=nlambda)
fit.PL=ADMM_Poisson(X=Data$z,y=Data$y,v=Data$v,H=H,lambda.list=lambda.list)
beta.hat=Model_Select(fit.PL,Data,method='MISE')

### illustrate the results
Plot_Coef_Linnet(beta.hat=beta.hat, Data=Data, domain, np=5) 


### Using Logistic Regession likelihood based SVCI
nlambda=20;lambda.list=10^seq(-1,1,length=nlambda)
fit.PL=ADMM_Logistic(X=Data$z,delta=Data$delta,H=H,lambda.list=lambda.list,v=c(rep(0,Data$n.data),rep(1,Data$n.dummy)))
beta.hat=Model_Select(fit.PL,Data,method='MISE')

### illustrate the results
Plot_Coef_Linnet(beta.hat=beta.hat, Data=Data, domain, np=5) 


