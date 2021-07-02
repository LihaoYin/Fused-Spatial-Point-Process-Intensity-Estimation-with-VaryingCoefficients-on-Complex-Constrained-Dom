
Covariates <- function(Phi, np=50){
  coords.grid = expand.grid(seq(1/(2*np),1-1/(2*np),1/np),seq(1/(2*np),1-1/(2*np),1/np))
  dist.grid = fields::rdist(coords.grid)
  corr = exp(-dist.grid/Phi)
  z = mvrnorm(n=1,rep(0,dim(corr)[1]),corr)
  z = matrix(z,np,np)
  return(z)
}

Map_Beta <- function(x, y, ClustType, a, type){
  mclust=ClustType(x,y,type=type)
  beta=a[mclust]
  return(beta)
}

Map_Covariate <- function(x,y,z){
  np=nrow(z)
  idy=cut(y,breaks=seq(0,1,length=(np+1)),labels=FALSE)
  idx=cut(x,breaks=seq(0,1,length=(np+1)),labels=FALSE)
  return(z[cbind(idx,idy)])
}

Sim_Planar <- function(window){
  n = 55000; np = nrow(z1.pattern)
  lx=window[2,1];ly=window[2,2]
  lon=runif(n)*lx
  lat=runif(n)*ly
  area=lx*ly
  
  beta1=Map_Beta(lon/lx,lat/ly,ClustType,a1,type=1)
  beta2=Map_Beta(lon/lx,lat/ly,ClustType,a2,type=2)
  beta0=Map_Beta(lon/lx,lat/ly,ClustType,a0,type=3)
  z1=Map_Covariate(lon/lx,lat/ly,z1.pattern)
  z2=Map_Covariate(lon/lx,lat/ly,z2.pattern)
  intensity=exp(z1*beta1+z2*beta2+beta0)
  
  while(max(intensity/(n/area)) > 1){
    n=ceiling(1.5*max(intensity*area)+1)
    lon=runif(n)*lx
    lat=runif(n)*ly
    beta1=Map_Beta(lon/lx,lat/ly,ClustType,a1,type=1)
    beta2=Map_Beta(lon/lx,lat/ly,ClustType,a2,type=2)
    beta0=Map_Beta(lon/lx,lat/ly,ClustType,a0,type=3)
    z1=Map_Covariate(lon/lx,lat/ly,z1.pattern)
    z2=Map_Covariate(lon/lx,lat/ly,z2.pattern)
    intensity=exp(z1*beta1+z2*beta2+beta0)
  }
  
  intensity[intensity / (n/area) < runif(n)] <- NA
  lon <- lon[-which(is.na(intensity))]
  lat <- lat[-which(is.na(intensity))]
  return(cbind(lon,lat))
}


Add_Dummy_Planar<-function(coords.data, window, qnx=NULL, qny=NULL, n.dummy=NULL, method='Poisson'){
  if(is.null(qny)){
    qny=qnx
  }
  lx=window[2,1];ly=window[2,2];area=lx*ly
  if(is.null(n.dummy)){
    n.dummy=qnx*qny;
  }
  n.data=nrow(coords.data);n.comb=n.dummy+n.data;
  pp=as.ppp(coords.data, as.vector(window))
  pixel=t(quadratcount(pp, nx=qnx, ny=qny))[,qny:1];
  v.dummy=(area/n.dummy)/(c(pixel)+1);
  idy=cut(coords.data[,2],breaks=seq(window[1,1],window[2,1],length=(qny+1)),labels=FALSE)
  idx=cut(coords.data[,1],breaks=seq(window[1,1],window[2,1],length=(qnx+1)),labels=FALSE)
  v.data=(area/n.dummy)/(pixel[cbind(idx,idy)]+1);
  v.comb=c(v.data,v.dummy)
  y.comb=c(1/v.data,rep(0,n.dummy))
  coords.dummy=expand.grid(window[1,1]+seq(window[2,1]/(2*qnx),window[2,1]-(window[2,1]/(2*qnx)),length=qnx),
                           window[1,2]+seq(window[2,2]/(2*qny),window[2,2]-(window[2,2]/(2*qny)),length=qny));

  colnames(coords.dummy)=c('lon','lat')
  coords.comb=rbind(coords.data,coords.dummy) 
  
  beta1=Map_Beta(coords.comb[,1]/lx,coords.comb[,2]/ly,ClustType,a1,type=1)
  beta2=Map_Beta(coords.comb[,1]/lx,coords.comb[,2]/ly,ClustType,a2,type=2)
  beta0=Map_Beta(coords.comb[,1]/lx,coords.comb[,2]/ly,ClustType,a0,type=3)
  z1=Map_Covariate(coords.comb[,1]/lx,coords.comb[,2]/ly,z1.pattern)
  z2=Map_Covariate(coords.comb[,1]/lx,coords.comb[,2]/ly,z2.pattern)

  beta.comb=cbind(beta1,beta2,beta0)
  z0=rep(1,n.comb)
  z=cbind(z1,z2,z0)
  return(list(n.data=n.data,n.dummy=n.dummy,n.comb=n.comb,
              coords=coords.comb,z=z,beta=beta.comb,
              v=v.comb,y=y.comb,delta=n.dummy/area))
}


Edgeset_Planar <- function(coords,type='KNN',para=5){
  ### Generate Edgeset of the Data on Planar Windows
  ## type='DT': Delaunay triangulation;
  ## type='KNN': K-nearestneighbor graph
  ## type='MST-D': Minimum spanning tree from delaunay triangulation graphs
  ## type='MST-K': Minimum spanning tree from KNN graphs
  n=nrow(coords)
  if(type=='DT'){
    g.D = deldir(x=coords[,1],y=coords[,2])
    w=rdist.vec(cbind(g.D$delsgs$x1,g.D$delsgs$y1),cbind(g.D$delsgs$x2,g.D$delsgs$y2))
    g.D = graph_from_edgelist(cbind(g.D$delsgs$ind1,g.D$delsgs$ind2))
    E(g.D)$weight <- w
    return(g.D)
  }
  
  if(type=='MST-D'){
    g.D = deldir(coords)
    w=rdist.vec(cbind(g.D$delsgs$x1,g.D$delsgs$y1),cbind(g.D$delsgs$x2,g.D$delsgs$y2))
    g.D = graph_from_edgelist(cbind(g.D$delsgs$ind1,g.D$delsgs$ind2))
    E(g.D)$weight <- w
    g.MST <- mst(g.D, weights = w)
    return(g.MST)
  }
  
  g.KNN = nng(coords,k=para)
  g.KNN.dataframe <- igraph::as_data_frame(g.KNN)
  from <- g.KNN.dataframe$from
  to <- g.KNN.dataframe$to
  w <- rep(0,length(from))
  for(i in 1:length(from)){
    w[i]=fields::rdist(coords[c(from[i],to[i]),])[1,2]
  }
  E(g.KNN)$weight <- w
  
  if(type=='KNN'){
    return(g.KNN)
  }
  
  if(type=='MST-K'){
    g.MST <- mst(g.KNN, weights = w)
    return(g.MST)
  }
}

Get_H <- function(g){
  np <- length(V(g))
  g.dataframe <- igraph::as_data_frame(g)
  from <- g.dataframe$from
  to <- g.dataframe$to

  nv <- length(from)
  idi <- c(1:nv,1:nv)
  idj <- c(from,to)
  x <- c(rep(1,nv),rep(-1,nv))
  H <- sparseMatrix(i=idi, j=idj, x=x, dims=c(nv,np))
  
  return(H)
}

lasso_shrinkage<-function(a, kappa, w){
  
  sig=sign(a)
  y=abs(a)-kappa*w
  y[which(y<0)]=0
  y=sig*y
  return(y);
}


ADMM_Poisson<-function(X, y, H, v, lambda.list,
                       reltol=1e-04, abstol=1e-02/nrow(X), maxiter=100, rho=1,
                       w=NULL, init=NULL, beta0=NULL){
  # 1. get parameters
  n=dim(X)[1];p=dim(X)[2];m=dim(H)[1];B=m/n;
  if(is.null(w)){
    w=1
  }
  
  # 2. set ready
  if(is.null(init)){
    beta=matrix(rnorm(n*p)/10,n,p)
  }else{
    beta=init
  }
  
  # 3. precompute static variables for x-update and factorization
  XX = X^2
  HtH = t(H)%*%H
  D = Diagonal(n)+rho*HtH
  U = Cholesky(D, super=TRUE)
  
  # 4. iteration
  sqrtn = sqrt(p*n);
  lambda.list.sort=sort(lambda.list,decreasing=TRUE)
  nlambda=length(lambda.list)
  output=list()
  output$lambda=lambda.list.sort
  output$x=matrix(0,n*p,nlambda)
  output$k=rep(0,nlambda)
  
  for(j in 1:length(lambda.list.sort)){
    lambda=lambda.list.sort[j]
    z=H%*%beta
    u=matrix(0,m,p)
    Hbeta = z
    
    h_r_norm=rep(0,maxiter)
    h_s_norm=rep(0,maxiter)
    h_eps_pri=rep(0,maxiter)
    h_eps_dual=rep(0,maxiter)
    
    for(k in 1:maxiter){
      L = max(v*exp(rowSums(X*beta)))

      # 4-1. update 'x'
      for(i in 1:p){
        expXb = exp(rowSums(X*beta))
        q = beta[,i]+v*X[,i]*(y-expXb)/L+rho*t(H)%*%(z[,i]-u[,i])
        beta[,i]=as.vector(solve(D,q))
      }
      
      # 4-2. update 'z'
      Hbeta = H%*%beta
      z_old = z;
      z = lasso_shrinkage(Hbeta + u, lambda/rho/B/L, w);
      
      # 4-3. update 'u'
      u = u + Hbeta - z;
      
      # 4-4. dianostics, reporting
      beta_norm = mean(abs(beta))
      z_norm = mean(abs(-z))
      h_r_norm[k] = mean(abs(Hbeta-z));
      h_s_norm[k] = mean(abs(-rho*(z-z_old)));
      if (beta_norm>z_norm){
        h_eps_pri[k] = sqrtn*abstol + reltol*beta_norm;
      } else {
        h_eps_pri[k] = sqrtn*abstol + reltol*z_norm;
      }
      h_eps_dual[k] = sqrtn*abstol + reltol*mean(abs(rho*u));
      
      # 4-5. termination
      if ((h_r_norm[k] < h_eps_pri[k])&(h_s_norm[k]<h_eps_dual[k])){
        break;
      }
      
      if(!is.null(beta0)){
        print(k)
        print(mean((beta-as.vector(beta0))^2))
      }
    }
    
    output$x[,j] = as.vector(beta); 
    output$df[j]=length(unique(beta))
    output$k[j]=k
  }
  
  # 5. report results
  return(output);  
}


ADMM_Logistic<-function(X, delta, H, v, lambda.list, L=1,
                        reltol=1e-04, abstol=1e-02/nrow(X), maxiter=100, rho=1,
                        w=NULL, init=NULL, beta0=NULL){
  # 1. get parameters
  n=dim(X)[1];p=dim(X)[2];m=dim(H)[1];B=m/n;
  if(is.null(w)){
    w=1
  }
  
  # 2. set ready
  if(is.null(init)){
    beta=matrix(rnorm(n*p)/10,n,p)
  }else{
    beta=init
  }
  
  # 3. precompute static variables for x-update and factorization
  XX = X^2
  HtH = t(H)%*%H
  D = Diagonal(n)+rho*HtH
  U = Cholesky(D, super=TRUE)
  
  # 4. iteration
  sqrtn = sqrt(p*n);
  lambda.list.sort=sort(lambda.list,decreasing=TRUE)
  nlambda=length(lambda.list)
  output=list()
  output$lambda=lambda.list.sort
  output$x=matrix(0,n*p,nlambda)
  output$k=rep(0,nlambda)
  
  for(j in 1:length(lambda.list.sort)){
    lambda=lambda.list.sort[j]
    z=H%*%beta
    u=matrix(0,m,p)
    Hbeta = z

    h_r_norm=rep(0,maxiter)
    h_s_norm=rep(0,maxiter)
    h_eps_pri=rep(0,maxiter)
    h_eps_dual=rep(0,maxiter)
    
    for(k in 1:maxiter){
      # 4-1. update 'x'
      for(i in 1:p){
        expXb = exp(rowSums(X*beta))
        q = delta/(expXb+delta)-v
        q =beta[,i]+X[,i]*q/L+rho*t(H)%*%(z[,i]-u[,i])
        beta[,i]=as.vector(solve(U,q))
      }
      
      # 4-2. update 'z'
      Hbeta = H%*%beta
      z_old = z;
      z = lasso_shrinkage(Hbeta + u, lambda/rho/B/L, w);
      
      # 4-3. update 'u'
      u = u + Hbeta - z;
      
      # 4-4. dianostics, reporting
      beta_norm = mean(abs(beta))
      z_norm = mean(abs(-z))
      h_r_norm[k] = mean(abs(Hbeta-z));
      h_s_norm[k] = mean(abs(-rho*(z-z_old)));
      if (beta_norm>z_norm){
        h_eps_pri[k] = sqrtn*abstol + reltol*beta_norm;
      } else {
        h_eps_pri[k] = sqrtn*abstol + reltol*z_norm;
      }
      h_eps_dual[k] = sqrtn*abstol + reltol*mean(abs(rho*u));
      
      # 4-5. termination
      if ((h_r_norm[k] < h_eps_pri[k])&(h_s_norm[k]<h_eps_dual[k])){
        break;
      }
      
      if(!is.null(beta0)){
        print(k)
        print(mean((beta-as.vector(beta0))^2))
      }
    }
    
    output$x[,j] = as.vector(beta); 
    output$df[j]=length(unique(beta))
    output$k[j]=k
  }
  
  # 5. report results
  return(output);  
}

SVCI <- function(Data, g, lambda, method='Poisson', Adaptive=TRUE, weight=NULL){
  H <- Get_H(g); delta=Data$delta
  if(isFALSE(1:10)){weight=1}
  if(is.null(weight)){
    weight=exp(-E(g)$weight)
  }
  if(method=='Poisson'){
    fit=ADMM_Poisson(X=Data$z,y=Data$y,v=Data$v,H=H,w=weight,
                     lambda.list=lambda.list)
  }else if(method=='Logistic'){
    fit=ADMM_Logistic(X=Data$z,delta=delta,H=H,v=c(rep(0,Data$n.data),rep(1,Data$n.dummy)),
                      lambda.list = lambda.list,L=0.5,w=weight)
    
  }
  n.comb=Data$n.comb; p=ncol(Data$z)
  BIC1=c();BIC2=c()
  for (i in 1:length(fit$lambda)){
    tmp=matrix(fit$x[,i],n.comb,p)
    BIC1=c(BIC1,Cal_BIC(tmp,Data,rd=1))
    BIC2=c(BIC2,Cal_BIC(tmp,Data,rd=2))
  }
  it=ceiling((which(BIC1==min(BIC1))+which(BIC2==min(BIC2)))/2);
  beta.hat=matrix(fit$x[,it],n.comb,p)
  
  if(Adaptive){
    w_ad=matrix(0,dim(H)[1],0)
    for(i in 1:p){
      #tmp=scales::rescale(as.vector(abs(H%*%beta.hat[,i])),to=c(0.1,1))^2
      w_ad=cbind(w_ad,weight/as.vector(abs(H%*%beta.hat[,i]))^2)
    }
    return(SVCI(Data,g,lambda,method=method,Adaptive=FALSE,weight=w_ad))
  }else{
    return(beta.hat)
  }
}

Deviance <- function(beta, y, X, v, n.data){
  mu = exp(rowSums(X*beta))
  logL = sum(v*(y-mu))-sum((v*y*log(y/mu))[1:n.data])
  return(logL)
}

Cal_BIC <- function(beta, Data,rd=2){
  tmp=round(beta,digits=rd)
  p=ncol(Data$z); k=0;
  for(j in 1:p){
    k=k+length(unique(tmp[,j]))
  }
  L=Deviance(beta,Data$y,Data$z,Data$v,Data$n.data)
  BIC=-L+k*log(Data$n.comb)
  return(BIC)
}

Plot_Coef_Planar <- function(beta.hat, beta.true, coords, window, np){
  p=ncol(beta.hat); q=list()
  coords.dummy=expand.grid(window[1,1]+seq(window[2,1]/(2*np),window[2,1]-(window[2,1]/(2*np)),length=np),
                           window[1,2]+seq(window[2,2]/(2*np),window[2,2]-(window[2,2]/(2*np)),length=np));
  idx=cut(coords[,1],breaks=seq(window[1,1],window[2,1],length=(np+1)),labels=FALSE)
  idy=cut(coords[,2],breaks=seq(window[1,2],window[2,2],length=(np+1)),labels=FALSE)
  for(j in 1:p){
    output=data.frame('beta'=beta.true[,j],'lon'=coords.dummy[,1],'lat'=coords.dummy[,2])
    q[[j]]=ggplot(output, aes(lon, lat)) + ggtitle(paste('true beta',as.character(j),' ')) +
      geom_raster(aes(fill = beta)) +
      scale_fill_gradientn(colours = rainbow(2),limits=c(-2,2),name=NULL)+ 
      theme(legend.position = "none")
    
    output=matrix(0,np,np)
    for(s in 1:np){
      for(t in 1:np){
        output[s,t]=mean(beta.hat[which((idx==s)&(idy==t)),j])
      }
    }
    output=data.frame('beta'=as.vector(output),'lon'=coords.dummy[,1],'lat'=coords.dummy[,2])
    q[[j+p]]=ggplot(output, aes(lon, lat)) + ggtitle(paste('estimate',as.character(j),' ')) +
      geom_raster(aes(fill = beta)) +
      scale_fill_gradientn(colours = rainbow(2),limits=c(-2,2),name=NULL)+ 
      theme(legend.position = "none")
  }
  grid.arrange(grobs=q,nrow=2, ncol=3)
}

Plot_Net <- function(domain){
  vtx = domain$vertices
  edges <- data.frame(vtx$x[domain$from],vtx$y[domain$from],vtx$x[domain$to],vtx$y[domain$to])
  colnames(edges) <- c("lon","lat","X2","Y2")
  ggplot() + 
    geom_segment(aes(x=lon, y=lat, xend = X2, yend = Y2), data=edges, size = 0.5, colour="black") 
}


ClustType_Linnet<-function(ends,window,type=1){
  a=c(-1,1)
  X1=window$xrange[1]+0.35*(window$xrange[2]-window$xrange[1])
  X2=window$xrange[1]+0.65*(window$xrange[2]-window$xrange[1])
  y1=window$yrange[1]+0.35*(window$yrange[2]-window$yrange[1])
  y2=window$yrange[1]+0.65*(window$yrange[2]-window$yrange[1])
  if(type==3){
    if((ends$x0 <= X2) & (ends$x0 >= X1) & (ends$y0 <= y2) & (ends$y0 >= y1)&(ends$x1 <= X2) & (ends$x1 >= X1) & (ends$y1 <= y2) & (ends$y1 >= y1)){
      return(a[2])
    }else{
      return(a[1])
    }
  }
  if(type==2){
    if((ends$x0<= X1)&(ends$x1<= X1)){
      return(a[2])
    }else{
      return(a[1])
    }
  }
  if(type==1){
    if((ends$y0<= y1)&(ends$y1<= y1)){
      return(a[2])
    }else{
      return(a[1])
    }
  }
}

Sim_Linnet <- function(domain){
  nl=domain$lines$n;dimyx=nrow(z1.pattern);N_PER_LEN=1
  label=matrix(0,nl,p);
  len=sqrt((domain$lines$ends$x1-domain$lines$ends$x0)^2+
             (domain$lines$ends$y1-domain$lines$ends$y0)^2)
  
  x=c();y=c();seg=c();tp=c()
  for(i in 1:nl){
    intensity=1
    n=max(ceiling(N_PER_LEN*len[i]),10)
    tps=runif(n)
    s1=domain$lines$ends$x0[i]+tps*(domain$lines$ends$x1[i]-domain$lines$ends$x0[i])
    s2=domain$lines$ends$y0[i]+tps*(domain$lines$ends$y1[i]-domain$lines$ends$y0[i])
    idy=cut(s2,breaks=seq(domain$window$yrange[1],domain$window$yrange[2],length=(dimyx+1)),labels=FALSE)
    idx=cut(s1,breaks=seq(domain$window$xrange[1],domain$window$xrange[2],length=(dimyx+1)),labels=FALSE)
    
    label[i,1]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=1)
    label[i,2]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=2)
    label[i,3]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=3)

    intensity=label[i,3]+z1.pattern[cbind(idx,idy)]*label[i,1]+z2.pattern[cbind(idx,idy)]*label[i,2]
    intensity=exp(intensity)
    
    while((max(intensity/(n/len[i]))>=1)){
      n=ceiling(n*(max(intensity/(n/len[i]))+1))
      tps=runif(n)
      s1=domain$lines$ends$x0[i]+tps*(domain$lines$ends$x1[i]-domain$lines$ends$x0[i])
      s2=domain$lines$ends$y0[i]+tps*(domain$lines$ends$y1[i]-domain$lines$ends$y0[i])
      idy=cut(s2,breaks=seq(domain$window$yrange[1],domain$window$yrange[2],length=(dimyx+1)),labels=FALSE)
      idx=cut(s1,breaks=seq(domain$window$xrange[1],domain$window$xrange[2],length=(dimyx+1)),labels=FALSE)
      
      label[i,1]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=1)
      label[i,2]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=2)
      label[i,3]=ClustType_Linnet(domain$lines$ends[i,],domain$window,type=3)
      
      intensity=label[i,3]+z1.pattern[cbind(idx,idy)]*label[i,1]+z2.pattern[cbind(idx,idy)]*label[i,2]
      intensity=exp(intensity)
    }
    
    tps[intensity/ (n/len[i]) < runif(n)] <- NA
    s1 <- s1[-which(is.na(tps))]
    s2 <- s2[-which(is.na(tps))]
    tps <- tps[-which(is.na(tps))]
    x=c(x,s1);y=c(y,s2);tp=c(tp,tps);seg=c(seg,rep(i,length(s1)))
  }
  return(list(Data=data.frame(x=x,y=y,seg=seg,tp=tp),labels=label))
}

Add_Dummy_Linnet<-function(Data,linnet,n.dummy){
  x=c();y=c();seg=c();tp=c();
  nl=linnet$lines$n;nd=nrow(Data$Data);nv=linnet$vertices$n;
  v.data=rep(0,nd);y.data=rep(0,nd);v.dummy=c();delta=c();
  from=c();to=c()
  vtx.list=vector('list',nv)
  np=round(n.dummy/nl)
  len=sqrt((linnet$lines$ends$x1-linnet$lines$ends$x0)^2+
             (linnet$lines$ends$y1-linnet$lines$ends$y0)^2)
  id.dummy=nd
  for(i in 1:nl){
    l.seg=len[i]/np;id.dummy=(id.dummy[length(id.dummy)]+1):(id.dummy[length(id.dummy)]+np)
    tp.dummy=seq(1/(2*np),1-1/(2*np),1/np)
    id.data=which(Data$Data$seg==i)
    tp.data=Data$Data$tp[id.data]
    n.count=hist(x=tp.data,breaks=seq(0,1,length.out=np+1),plot=FALSE)$counts
    v.dummy=c(v.dummy,l.seg/(n.count+1))
    idx=cut(tp.data,breaks=seq(0,1,length=(np+1)),labels=FALSE)
    v.data[id.data]=l.seg/(n.count[idx]+1)
    y.data[id.data]=1/v.data[id.data]
    x=c(x,linnet$lines$ends$x0[i]+tp.dummy*(linnet$lines$ends$x1[i]-linnet$lines$ends$x0[i]))
    y=c(y,linnet$lines$ends$y0[i]+tp.dummy*(linnet$lines$ends$y1[i]-linnet$lines$ends$y0[i]))
    seg=c(seg,rep(i,np))
    tp=c(tp,tp.dummy)
    
    tp.comb=c(tp.data,tp.dummy)
    id.comb=c(id.data,id.dummy)
    n.comb=length(tp.comb)
    odr=order(tp.comb)
    id.odr=id.comb[odr]
    from=c(from,id.odr[-n.comb])
    to=c(to,id.odr[-1])
    
    v1=linnet$from[i];v2=linnet$to[i]
    vtx.list[[v1]]=c(vtx.list[[v1]],id.odr[1])
    vtx.list[[v2]]=c(vtx.list[[v2]],id.odr[n.comb])
    
    delta=c(delta,np/len[i])
  }
  
  for(i in 1:nv){
    cnd=vtx.list[[i]]
    if(length(cnd)>1){
      g=make_full_graph(length(cnd),directed=FALSE)
      g=igraph::as_data_frame(g)
      from=c(from,cnd[g$from])
      to=c(to,cnd[g$to])
    }
  }
  
  output=list()
  output$n.data=nrow(Data$Data);output$n.dummy=length(x);
  output$n.comb=output$n.data+output$n.dummy
  output$coords=cbind(c(Data$Data$x,x),c(Data$Data$y,y))
  output$coords.dummy=data.frame(x=x,y=y,seg=seg,tp=tp)
  output$coords.data=Data$Data
  output$y=c(y.data,rep(0,length(v.dummy)))
  output$v=c(v.data,v.dummy)
  output$from=from
  output$to=to
  
  labels=Data$labels; output$labels=labels
  beta=matrix(0,output$n.comb,0)
  for(j in 1:ncol(labels)){
    beta=cbind(beta,labels[c(Data$Data$seg,seg),j])
  }
  output$beta=beta
  
  dimyx=nrow(z1.pattern)
  idy=cut(output$coords[,2],breaks=seq(linnet$window$yrange[1],linnet$window$yrange[2],length=(dimyx+1)),labels=FALSE)
  idx=cut(output$coords[,1],breaks=seq(linnet$window$xrange[1],linnet$window$xrange[2],length=(dimyx+1)),labels=FALSE)
  output$z=cbind(z1.pattern[cbind(idx,idy)],z2.pattern[cbind(idx,idy)],z0=rep(1,output$n.comb))
  
  output$delta=delta[c(output$coords.data$seg,output$coords.dummy$seg)]
  return(output)
}

Plot_Coef_Linnet <- function(beta.hat, Data, domain, np=5){
  coords=rbind(Data$coords.data,Data$coords.dummy);labels=Data$labels
  p=ncol(beta.hat); q=vector(mode='list',length=2*p)
  x=domain$vertices$x;y=domain$vertices$y;
  nl=domain$lines$n;nv=domain$vertices$n;
  from=c();to=c();seg=c()
  beta.es=matrix(0,0,p);beta.true=matrix(0,0,p)
  
  len=sqrt((domain$lines$ends$x1-domain$lines$ends$x0)^2+
             (domain$lines$ends$y1-domain$lines$ends$y0)^2)
  
  id.dummy=(nv-np+2):nv
  for(i in 1:nl){
    coords.tmp=coords[which(coords$seg==i),]
    beta.tmp=beta.hat[which(coords$seg==i),]
    l=len[i]/np;id.dummy=id.dummy+np-1
    tp.dummy=seq(1/np,1-1/np,1/np)
    
    x=c(x,domain$lines$ends$x0[i]+tp.dummy*(domain$lines$ends$x1[i]-domain$lines$ends$x0[i]))
    y=c(y,domain$lines$ends$y0[i]+tp.dummy*(domain$lines$ends$y1[i]-domain$lines$ends$y0[i]))
    
    v1=domain$from[i];v2=domain$to[i]
    from=c(from,c(v1,id.dummy))
    to=c(to,c(id.dummy,v2))
    seg=c(seg,rep(i,np))
    
    for(k in 1:np){
      tmp=beta.tmp[which((coords.tmp$tp<(k/np))&(coords.tmp$tp>((k-1)/np))),]
      tmp=matrix(tmp,length(tmp)/3,3)
      beta.es=rbind(beta.es,colMeans(tmp))
      
      beta.true=rbind(beta.true,labels[i,])
    }
  }
  
  for(k in 1:p){
    edges <- data.frame(x[from],y[from],x[to],y[to],beta.true[,k])
    colnames(edges) <- c("lon","lat","X2","Y2","beta")
    q[[k]]=ggplot() + ggtitle(paste('true beta',as.character(k),' '))+
      geom_segment(aes(x=lon, y=lat, xend = X2, yend = Y2, colour=beta),
                   data=edges, size = 0.8)+
      scale_color_gradientn(colours = rainbow(2), guide = "colourbar",limits=c(-2,2),name=NULL,values=c(0,.3,1))+ 
      theme(legend.position = "none",panel.background=element_rect(fill="white", colour = "black"))
    
    edges <- data.frame(x[from],y[from],x[to],y[to],beta.es[,k])
    colnames(edges) <- c("lon","lat","X2","Y2","beta")
    q[[k+p]]=ggplot() + ggtitle(paste('estimate',as.character(k),' '))+
      geom_segment(aes(x=lon, y=lat, xend = X2, yend = Y2, colour=beta),
                   data=edges, size = 0.8)+
      scale_color_gradientn(colours = rainbow(2), guide = "colourbar",limits=c(-2,2),name=NULL,values=c(0,.3,1))+ 
      theme(legend.position = "none",panel.background=element_rect(fill="white", colour = "black"))
  }
  
  grid.arrange(grobs=q,nrow=2, ncol=3)
}

Model_Select <- function(fit,Data,method='MISE'){
  err=c();n.data=Data.PL$n.data;n.comb=Data$n.comb;p=ncol(Data$z)
  for (i in 1:length(fit$lambda)){
    temp=matrix(fit$x[,i],Data$n.comb,p)
    err=c(err,mean((temp[(1+n.data):n.comb,]-Data$beta[(1+n.data):n.comb,])^2))
  }
  it=which(err==min(err))
  return(matrix(fit$x[,it],n.comb,p))
}

