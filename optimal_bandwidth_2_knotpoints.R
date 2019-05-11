trun=function(prediktor,knot,orde)
{
  prediktor[prediktor<knot]=knot
  b=(prediktor-knot)^orde
  return(b)
}
mp<-function(x,eps=1e-006)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1])%*%t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)]%*%diag(1/diago)%*%t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}
Dataset=data.matrix(Dataset)
spline=function(data)
{
  #cukup hanya mengganti datt=5-7
  datt=5
  daty=8
  bb=min(unique(Dataset[,datt]))
  ba=max(unique(Dataset[,datt]))
  p=as.numeric(readline("inputkan orde = "))
  jk=2
  w=combn(bb:ba,jk)
  d=ncol(w)
  w=matrix(c(w),d,jk,TRUE)
  n=nrow(data)
  prediktor=data[,datt]
  dataurut=data[order(prediktor),1:(daty+1)]
  x=dataurut[,datt]
  y=dataurut[,daty:(daty+1)]
  X=cbind(1,dataurut[,(1:4)],x)
  M=ncol(y)
  xpar=cbind(1,data[,1:4])
  lebarx=ncol(xpar)
  GCV=rep(0,d)
  MSE=rep(0,d)
  R2=rep(0,d)
  iden=diag(1,n)
  iden2=diag(1,M*n)
  OLS=matrix(0,n,2)
  ytopi=matrix(0,n,2)
  res=matrix(0,n,2)
  OLS=mp(t(X)%*%X)%*%t(X)%*%y
  ytopi=X%*%OLS
  res=y-ytopi
  sigmaht=matrix(0,2,2)
  W=matrix(0,2*n,2*n)
  for (i in 1:n)
  {
    sigmaht[1,1]=sigmaht[1,1]+res[i,1]*res[i,1]
    sigmaht[1,2]=sigmaht[1,2]+res[i,1]*res[i,2]
    sigmaht[2,1]=sigmaht[2,1]+res[i,2]*res[i,1]
    sigmaht[2,2]=sigmaht[2,2]+res[i,2]*res[i,2]
  }
  W[1:n,1:n]=sigmaht[1,1]*iden
  W[(n+1):(M*n),1:n]=sigmaht[1,2]*iden
  W[(n+1):(n*M),(n+1):(n*M)]=sigmaht[2,2]*iden
  W[1:n,(n+1):(M*n)]=W[(n+1):(M*n),1:n]
  W=mp(W)
  ye=matrix(0,M*n)
  for (i in 1:n)
  {
    ye[i]=y[i,1]
    ye[i+n]=y[i,2]
  }
  cat("======================================================================= \n")
  cat("      TITIK \n")
  cat(" KNOT[1]   KNOT[2]          MSE               GCV            R2 \n")
  cat("======================================================================= \n")
  for(j in 1:d)
  {
    v1=matrix(0,n,p+1)
    for(i in 1:p)
    {
      v1[,i]=x^(i-1)
      v1[,p+1]=x^(p)
    }
    v2=matrix(0,n,jk)
    for(i in 1:jk)
    {
      v2[,i]=trun(x,w[j,i],p)
    }
    Te=cbind(v1[,-1],v2)
    Lte=ncol(Te)
    Xe=cbind(xpar,Te)
    X=matrix(0,M*n,M*(Lte+lebarx))
    for (i in 1:n)
    {
      for (k in 1:(Lte+lebarx))
      {
        X[i,k]=Xe[i,k]
        X[i+n,k+Lte+lebarx]=Xe[i,k]
      }
    }
    betatopi=mp(t(X)%*%W%*%X)%*%t(X)%*%W%*%ye
    ytopi=X%*%betatopi
    H=X%*%mp(t(X)%*%W%*%X)%*%t(X)%*%W
    MSE[j]=(t(ye-ytopi)%*%(ye-ytopi))/(n*M)
    GCV[j]=MSE[j]/(1-((1/M*n)*sum(diag(H))))^2
    JKT=t(ye-(mean(ye)))%*%(ye-(mean(ye)))
    JKG=t(ye-ytopi)%*%(ye-ytopi)
    R2[j]=1-(JKG/JKT)
    cat("  ",format(w[j,1]),"      ",format(w[j,2]),"\t\t",format(MSE[j]),"\t ",format(GCV[j]),"\t ",format(R2[j]),"\n")
  }
  s<-matrix(c(w[,1],w[,2],MSE,GCV,R2),d,5)
  GCVmin<-min(GCV)
  knot<-rep(0,jk)
  knot[1]<-s[s[,4]==GCVmin,1]
  knot[2]<-s[s[,4]==GCVmin,2]
  mse<-s[s[,4]==GCVmin,3]
  rr<-s[s[,4]==GCVmin,5]
  cat("======================================================================= \n")
  cat(" Berdasarkan data diatas dapat diperoleh : \n")
  cat(" GCV Minimum = ",GCVmin,"\n")
  cat(" Pada : titik knot [1] =",knot[1],"\n")
  cat("        titik knot [2] =",knot[2],"\n \n")
  knot=seq(1,d,1)
  plot(knot,GCV,type="l",lwd=2)
}
spline(Dataset)
