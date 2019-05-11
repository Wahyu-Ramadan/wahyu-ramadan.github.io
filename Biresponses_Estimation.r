#fungsi truncated
trun=function(prediktor,knot,orde)
{
  prediktor[prediktor<knot]=knot
  b=(prediktor-knot)^orde
  return(b)
}

#fungsi invers
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

#syntax agar data terbaca sebagai matrix
Dataset=data.matrix(Dataset)

#fungsi estimasi
spline=function(data)
{
  #kerakteristik data yg digunakan: x=1:5, t=6:9, y=10,11 (datx=data parametrik, datt=data nonparametrik, daty=respon)
  n=nrow(data)
  datx=1
  datt=5
  daty=8
  iden=diag(1,n)
  x=data[,datx:(daty-1)]
  xp=x[,-(datt:(daty-1))]
  xp=cbind(1,xp)
  
  #te merupakan komponen matriks prediktor nonparametrik
  te=data[,(datt:(daty-1))]
  n_nonp=ncol(te)
  y=data[,daty:(daty+1)]
  X=cbind(1,x)
  M=ncol(y)
  pe=datt-datx
  OLS=matrix(0,n,2)
  ytopi=matrix(0,n,2)
  res=matrix(0,n,2)
  OLS=mp(t(X)%*%X)%*%t(X)%*%y
  ytopi=X%*%OLS
  res=y-ytopi
  sigmaht=matrix(0,2,2)
  
  #matriks y birespon
  
  ye=matrix(0,M*n)
  for (i in 1:n)
  {
    ye[i]=y[i,1]
    ye[i+n]=y[i,2]
  }
  
  #matriks pembobot w
  
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
  
  et=ncol(te)
  p=rep(0,n_nonp)
  k=rep(0,n_nonp)
  w=matrix(0,3,3)
  w=matrix(c(2,3,4,5,6,17,85,86,87),3,3,TRUE) #masukan matriks dengan hasil titik knot 
  #p adalah orde, k adalah jumlah knot
  
  p=c(2,2,1)
  k=c(3,3,3)
  
  #matriks variabel nonparametrik
  #untuk t1
  
  v1<-matrix(0,n,p[1]+1)
  for(i in 1:p[1])
  {
    v1[,i]=te[,1]^(i-1)
    v1[,p[1]+1]=te[,1]^(p[1])
  }
  v2<-matrix(0,n,k[1])
  for(j in 1:k[1])
  {
    v2[,j]<-trun(te[,1],w[1,j],p[1])
  }
  T1=cbind(v1[,-1],v2)
  
  #untuk t2
  
  v1<-matrix(0,n,p[2]+1)
  for(i in 1:p[2])
  { 
    v1[,i]<-te[,2]^(i-1)
    v1[,p[2]+1]=te[,2]^(p[2])
  }
  v2<-matrix(0,n,k[2])
  for(j in 1:k[2])
  {
    v2[,j]<-trun(te[,2],w[2,j],p[2])
  }
  T2=cbind(v1[,-1],v2)
  
  #untuk t3
  
  v1<-matrix(0,n,p[3]+1)
  for(i in 1:p[3])
  {
    v1[,i]<-te[,3]^(i-1)
    v1[,p[3]+1]=te[,3]^(p[3])
  }
  v2<-matrix(0,n,k[3])
  for(j in 1:k[3])
  {
    v2[,j]<-trun(te[,3],w[3,j],p[3])
  }
  T3=cbind(v1[,-1],v2)
  
  Tsem=cbind(T1,T2,T3)
  c1=cbind(xp,Tsem)
  pp=ncol(xp)
  pe=sum(p)
  ka=sum(k)
  
  #matriks C (Variabel gabunga/n variabel nonpar & para)
  print(pe+ka+pp)
  print(n)
  C=matrix(0,M*n,((pe+ka+pp)*M))
  for (i in 1:n)
  {
    for (j in 1:(pe+ka+pp))
    {
      C[i,j]=c1[i,j]
      C[i+n,j+pe+ka+pp]=c1[i,j]
    }
  }
  
  #Estimasi
  
  delta<-mp(t(C)%*%W%*%C)%*%t(C)%*%W%*%ye
  ytopi<-C%*%delta
  error=ye-ytopi
  Asemi<-C%*%mp(t(C)%*%W%*%C)%*%t(C)%*%W
  MSE<-(t(ye-ytopi)%*%(ye-ytopi))/(M*n)
  JKT<-t(ye-(mean(ye)))%*%(ye-(mean(ye)))
  JKG<-t(ye-ytopi)%*%(ye-ytopi)
  R2<-1-(JKG/JKT)
  delta=round(delta,2)
  delta2=matrix(c(delta),19,2,FALSE)
  colnames(delta2)<-c("Y1","Y2")
  rownames(delta2)<-c("delta1","delta2","delta3","delta4","delta5","delta6","delta7","delta8","delta9","delta10","delta11","delta12","delta13","delta14","delta15","delta16","delta17","delta18","delta19")
  cat("\n\n    NILAI PARAMETER ")
  cat("\n====================================================\n")
  print(delta2)
  cat("\n\n    PENDEKATAN SEMIPARAMETRIK ")
  cat("\n====================================================\n")
  cat("    MSE\t\t        R2\t\t      \n")
  cat("====================================================\n")
  cat(MSE,"\t     ",R2,"\t     \n\n\n\n")
  respon_1=ye[1:140]
  respon_2=ye[141:280]
  ytopi1=ytopi[1:140]
  ytopi2=ytopi[141:280]
  win.graph()
  Unit_Observasi=seq(1,length(respon_1),1)
  plot(Unit_Observasi,respon_1,main = "Plot estimasi nilai UNBK Bahasa Inggris",col="red",type="p")
  lines(Unit_Observasi,ytopi1,main = "Plot estimasi nilai UNBK Bahasa Inggris",col="blue",type="l",lwd=3)
  win.graph()
  plot(Unit_Observasi,respon_2,main = "Plot estimasi nilai UNBK KK",col="red",type="p")
  lines(Unit_Observasi,ytopi2,main = "Plot estimasi nilai UNBK KK",col="blue",type="l",lwd=3)
  
  #pendekatan paramterik
  d=cbind(xp,te)
  pe=ncol(d)
  D=matrix(0,M*n,((pe)*M))
  for (i in 1:n)
  {
    for (j in 1:(pe))
    {
      D[i,j]=d[i,j]
      D[i+n,j+pe]=d[i,j]
    }
  }
  delta<-mp(t(D)%*%D)%*%t(D)%*%ye
  ytopi2<-D%*%delta
  error=ye-ytopi
  Asemi<-D%*%mp(t(D)%*%D)%*%t(D)
  MSE<-(t(ye-ytopi2)%*%(ye-ytopi2))/(M*n)
  JKT<-t(ye-(mean(ye)))%*%(ye-(mean(ye)))
  JKG<-t(ye-ytopi2)%*%(ye-ytopi2)
  R2<-1-(JKG/JKT)
  cat("    PENDEKATAN PARAMETRIK ")
  cat("\n====================================================\n")
  cat("    MSE\t\t        R2\t\t      \n")
  cat("====================================================\n")
  cat(MSE,"\t     ",R2,"\t     \n\n")
  
  #regresi parsial
  delta<-mp(t(C)%*%W%*%C)%*%t(C)%*%W%*%ye
  E=C
  cc=ncol(E)
  for(i in 1:(cc/2))
  {
    E[1:140,i]=C[1,i]
    E[141:280,(i+19)]=C[1,i]
  }
  
  #parsial rapot
  G=E
  G[,5]=C[,5]
  G[,24]=C[,24]
  ytopi<-G%*%delta
  respon_1=ye[1:140]
  respon_2=ye[141:280]
  ytopi1=ytopi[1:140]
  ytopi2=ytopi[141:280]
  win.graph()
  Unit_Observasi=seq(1,length(respon_1),1)
  plot(Unit_Observasi,respon_1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh nilai rapot",col="red",type="p")
  lines(Unit_Observasi,ytopi1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh nilai rapot",col="blue",type="l",lwd=3)
  win.graph()
  plot(Unit_Observasi,respon_2,main = "Plot estimasi parsial nilai UNBK KK oleh nilai rapot",col="red",type="p")
  lines(Unit_Observasi,ytopi2,main = "Plot estimasi parsial nilai UNBK KK oleh nilai rapot",col="blue",type="l",lwd=3)
  
  #parsial jarak
  G=E
  G[,6:10]=C[,6:10]
  G[,25:29]=C[,25:29]
  ytopi<-G%*%delta
  respon_1=ye[1:140]
  respon_2=ye[141:280]
  ytopi1=ytopi[1:140]
  ytopi2=ytopi[141:280]
  win.graph()
  Unit_Observasi=seq(1,length(respon_1),1)
  plot(Unit_Observasi,respon_1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh jarak sekolah",col="red",type="p")
  lines(Unit_Observasi,ytopi1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh jarak sekolah",col="blue",type="l",lwd=3)
  win.graph()
  plot(Unit_Observasi,respon_2,main = "Plot estimasi parsial nilai UNBK KK oleh jarak sekolah",col="red",type="p")
  lines(Unit_Observasi,ytopi2,main = "Plot estimasi parsial nilai UNBK KK oleh jarak sekolah",col="blue",type="l",lwd=3)
  
  #parsial pend.ortu
  G=E
  G[,11:15]=C[,11:15]
  G[,30:34]=C[,30:34]
  ytopi<-G%*%delta
  respon_1=ye[1:140]
  respon_2=ye[141:280]
  ytopi1=ytopi[1:140]
  ytopi2=ytopi[141:280]
  win.graph()
  Unit_Observasi=seq(1,length(respon_1),1)
  plot(Unit_Observasi,respon_1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh lamanya pendidikan orang tua",col="red",type="p")
  lines(Unit_Observasi,ytopi1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh lamanya pendidikan orang tua",col="blue",type="l",lwd=3)
  win.graph()
  plot(Unit_Observasi,respon_2,main = "Plot estimasi parsial nilai UNBK KK oleh lamanya pendidikan orang tua",col="red",type="p")
  lines(Unit_Observasi,ytopi2,main = "Plot estimasi parsial nilai UNBK KK oleh lamanya pendidikan orang tua",col="blue",type="l",lwd=3)
  
  #parsial US
  G=E
  G[,16:19]=C[,16:19]
  G[,35:38]=C[,35:38]
  ytopi<-G%*%delta
  respon_1=ye[1:140]
  respon_2=ye[141:280]
  ytopi1=ytopi[1:140]
  ytopi2=ytopi[141:280]
  win.graph()
  Unit_Observasi=seq(1,length(respon_1),1)
  plot(Unit_Observasi,respon_1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh nilai ujian sekolah",col="red",type="p")
  lines(Unit_Observasi,ytopi1,main = "Plot estimasi parsial nilai UNBK Bahasa Inggris oleh nilai ujian sekolah",col="blue",type="l",lwd=3)
  win.graph()
  plot(Unit_Observasi,respon_2,main = "Plot estimasi parsial nilai UNBK KK oleh nilai ujian sekolah",col="red",type="p")
  lines(Unit_Observasi,ytopi2,main = "Plot estimasi parsial nilai UNBK KK oleh nilai ujian sekolah",col="blue",type="l",lwd=3)
  
}
