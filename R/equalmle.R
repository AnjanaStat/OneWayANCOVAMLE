#' This provides MLEs for the parameters of a one-way ANCOVA model when the treatment effects are identical.
#'
#' More detailed description
#'
#' @param data1 a real data set
#' @param data2 a real data set
#' @param k a positive integer
#' @param q a positive integer
#'
#' @return numeric vector
#'
#' @examples
#' k=4;q=3
#' N=c(20,10,10,20,20,10,10,20,20,10,10,20);S=c(1,1,2,1,1,2,3,1,2,4,6,3,2,4)
#' g=NULL
#' for(i in 1:(k*q))
#' {
#'  g[[i]]=rnorm(N[i],0,sqrt(S[i]))
#' }
#' X=g
#' G2=NULL
#' N=c(20,10,10,20);a=c(1,2,3,4)
#' for(i in 1:k)
#' {
#'  G2[[i]]=rnorm(N[i],a[i],sqrt(S[i]))
#' }
#' Y=G2
#' equalmle(Y,X,k,q)
#' @export
equalmle<-function(data1,data2,k,q)
{
  Y=lapply(data1, function(col)col[!is.na(col)])
  X=lapply(data2, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(Y,length)))
  yM=unlist(rbind(lapply(Y,mean)))
  xM=unlist(rbind(lapply(X,mean)))
  t15=NULL;tm5=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=N[i]*xM[(h-1)*k+i]
      tm5[i]=t1
    }
    t15[h]=sum(tm5)/sum(N)
  }
  t16=N*yM
  totY=sum(t16)/sum(N)
  t17=NULL
  tm6=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(i in 1:k)
      {
        t1=(X[[(h-1)*k+i]]-t15[h])*(X[[(l-1)*k+i]]-t15[l])
        tm6[i]=sum(t1)
      }
      i2=sum(tm6)
      t17[l+(h-1)*q]=i2
    }
  }
  sxx=matrix(t17,nrow=q,ncol=q,byrow=TRUE)
  t18=NULL;tm7=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=(X[[(h-1)*k+i]]-t15[h])*(Y[[i]]-totY)
      tm7[i]=sum(t1)
    }
    i2=sum(tm7)
    t18[h]=i2
  }
  sxy=matrix(t18,nrow=q,ncol=1,byrow=TRUE)
  b01=solve(sxx)%*%sxy
  a01=totY-sum(b01*t15)
  t19=NULL
  t20=NULL
  for(i in 1:k)
  {
    t19=NULL
    for(j in 1:N[i])
    {
      t19[j]=0
    }
    for(h in 1:q)
    {
      t1=b01[h]*(X[[(h-1)*k+i]]-xM[(h-1)*k+i])
      t19=t19+t1
    }
    i2=(Y[[i]]-yM[i])-t19
    t20[i]=sum(i2^2)/(N[i]-q-1)
  }
  B01=b01
  S2=t20
  ###Algorithm (MLEs under null parameter space)
  repeat
  {
    U=NULL
    for(i in 1:k)
    {
      U[i]=N[i]/S2[i]
    }
    tm8=NULL;tm9=NULL
    for(i in 1:k)
    {
      for(h in 1:q)
      {
        t1=B01[h]*xM[(h-1)*k+i]
        tm8[h]=t1
      }
      i3=yM[i]-sum(tm8)
      tm9[i]=i3
    }
    an1=sum(U*tm9)/sum(U)
    t21=NULL
    t22=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(i in 1:k)
        {
          t1=X[[(h-1)*k+i]]*X[[(l-1)*k+i]]
          t21[i]=sum(t1)/S2[i]
        }
        i2=sum(t21)
        t22[l+(h-1)*q]=i2
      }
    }
    pxx=matrix(t22,nrow=q,ncol=q,byrow=TRUE)
    t23=NULL;t24=NULL
    for(h in 1:q)
    {
      for(i in 1:k)
      {
        t1=(Y[[i]]-an1)*X[[(h-1)*k+i]]
        t23[i]=sum(t1)/S2[i]
      }
      i2=sum(t23)
      t24[h]=i2
    }
    qxy=matrix(t24,nrow=q,ncol=1,byrow=TRUE)
    bn1=solve(pxx)%*%qxy
    t25=NULL
    t26=NULL
    for(i in 1:k)
    {
      t25=NULL
      for(j in 1:N[i])
      {
        t25[j]=0
      }
      for(h in 1:q)
      {
        t1=bn1[h]*X[[(h-1)*k+i]]
        t25=t25+t1
      }
      i2=(Y[[i]]-an1)-t25
      t26[i]=sum(i2^2)/N[i]
    }
    dif3=abs(a01-an1);dif4=max(abs(B01-bn1))
    if(dif3<=0.00001&dif4<=0.00001)
    {
      break
    }
    a01=an1;B01=bn1;S2=t26
  }
  value=c(a01,bn1,S2)
  return(value)
}
