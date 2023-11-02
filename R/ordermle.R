#' This provides MLEs for the parameters of a one-way ANCOVA model when the treatment effects are monotonically increasing.
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
#' ordermle(Y,X,k,q)
#' @export
ordermle<-function(data1,data2,k,q)
{
  Y=lapply(data1, function(col)col[!is.na(col)])
  X=lapply(data2, function(col)col[!is.na(col)])
  N=unlist(rbind(lapply(Y,length)))
  yM=unlist(rbind(lapply(Y,mean)))
  xM=unlist(rbind(lapply(X,mean)))
  t2=NULL
  tm1=NULL
  for(h in 1:q)
  {
    for(l in 1:q)
    {
      for(i in 1:k)
      {
        t1=(X[[(h-1)*k+i]]-xM[(h-1)*k+i])*(X[[(l-1)*k+i]]-xM[(l-1)*k+i])
        tm1[i]=sum(t1)
      }
      i2=sum(tm1)
      t2[l+(h-1)*q]=i2
    }
  }
  Sxx=matrix(t2,nrow=q,ncol=q,byrow=TRUE)
  t3=NULL
  tm2=NULL
  for(h in 1:q)
  {
    for(i in 1:k)
    {
      t1=(X[[(h-1)*k+i]]-xM[(h-1)*k+i])*(Y[[i]]-yM[i])
      tm2[i]=sum(t1)
    }
    i2=sum(tm2)
    t3[h]=i2
  }
  Sxy=matrix(t3,nrow=q,ncol=1,byrow=TRUE)
  b0=solve(Sxx)%*%Sxy
  t4=NULL
  t5=NULL
  for(i in 1:k)
  {
    for(h in 1:q)
    {
      t1=b0[h]*xM[(h-1)*k+i]
      t4[h]=t1
    }
    t5[i]=sum(t4)
  }
  a0=yM-t5
  t6=NULL
  t7=NULL
  for(i in 1:k)
  {
    t6=NULL
    for(j in 1:N[i])
    {
      t6[j]=0
    }
    for(h in 1:q)
    {
      t1=b0[h]*(X[[(h-1)*k+i]]-xM[(h-1)*k+i])
      t6=t6+t1
    }
    i2=(Y[[i]]-yM[i])-t6
    t7[i]=sum(i2^2)/(N[i]-q-1)
  }
  B0=b0
  S1=t7
  ###Algorithm (MLEs under full parameters space)
  repeat
  {
    t9=NULL;W=NULL;t8=NULL
    for(i in 1:k)
    {
      for(h in 1:q)
      {
        t1=B0[h]*xM[(h-1)*k+i]
        t8[h]=t1
      }
      i3=yM[i]-sum(t8)
      t9[i]=i3
      W[i]=N[i]/S1[i]
    }
    An=Iso::pava(t9,W)
    t10=NULL
    t11=NULL;tm3=NULL
    for(h in 1:q)
    {
      for(l in 1:q)
      {
        for(i in 1:k)
        {
          t1=(X[[(h-1)*k+i]]*X[[(l-1)*k+i]])
          tm3[i]=sum(t1)/S1[i]
        }
        i4=sum(tm3)
        t11[l+(h-1)*q]=i4
      }
    }
    Pxx=matrix(t11,nrow=q,ncol=q,byrow=TRUE)
    t12=NULL;tm4=NULL
    for(h in 1:q)
    {
      for(i in 1:k)
      {
        t1=X[[(h-1)*k+i]]*(Y[[i]]-An[i])
        tm4[i]=sum(t1)/S1[i]
      }
      i2=sum(tm4)
      t12[h]=i2
    }
    Qxy=matrix(t12,nrow=q,ncol=1,byrow=TRUE)
    bn=solve(Pxx)%*%Qxy
    t13=NULL
    t14=NULL
    for(i in 1:k)
    {
      t13=NULL
      for(j in 1:N[i])
      {
        t13[j]=0
      }
      for(h in 1:q)
      {
        t1=bn[h]*X[[(h-1)*k+i]]
        t13=t13+t1
      }
      i2=(Y[[i]]-An[i])-t13
      t14[i]=sum(i2^2)/N[i]
    }
    dif1=max(abs(a0-An));dif2=max(abs(B0-bn))
    if(dif1<=0.00001&dif2<=0.00001)
    {
      break
    }
    a0=An;B0=bn;S1=t14
  }
  mle=c(An,B0,S1)
  return(mle)
}
