library(MASS)
library(ghyp)
start.time <- Sys.time()
set.seed(1)
  n=1000;#Size
  d=10; #Dim
  g=100
  tau = 0.5
  beta0=matrix(0,d,1);
  B0=diag(d);
  B0=diag(d);
  #初始x分布
  x=t(mvrnorm(n, beta0, B0));
  m=10
  sample = n/m #Sample Size
  beta_true=c(1,3,2,10,4,3,-1,4,5,0) #beta
  y=matrix(0,n,1);
  Y=matrix(0,n,1);
  delta0=1;
  alpha0=1;
  M=1000;
  beta = matrix(0,d,M);
  sigma = vector(length = M);
  s=matrix(0,d,M)
  s1=matrix(0,d,M)
  v=matrix(0,n,M); 
  iota=matrix(0,d,d);
  zeta=matrix(0,d,1);
  sigma_true=1;
  y=matrix(0,n,1);
  Y=matrix(0,n,1);
  #tau的赋值
  Beta1=matrix(0,length(tau),d);
  Var1=matrix(0,length(tau),d);
  sigma1=matrix(0,length(tau),1);

  x=t(mvrnorm(n, beta0, B0));

  for (i in 1:n){
    tmp0=rexp(1,(tau*(1-tau))*1/sigma_true)
    y[i]=t(x[ ,i])%*%beta_true+(1-2*tau)*tmp0+rnorm(1,0,1)*sqrt(2*tmp0*sigma_true)
  }
  sigma[1] = 1/rgamma(1,shape=delta0,scale=1/alpha0);
  v[,1]= rexp(n,(tau*(1-tau))*1/sigma[1]);
  beta[,1]=mvrnorm(1, beta0, B0);
  s[,1]=rgamma(1,shape=delta0,scale=1/alpha0);
  s1[,1]=rgamma(1,shape=delta0,scale=1/alpha0);
  Y<-y
  X<-x
  k=1
  for (j in 1:(M-1)){
    #v update
    psi=1/(2*sigma[j]);
    chi=0;
    for(i in ((k-1)*sample+1):(sample*k)){
      chi=((y[i]-t(x[ ,i])%*%beta[ ,j])^2+beta[ ,j]%*%x[ ,i]%*%t(x[ ,i])%*%beta[ ,j]/g)/(2*sigma[j])
      v[i,j+1]=rgig(n=1, lambda = 0, psi=psi, chi=chi)
    } 
    #x,y update
    for (i in ((k-1)*sample+1):(sample*k)){
      Y[i]=1/sqrt(2)*(y[i]-(1-2*tau)*v[i,j+1])
    }
    X[,((k-1)*sample+1):(k*sample)]=1/sqrt(2)*x[,((k-1)*sample+1):(k*sample)]
    tmp1=matrix(0,d,d);
    tmp2=matrix(0,d,1);
    for (i in ((k-1)*sample+1):(sample*k)){
      tmp1=tmp1+X[ ,i]%*%t(X[ ,i])/v[i,j+1];
      tmp2=tmp2+X[ ,i]*Y[i]/v[i,j+1];
    }
    #NIG
    mu = solve((1+1/g)*(tmp1))%*%(tmp2)
    Lambda = (1+1/g)*((tmp1))
    a=3*sample/2
    Z<-diag((v[((k-1)*sample+1):(sample*k),j+1])^(-1))
    b=1/2*Y[((k-1)*sample+1):(sample*k)]%*%Z%*%Y[((k-1)*sample+1):(sample*k)]-1/2*t(mu)%*%Lambda%*%mu+tau*(1-tau)*sum(v[((k-1)*sample+1):(sample*k),j+1])
    ksi=d/2+a;
    phi =b+1/2*t(beta[ ,j]-mu)%*%Lambda%*%(beta[ ,j]-mu)
    sigma[j+1]=1/rgamma(1,shape=ksi,scale = 1/phi);
    #beta update
    iota=sigma[j+1]*solve(Lambda)
    zeta=mu
    beta[ ,j+1] = mvrnorm(1,zeta,iota);
  }
  
  for (i in 1:d){
    Beta1[1,i]=mean(beta[i,c((0.9*M):M)]);
  }
  sigma1=mean(sigma[c((0.9*M):M)]);
  
  sigma[1]=sigma1
  beta[ ,1]=Beta1
  ## Gibbs sampler
  for(k in 2:m){
    mu0=mu
    Lambda0=Lambda
    a0=a
    b0=b
    for (j in 1:(M-1)){
      #v update
      psi=1/(2*sigma[j]);
      chi=0;
      for(i in ((k-1)*sample+1):(sample*k)){
        chi=((y[i]-t(x[ ,i])%*%beta[ ,j])^2+beta[ ,j]%*%x[ ,i]%*%t(x[ ,i])%*%beta[ ,j]/g)/(2*sigma[j])
        v[i,j+1]=rgig(n=1, lambda = 0, psi=psi, chi=chi)
      } 
      #x,y的update
      for (i in ((k-1)*sample+1):(sample*k)){
        Y[i]=1/sqrt(2)*(y[i]-(1-2*tau)*v[i,j+1])
      }
      X[,((k-1)*sample+1):(k*sample)]=1/sqrt(2)*x[,((k-1)*sample+1):(k*sample)]
      Z<-diag((v[((k-1)*sample+1):(sample*k),j+1])^(-1))
      tmp1=matrix(0,d,d);
      tmp2=matrix(0,d,1);
      for (i in ((k-1)*sample+1):(sample*k)){
        tmp1=tmp1+X[ ,i]%*%t(X[ ,i])/v[i,j+1];
        tmp2=tmp2+X[ ,i]*Y[i]/v[i,j+1];
      }
      #NIG
      Lambda= Lambda0+tmp1+1/g*(tmp1)
      mu = ginv(Lambda)%*%(Lambda0%*%mu0+(tmp1)%*%ginv(tmp1)%*%(tmp2))
      a= a0+(sample-k-2)/2+(k+2)/2+sample
      tmp3=0;
      for(i in ((k-1)*sample+1):(sample*k)){
        tmp3=tmp3+(Y[i]-t(X[ ,i])%*%(ginv(tmp1)%*%(tmp2)))^2/v[i,j+1];
      }
      b=b0+1/2*tmp3+1/2*t(mu0-mu)%*%Lambda0%*%(mu0-mu)+1/2*t(mu)%*%(tmp1/g)%*%mu+1/2*t(ginv(tmp1)%*%(tmp2)-mu)%*%(tmp1)%*%(ginv(tmp1)%*%(tmp2)-mu)+tau*(1-tau)*sum(v[((k-1)*sample+1):(sample*k),j+1])
      #sigma update
      ksi=d/2+a;
      phi =b+1/2*t(beta[ ,j]-mu)%*%Lambda%*%(beta[ ,j]-mu)
      sigma[j+1]=1/rgamma(1,shape=ksi,scale = 1/phi);
      #beta update
      iota=sigma[j+1]*solve(Lambda)
      zeta=mu
      beta[ ,j+1] = mvrnorm(1,zeta,iota);
    }
  }
  for (i in 1:d){
    Beta1[1,i]=mean(beta[i,c((0.9*M):M)]);
  }
  end.time <- Sys.time()
  time_taken <- end.time - start.time
  print(paste("Running time：", time_taken))