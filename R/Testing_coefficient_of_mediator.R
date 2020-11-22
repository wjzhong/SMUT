# Testing_coeffcient_of_mediator is to test the fixed effects(theta) in the outcome model
# outcome model: outcome ~ intercept + G*gamma + mediator*theta + covariates*tau + error
# the p value is for testing whether theta is zero
# null hypothesis is H_0: theta = 0
# alternative hypothesis is H_a: theta != 0

Testing_coefficient_of_mediator=function(G,mediator,outcome,covariates=NULL,outcome_type="continuous",method="score",
              approxi=TRUE, debug=FALSE){
  G = G + matrix(0, nrow = nrow(G), ncol = ncol(G))
  
  if (!is.null(covariates)){
    covariates=as.matrix(covariates)
  }
  
  
  #initial values

  intercept=runif(1,-1,1)
  mu=runif(1,-1,1)
  sigma2_gamma=runif(1,0,1)
  sigma2_error=runif(1,0,1)
  if (!is.null(covariates)){
    tau=runif(ncol(covariates),-1,1)
  }
  
  max_iteration = 100
  iteration = 0
  ll=rep(NA,max_iteration)
  r=NULL
  d1=d2=1
  delta_ll=1
  D=NULL
  V_inv=NULL
  Q_null=NULL
  gamma_hat=NULL
  
  # EM algorithm
  while( (iteration<max_iteration) & ( max(abs(d1),abs(d2))>1e-6)  ){
    iteration=iteration+1
    if (is.null(covariates)){
      previous_param=list("D"=D,"V_inv"=V_inv,"gamma_hat"=gamma_hat,"mu"=mu,"sigma2_gamma"=sigma2_gamma,
                          "Q_null"=Q_null,"intercept"=intercept,"sigma2_error"=sigma2_error)
    }else{
      previous_param=list("D"=D,"V_inv"=V_inv,"gamma_hat"=gamma_hat,"mu"=mu,"sigma2_gamma"=sigma2_gamma,
                          "Q_null"=Q_null,"intercept"=intercept,"sigma2_error"=sigma2_error,"tau"=tau)
    }
    

    # E step
    D=diag(rep(sigma2_gamma,ncol(G)),nrow = ncol(G), ncol = ncol(G) )
    V=eigenMapMatMult(eigenMapMatMult(G, D), t(G))+sigma2_error*diag(nrow(G))
    V_inv = chol2inv(chol(V))
    if (is.null(covariates)){
      gamma_hat =  D %*% (t(G) %*%  (V_inv %*% (outcome - intercept - G %*% (rep(mu,ncol(G))  )  ))) +  mu
    }else{
      gamma_hat =  D %*% (t(G) %*%  (V_inv %*% (outcome - intercept - G %*% (rep(mu,ncol(G))  )  - covariates %*% tau ))) +  mu
    }
    

    # M step

    mu=sum(t(1/(ncol(G))*t(gamma_hat)  ) )
    
    if (is.null(covariates)){
      intercept=mean(outcome-G%*%gamma_hat)
    }else{
      intercept=mean(outcome -covariates %*% tau - G%*%gamma_hat)
      tau=MASS::ginv( t(covariates) %*% covariates ) %*% t(covariates) %*% ( outcome - intercept - G%*%gamma_hat )
    }
    
    tGV_invG=eigenMapMatMult(eigenMapMatMult(t(G),V_inv),G)

    if (is.null(covariates)){
      
      d1=(sigma2_gamma^2/ (ncol(G)) )*(t(outcome-intercept-G%*%(rep(mu,ncol(G))))%*%V_inv%*%G%*%t(G)%*%V_inv%*%(outcome-intercept-G%*%(rep(mu,ncol(G))  ))-sum(diag( (tGV_invG) )))
      sigma2_gamma=sigma2_gamma+d1
      
      d2=(sigma2_error^2/ (nrow(G)) )*(t(outcome-intercept-G%*%(rep(mu,ncol(G))))%*%V_inv%*%V_inv%*%(outcome-intercept-G%*%(rep(mu,ncol(G))  ))-sum(diag( V_inv )))
      sigma2_error=as.numeric(sigma2_error+d2)
      
    }else{
      
      d1=(sigma2_gamma^2/ (ncol(G)) )*(t(outcome-intercept-G%*%(rep(mu,ncol(G))) - covariates %*% tau )%*%V_inv%*%G%*%t(G)%*%V_inv%*%(outcome-intercept-G%*%(rep(mu,ncol(G))) - covariates %*% tau )-sum(diag( (tGV_invG) )))
      sigma2_gamma=sigma2_gamma+d1
      
      d2=(sigma2_error^2/ (nrow(G)) )*(t(outcome-intercept-G%*%(rep(mu,ncol(G)))  - covariates %*% tau )%*%V_inv%*%V_inv%*%(outcome-intercept-G%*%(rep(mu,ncol(G))) - covariates %*% tau )-sum(diag( V_inv )))
      sigma2_error=as.numeric(sigma2_error+d2)
      
    }
    


    D_inv=diag(rep(sigma2_gamma^(-1),ncol(G)),nrow = ncol(G), ncol = ncol(G) )
    if (is.null(covariates)){
      Q_null=-0.5*nrow(G)*log(2*3.1415926*sigma2_error)-1/(2*sigma2_error)*t(outcome-intercept-G%*%gamma_hat)%*%(outcome-intercept-G%*%gamma_hat)-0.5*sum(log(2*3.1415926*diag(D)))-0.5*t(gamma_hat-mu)%*%D_inv%*%(gamma_hat-mu)
      
    }else{
      Q_null=-0.5*nrow(G)*log(2*3.1415926*sigma2_error)-1/(2*sigma2_error)*t(outcome-intercept-G%*%gamma_hat - covariates %*% tau)%*%(outcome-intercept-G%*%gamma_hat - covariates %*% tau)-0.5*sum(log(2*3.1415926*diag(D)))-0.5*t(gamma_hat-mu)%*%D_inv%*%(gamma_hat-mu)
      
    }
    
    Q_null=as.numeric(Q_null)

    ll[iteration]=Q_null
    if (is.null(covariates)){
      current_param=list("D"=D,"V_inv"=V_inv,"gamma_hat"=gamma_hat,"mu"=mu,"sigma2_gamma"=sigma2_gamma,
                         "Q_null"=Q_null,"intercept"=intercept,"sigma2_error"=sigma2_error)
    }else{
      current_param=list("D"=D,"V_inv"=V_inv,"gamma_hat"=gamma_hat,"mu"=mu,"sigma2_gamma"=sigma2_gamma,
                         "Q_null"=Q_null,"intercept"=intercept,"sigma2_error"=sigma2_error,"tau"=tau)
    }
    


    if (iteration>1){
      delta_ll=ll[iteration]-ll[iteration-1]
      if(is.na(delta_ll)){
        delta_ll=-1
      }
    }
    
    if (debug){
      cat("iteration = ",iteration,"\n")
      cat("Q_null = ",Q_null,"\n")
      cat("intercept = ",intercept,"\n")
      if (!is.null(covariates)){
        cat("tau = ",tau,"\n")
      }
      cat("sigma2_error = ",sigma2_error,"\n")
      cat("sigma2_gamma = ",sigma2_gamma,"\n")
      cat("mu = ",mu,"\n")
      cat("max d =",max(abs(d1),abs(d2)),"\n")
      cat("delta_ll =",delta_ll,"\n")
      cat("\n")
    }
    
    if (delta_ll<0){
      current_param=previous_param
    }

  }

  intercept=current_param$intercept
  if (!is.null(covariates)){
    tau=current_param$tau
  }
  mu=current_param$mu
  V_inv=current_param$V_inv
  gamma_hat=current_param$gamma_hat

  # score test
  if (is.null(covariates)){
    first_derivative = as.numeric(t(outcome - intercept - G %*% gamma_hat) %*% V_inv %*% mediator)
  }else{
    first_derivative = as.numeric(t(outcome - intercept - G %*% gamma_hat - covariates %*% tau) %*% V_inv %*% mediator)
  }
  

  Info11 = t(mediator) %*% V_inv %*% mediator
  if (approxi==T){
    Info=Info11
  }else{
    if (is.null(covariates)){
      Info12 = c( t(rep(1,nrow(G)))%*%V_inv%*%mediator,t(G%*%rep(1,ncol(G)))%*%V_inv%*%mediator,0,0 )
      Info22=matrix(0,4,4)
      Info22[1,1]=t(rep(1,nrow(G)))%*%V_inv%*%rep(1,nrow(G))
      Info22[1,2]=Info22[2,1]=t(G%*%rep(1,ncol(G)))%*%V_inv%*%rep(1,nrow(G))
      Info22[2,2]=t(G%*%rep(1,ncol(G)))%*%V_inv%*%(G%*%rep(1,ncol(G)))
      V_invGtG=eigenMapMatMult(eigenMapMatMult(V_inv,G),t(G) )
      Info22[3,3]=0.5*sum(diag( eigenMapMatMult(V_invGtG,V_invGtG) ))
      Info22[3,4]=Info22[4,3]=0.5*sum(diag( eigenMapMatMult(V_inv,V_invGtG) ))
      Info22[4,4]=0.5*sum(diag( eigenMapMatMult(V_inv,V_inv) ))
    }else{
      Info12 = c( t(rep(1,nrow(G)))%*%V_inv%*%mediator,t(G%*%rep(1,ncol(G)))%*%V_inv%*%mediator,0,0, t(covariates) %*% V_inv %*% mediator )
      Info22=matrix(0,4+ncol(covariates),4+ncol(covariates))
      Info22[1,1]=t(rep(1,nrow(G)))%*%V_inv%*%rep(1,nrow(G))
      Info22[1,2]=Info22[2,1]=t(G%*%rep(1,ncol(G)))%*%V_inv%*%rep(1,nrow(G))
      Info22[2,2]=t(G%*%rep(1,ncol(G)))%*%V_inv%*%(G%*%rep(1,ncol(G)))
      V_invGtG=eigenMapMatMult(eigenMapMatMult(V_inv,G),t(G) )
      Info22[3,3]=0.5*sum(diag( eigenMapMatMult(V_invGtG,V_invGtG) ))
      Info22[3,4]=Info22[4,3]=0.5*sum(diag( eigenMapMatMult(V_inv,V_invGtG) ))
      Info22[4,4]=0.5*sum(diag( eigenMapMatMult(V_inv,V_inv) ))
      Info22[(4+1):(4+ncol(covariates)),(4+1):(4+ncol(covariates))]=t(covariates) %*% V_inv %*% covariates
      Info22[1,(4+1):(4+ncol(covariates))]=Info22[(4+1):(4+ncol(covariates)),1]=t(rep(1,nrow(G))) %*% V_inv %*% covariates
      Info22[2,(4+1):(4+ncol(covariates))]=Info22[(4+1):(4+ncol(covariates)),2]=t(G%*%rep(1,ncol(G))) %*% V_inv %*% covariates
    }
    
    Info=Info11-t(Info12)%*%MASS::ginv(Info22)%*%Info12
  }

  score_statistic = abs( (first_derivative)^2/(Info))
  pvalue=as.numeric(1-pchisq(score_statistic,df=1))
  return(pvalue)
}
