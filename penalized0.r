 source(paste("penalized_0.9-45/penalized/R/checkinput.R", sep=''));
 source(paste("penalized_0.9-45/penalized/R/logit.R", sep=''));
 source(paste("penalized_0.9-45/penalized/R/core.R", sep=''));
# source(paste("penalized_0.9-45/penalized/R/penfit.R", sep=''));


####################################
# Fits the penalized regression model
####################################                                    
penalized0 <- function(response, penalized, GRPlen, unpenalized, lambda1=0, lambda2=0, positive = FALSE, data, fusedl=FALSE,  
  model = c("cox", "logistic", "linear", "poisson"), startbeta, startgamma, steps =1, epsilon = 1e-10, 
  maxiter, standardize = FALSE, trace = TRUE) {

  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive) 25 else Inf

  # Park and Hastie type steps?
  if (steps == "Park" || steps == "park") {
    steps <- 1
    park <- TRUE
  } else park <- FALSE

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame()) 
  
 if(standardize)  prep$X<-scale(prep$X) 
  
  
  # check for the presence of penalty parameters
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) && all(lambda2 == 0) && !any(prep$positive))
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)$fit
  
   # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl

  # make weights for lambda1 that is a vector
  if (length(lambda1) == pp && (!all(lambda1==0))) {
    wl1 <- c(numeric(pu), lambda1)
    lambda1 <- 1
  } else {
    wl1 <- 1
  }
  if (length(lambda2) == pp)
    lambda2 <- c(numeric(pu), lambda2)

  # If a steps argument is given, determine where to start
  if (park || steps > 1 && fusedl == FALSE) {
    if (pu > 0)
      lp <- drop(prep$X[,1:pu,drop=FALSE] %*% prep$nullgamma)
    else 
      lp <- numeric(n)
    chck <- (wl1 > 0) & c(rep(FALSE, pu), rep(TRUE, pp))
    gradient <- drop(crossprod(prep$X[,chck,drop=FALSE], fit(lp)$residuals))
    if (length(wl1)>1) {
      rel <- gradient / (wl1[chck] * prep$baselambda1[chck])
    } else {
      rel <- gradient / (wl1 * prep$baselambda1[chck])
    }
    from <- max(ifelse(prep$positive[chck],  rel, abs(rel)))
    if (from < lambda1) {
      warning("Chosen lambda1 greater than maximal lambda1: \"steps\" argument ignored")
      steps <- 1
      park <- FALSE
      from <- lambda1
    }
  } else {
    from <- lambda1
  }
  lambda1s <- seq(from, lambda1, length.out=steps)

  # fit the model for all lambdas
  beta <- prep$beta
  louts <- if (park) 4*pp else length(lambda1s)
  outs <- vector("list", louts)
  rellambda1 <- lambda1s[1]
  ready <- FALSE
  i <- 0
  while (!ready) {
    ready <- (rellambda1 == lambda1)
    i <- i+1 
   
   if(fusedl){
       out<-fusedlasso0(beta = beta, GRPlen, lambda1 = rellambda1 * wl1 * prep$baselambda1,
        lambda2 = lambda2 * prep$baselambda2, chr = prep$chr, positive = prep$positive, X = prep$X,
        fit = fit, trace = trace, epsilon = epsilon, maxiter = maxiter)
      } 
    if (trace) cat("\n")
    beta <- out$beta
    
    if (!ready) {
      if(!fusedl){
      if (park) {
        newpark <- .park(beta = beta, lambda = rellambda1 * wl1 * prep$baselambda1,
            lambda2 = 0, positive = prep$positive, X = prep$X, fit = out$fit)
        rellambda1 <- rellambda1 * (1-newpark$hh)
        if (rellambda1 < lambda1 || rellambda1 == Inf) {
          rellambda1 <- lambda1
          beta <- out$beta
        } else {
          beta <- newpark$beta
        }
        lambda1s <- c(lambda1s, rellambda1)
      } else {
        rellambda1 <- lambda1s[i+1]
        beta <- out$beta
      }
    } else {rellambda1 <- lambda1s[i+1]
            beta <- out$beta}
   }
   
    
    
    outs[[i]] <- out
  }
  
  
  if(length(outs)==1) 
    outs <- outs[[1]]

  outs
}

###############################################################################################
### Within group difference function
###############################################################################################
Diffw<-function(x, GRPlen)
{
    pos_end<-cumsum(GRPlen)  #last dummy variable positions;
	pos_begin<-c(1, pos_end[-length(pos_end)] +1)
 diff1<-NULL  
  for(i in 1:length(GRPlen))
  { 
   xx <- x[pos_begin[i]:pos_end[i]]
  diff1<-c(diff1, diff(xx),0) 
  }
  diff1<-diff1[-length(diff1)]
  return(diff1)
}
   
   
 

###############################################################################################

###Fused Lasso
####################

###############################################################################################
fusedlasso0 <- function(beta, GRPlen, chr, lambda1, lambda2,fit, X, positive, trace = FALSE,
  epsilon=1e-5, maxiter=Inf) {

  # Input:
  #  beta: a vector of length m (say) : starting values
  #  lambda1 and lambda2: vectors of length m
  #  X : data matrix
  #  positive: logical vector denoting the coefficients which should pe penalized
  #   Should return a list with at least:
  #     W:          The weights matrix, or its diagonal
  #     loglik:     The unpenalized loglikelihood, numeric)
  #     residuals:  The residuals
  # If the model includes an Intercept then starting with all the beta coefficients
  # equal to zero might slow down the convergence. Hence it's better to use a
  # warm start

  m=length(beta)
  n=nrow(X) 

   pos_end<-cumsum(GRPlen)  #last dummy variable positions;
 
# find regression coefficients free of L1-penalty or positivity restraint
  
  free <- (lambda1 == 0 & lambda2==0) & !positive

# initialize
  
  LL <- -Inf
  penalty <- penalty1 <- penalty2 <- Inf
  converged = FALSE
  active <- !logical(m)

  if(any(free)){direc=numeric(length(free[-which(free)]))}
  if(!any(free)){direc=numeric(length(beta))}
  if(all(free)){direc=numeric(length(free))}

  actd = (diff(direc))!=0   #actually not used for anything;
  checkd = matrix(0,1,2)
  
  nvar <- m
  oldmo=length(beta)
  
  tryNR <- FALSE
  NRfailed <- FALSE
  finish <- FALSE
  newfit <- TRUE
  
  retain <-0
  cumsteps <- 0
  iter <- 0
 #iterate

 while(!finish){
  #while(iter<110&!finish){
            nzb = (beta!=0)
            if(any(free)){direc=numeric(length(free[-which(free)]))}
            if(!any(free)){direc=numeric(length(beta))}
            if(all(free)){direc=numeric(length(free))}
    # calculate the local likelihood fit
  
if (newfit) {
        activeX <- X[,nzb, drop=FALSE]
        linpred <- drop(activeX %*% beta[nzb])
        lp=linpred
        localfit=fit(lp)
 # Check for divergence
      if (is.na(localfit$loglik)) {
        if (trace) {
          cat(rep("\b", trunc(log10(nvar))+1), sep ="")
          warning("Model does not converge: please increase lambda.", call.=FALSE)
        }
        converged <- FALSE
        break
      }
    grad <- drop(crossprod(X,localfit$residuals))
    oldLL <- LL
    oldpenalty <- penalty
    LL <- localfit$loglik
    penalty1=sum(lambda1[active]*abs(beta[active]))
    penalty2 = 0
  
if((!any(chr==0) && length(table(chr))<=1) || (any(chr ==0) && length(table(chr))<=2)){
   #   penalty2= sum(lambda2[1:(length(lambda2)-1)]*abs(((diff(beta)))))  # change to difference within groups;
	  penalty2= sum(lambda2[1:(length(lambda2)-1)]*abs(((Diffw(beta, GRPlen))))) 
      penalty=penalty1+penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }

if((!any(chr==0) && length(table(chr))>1) || (any(chr==0) && length(table(chr))>2)){
     chrn=chr[!free]
      
     for (chri in 1:length((as.numeric(names(table(chrn)))))){
      betac=beta[!free]
      lambda2c=lambda2[!free]
      beta1.1 = betac[which(chrn==(as.numeric(names(table(chrn))))[chri])]
      lambda2.1 = lambda2c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
      penalty21= sum(lambda2.1[1:(length(lambda2.1)-1)]*abs(((diff(beta1.1))))) #need to change;;
      penalty2= penalty2 + penalty21
      }
  
      penalty <- penalty1 + penalty2
      finishedLL <- (2 * abs(LL - oldLL) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      finishedpen <- (2 * abs(penalty - oldpenalty) / (2 * abs(LL - penalty) + 0.1) < epsilon)
      cumsteps <- 0
    }
}
 
    
# Calculate the penalized gradient from the likelihood gradient
    ###Initialize###
   
     

    #####  Core ########################
chrn=chr[!free]

if(!(all(free))){

       gradc=grad[!free]
       nzbc=nzb[!free]
       betac=beta[!free]
       lambda1c=lambda1[!free]
       lambda2c=lambda2[!free]
       positivec=positive[!free]

   for (chri in 1:length((as.numeric(names(table(chrn)))))){

     beta.1 = betac[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     grad.1 = gradc[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     lambda1.1 = lambda1c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     lambda2.1 = lambda2c[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     nzb.1 = nzbc[which(chrn==(as.numeric(names(table(chrn))))[chri])]
     positive.1 = positivec[which(chrn==(as.numeric(names(table(chrn))))[chri])]

     gdir <- getdirec0(grad=grad.1,nzb=nzb.1,beta=beta.1, GRPlen, lambda1=lambda1.1,lambda2=lambda2.1,positive=positive.1)
     direc[which(chrn==(as.numeric(names(table(chrn))))[chri])] = gdir$direc_i
        } 
     if(any(free)){direc=c(grad[which(free)],direc)}
	 
}  else {direc=grad}
 
     oldcheckd = checkd
     oldactive=active
     active= direc!=0 | nzb

   
     checkd = nonZero0(direc[!free])
     actd = (diff(direc[!free]))!= 0 #actually not used for anything;


     activdir=direc[active]
     activbeta=beta[active]
 
  # check if retaining the old fit of the model does more harm than good
    oldnvar <- nvar
    nvar <- sum(active)
    if ((oldLL - oldpenalty > LL - penalty) || (nvar > 1.1* oldnvar)) { retain <- 0.5 * retain     }
  # check convergence
    if(length(checkd)==length(oldcheckd)){
    finishednvar <- (!any(xor(active, oldactive)) &&  all(checkd==oldcheckd))     }
    if(length(checkd)!=length(oldcheckd)){
    finishednvar <- (!any(xor(active, oldactive)) &&  length(checkd)==length(oldcheckd))     }

    finishedn <- !any(xor(active, oldactive))
    if(any(free)){
    finish <- (finishedLL && finishedn && finishedpen) || (all(activdir == 0)) || (iter == maxiter)}
    if(!any(free)){
    finish <- (finishedLL && finishedn && finishedpen) || (all(activdir == 0)) || (iter == maxiter)}

 #  if(iter%%4==1) {
   # print(iter)  
   # print(finishedLL)
   # print(finishedn)
   # print(finishedpen)
  #  print(beta)
  #  print(direc)
# # print(active)
  # print(tryNR)
  #}
  if (!finish) {
      iter <- iter+1
  
      
              if (tryNR){
             NRs<- tryNRstep(beta=beta,direc=direc,active=active,X=X,newfit=newfit,localfit=localfit)
             beta = NRs$beta
             newfit = NRs$newfit
             NRfailed = NRs$NRfailed
                         }
 
 if (!tryNR || NRfailed) {
 
      if (newfit) {
          Xdir <- drop(X[,active,drop=F] %*% activdir)
          if (is.list(localfit$W)) {
            curve <- (sum(Xdir * Xdir * localfit$W$diagW) - drop(crossprod(crossprod(localfit$W$P, Xdir)))) / sum(activdir * activdir)
          } else if (length(localfit$W) > 1) {
            curve <- sum(Xdir * Xdir * localfit$W) / sum(activdir * activdir)
          } else {
            curve <- sum(Xdir * Xdir) / sum(activdir * activdir)
          }
           topt <- 1 / curve
       }
       
  
        # how far can we go in the calculated direction before finding a new zero?
        tedge <- numeric(length(beta))
        tedge[active] <- -activbeta / activdir
        tedge[tedge <= 0] <-tedge[free] <-  2 * topt
        if(length(lambda1)==1 && length(lambda2)==1 && lambda1==0 && lambda2==0){tedgeg=NULL}
        if((length(lambda1)!=1 || length(lambda2)!=1) && (sum(lambda1==0)!=length(lambda1) || sum(lambda2==0)!=length(lambda2))){
        #  tedgeg = -diff(beta)/diff(direc)    ###change to within group difference
		   tedgeg = -Diffw(beta, GRPlen)/Diffw(direc, GRPlen) 
           tedgeg[which(is.na(tedgeg))]=0
         # dchr = diff(chr)                   ###change to within group difference
		  dchr = Diffw(chr, GRPlen)  
          tedgeg[which(dchr!=0)]=0
          tedgeg[tedgeg <= 0] <-  2 * topt
          }


         if(!is.null(tedgeg)){
          mintedge <- min(min(tedge),min(tedgeg))}
          if(is.null(tedgeg)){mintedge <- min(tedge)}
    
      
        # recalculate beta
          if (mintedge + cumsteps < topt) {
          
          beta[active] <- activbeta + mintedge * activdir
          beta[tedge == mintedge] <- 0  # avoids round-off error
          if(!is.null(tedgeg)){
          tmin=which(tedgeg==mintedge)
          beta[tmin]=beta[tmin+1]        ####careful! what does it mean if min is at pos_end? 
           }
          cumsteps=cumsteps+mintedge
          newfit <- (cumsteps > retain * topt) || (nvar == 1)  
          NRfailed <- FALSE
          tryNR <- FALSE

          } else {
          beta[active] <- activbeta + (topt - cumsteps) * activdir
          #if(iter<=1){tryNR <- (cumsteps == 0) && finishednvar && !NRfailed && nvar < m}else{tryNR <- (cumsteps == 0)&& !NRfailed && nvar < n}
          tryNR <- (cumsteps == 0) && !NRfailed && finishednvar && nvar < m
          newfit <- TRUE
          }
         
        }

      }  else {
      converged <- (iter < maxiter)
      } 
	  
     }
 
  return(list(beta = beta, fit = localfit, penalty = c(L1 = penalty1, L2 = penalty2), iterations = iter, converged = converged))
}



######################################################################################
##function for getting the indices of non-zero coefficients and coefficients equal to each other#####
nonZero0 <- function(inp )
{ inp<-round(inp, 9)
 n <- length(inp)

 indices <- matrix(rep(0,n*2),n,2)
 j=1
 i=1
 while(i<=n)
 {
  if(inp[i] == 0)
  {
    i <- i+1
  }
  else if(i<n && (inp[i+1]==inp[i]))
  {
    indices[j,1] <- i
    while(i<n && inp[i+1]==inp[i])
    {
      i <- i+1
    }
    indices[j,2] <- i
    j <- j+1
    i <- i+1
  }
  else
  {
    indices[j,1] <- indices[j,2] <- i
    i <- i+1
    j <- j+1
  }

 }

 return(indices[(1:j-1),])

}
#########################################


 
########################################################################################################################
#Calculate next direction;
########################################################################################################################
getdirec0 <- function(grad,nzb,beta, GRPlen, lambda1,lambda2,positive){
 
    p<-length(beta)
    pos_end<-cumsum(GRPlen)  #last dummy variable positions;
	pos_begin<-c(1, pos_end[-length(pos_end)] +1)
	GRPi<-rep(1:length(GRPlen), GRPlen)
    ########################## Core ########################
	
     direc_i=numeric(p)
     i=1
     vmax = 0
     lam1 <-  lam2 <- P_grad<-cu_g <-numeric(p)
     oldvmax = vmax
     olddirec_i=direc_i
	 
     if(any(nzb)){  #any need to change to within group average? 
              nz <- nonZero0(beta)
              nz=rbind(nz,c(0,0))
	  temp1<-which(nz[,2]-nz[,1]>1)
	  if(length(temp1)>0) {
	 nz1<-NULL
	 for(ik in  temp1 )  
	 { tempnz<-nz[ik,1]:nz[ik,2]  
       temp<- tempnz[which(is.element(tempnz[1:(length(tempnz)-1)], pos_end))]
	   if(length(temp)>0)   
	    {  nz1<-rbind(nz1, c(nz[ik,1], temp[1]) )
		    if(length(temp)==1)  nz1<-rbind(nz1, c(temp+1, nz[ik,2])) else
	       if(length(temp)>1) { for(jk in 1:(length(temp)-1)) nz1<-rbind(nz1,  c(temp[jk]+1, temp[jk+1]))  
		   nz1<-rbind(nz1, c(temp[length(temp)]+1, nz[ik,2])) }
	    }
	 }
     nz2<-rbind(nz1, nz)
     nz<-nz2[!duplicated(nz2[,1]),]
     nz[1:(nrow(nz)-1),]<-nz[order(nz[1:(nrow(nz)-1),1]),]
	                      }
						  
	 for (ik in 1:(nrow(nz)-1)){
              cu_g[nz[ik,1]:nz[ik,2]]=(sum(grad[nz[ik,1]:nz[ik,2]]))
              P_grad[nz[ik,1]:nz[ik,2]]=(cu_g[nz[ik,1]:nz[ik,2]])/(length(nz[ik,1]:nz[ik,2]))
			  
        if(nz[ik,1]==nz[ik,2]){
          lam1[nz[ik,1]]=lambda1[nz[ik,1]]*sign(beta[nz[ik,1]])
				
          if(is.element(nz[ik,1], pos_begin))  lam2[nz[ik,1]]=-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]+1]-beta[nz[ik,1]])*-1
          if(is.element(nz[ik,1], pos_end))    lam2[nz[ik,1]]=-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1])
          if(!is.element(nz[ik,1], c(pos_begin, pos_end)))   lam2[nz[ik,1]]=((-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))-(lambda2[nz[ik,1]]*sign(beta[nz[ik,1]+1]-beta[nz[ik,1]])*-1))
                          }
				 
        if(nz[ik,1]!=nz[ik,2]){  
         lam1[nz[ik,1]:nz[ik,2]]=lambda1[nz[ik,1]:nz[ik,2]]*sign(beta[nz[ik,1]:nz[ik,2]])
		  
		  if(is.element(nz[ik,1], pos_begin)&&is.element(nz[ik,2], pos_end)) lam2[nz[ik,1]:nz[ik,2]]=0		  
          if(is.element(nz[ik,1], pos_begin) && !is.element(nz[ik,2], pos_end)) lam2[nz[ik,1]:nz[ik,2]]=(-lambda2[nz[ik,1]]*sign(beta[nz[ik,2]+1]-beta[nz[ik,1]])*-1)/(length(nz[ik,1]:nz[ik,2]))
          if(!is.element(nz[ik,1], pos_begin)&& is.element(nz[ik,2], pos_end)) lam2[nz[ik,1]:nz[ik,2]]=(-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))/(length(nz[ik,1]:nz[ik,2])) 
          if(!is.element(nz[ik,1], pos_begin)&& !is.element(nz[ik,2], pos_end))  lam2[nz[ik,1]:nz[ik,2]]=(((-lambda2[nz[ik,1]]*sign(beta[nz[ik,1]]-beta[nz[ik,1]-1]))-(lambda2[nz[ik,1]]*sign(beta[nz[ik,2]+1]-beta[nz[ik,1]])*-1)))/(length(nz[ik,1]:nz[ik,2]))
                           }
				              }
                 }
 direc_i= P_grad - lam1 + lam2
 olddirec_i=direc_i


  while(i<p){
       lam1_1 <- lam1_2  <- lam2_1 <- lam2_2 <-cu_g1 <- cu_g2 <- P_grad1 <- P_grad2<-numeric(p)
       direc_i <- numeric(p)
       direc_i= P_grad - lam1 + lam2

        
  ######## Penalized Gradient for nonzero coefficients##################
  if(i<p){

      if(nzb[i]){
       #nzb_m=c(nzb[i:p],F)
	   nzb_m=c(nzb[i: pos_end[(pos_end>=i)][1]],F)
       m=which(diff(nzb_m)==-1)
       m=i+(m-1)
       m=m[1]
                                  #############gives the length of the non-zero stretch######
      if(i==length(nzb)){m=length(nzb)}
      betam=beta[i:m]
        if(sum(diff(betam))==0){j=m}
        if(sum(diff(betam))!=0){
           j=which(diff(betam)!=0)
           j=i+(j-1)
           j=j[1]}
           m=j
               }

    if(!nzb[i]){
            #nzb_mo=c(nzb[i:p],T)
			nzb_mo=c(nzb[i: pos_end[(pos_end>=i)][1]],T)
            mo1=which(diff(nzb_mo)==1)[1]
            mo1=i+(mo1-1)
            m=mo1
            if(is.element(i, pos_end)) {m=i}
              }

      P_grad1=P_grad
      lam1_1=lam1
      lam2_1=lam2
      P_grad1[i:m]=0
      lam1_1[i:m]=0
      lam2_1[i:m]=0

if((length(i:m)==1) && (!nzb[i]) ){

       P_grad1[m]=grad[m]
       lam1_1[m]=lambda1[m]*sign(grad[m])
        if(!is.element(m, pos_end)){
        if(!is.element(m, pos_begin)){
           if(beta[m-1]!=0 & beta[m+1]!=0){
            lam2_1[m]=((-lambda2[m]*sign(beta[m]-beta[m-1]))-(lambda2[m]*sign(beta[m+1]-beta[m])*-1))
           }
           if(beta[m-1]==0 & beta[m+1]!=0){
                lam2_1[m]=((-lambda2[m]*sign(grad[m]-0))-(lambda2[m]*sign(beta[m+1]-beta[m])*-1))
              }
           if (beta[m-1]!=0 & beta[m+1]==0){
                 lam2_1[m]=((-lambda2[m]*sign(beta[m]-beta[m-1]))- (lambda2[m]*sign(0-grad[m])*-1))
               }
                 }

         if(is.element(m, pos_begin)){
           lam2_1[m]=-lambda2[m]*sign(beta[m+1]-beta[m])*-1
           }
          }


       if(is.element(m, pos_end)){
           if(beta[m-1]!=0){
           lam2_1[m]=-lambda2[m]*sign(beta[m]-beta[m-1])
              }
           if(beta[m-1]==0){
           lam2_1[m]=-lambda2[m]*sign(grad[m]-0)
           }
          }
        }


 if(length(i:m)!=1){
     if(nzb[i]) {

      cu_g1=cu_g2=cu_g
       cu_g1[i:m]=cumsum(grad[i:m])
       cu_g2[i:(m-1)]=cu_g1[m]-cu_g1[i:(m-1)]
 
   # ##################################################
      P=numeric(length(cu_g1))
      P[i:m]=1
      P[i:m]=cumsum(P[i:m])
      P_grad1[i:m]=(cu_g1[i:m])/P[i:m]
      P_grad2[i:m]=(cu_g2[i:m])/(P[m]-P)[i:m]
      P_grad2[c(pos_end, m)]=0
    ##########################################

    ###############Calculate lambda1 and lambda2 with signs ################

     if(is.element(i, pos_begin)){
      if (!is.element(m, pos_end)){
       lam2_1[i:m]=((-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1))/P[i:m]
       lam2_1[m]=(-lambda2[m]*sign(beta[m+1]-beta[m])*-1)/(P[m])
       }
       if(is.element(m, pos_end)){
       lam2_1[i:m]=(-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)/P[i:m]
       lam2_1[m]=0
       }
       }

      if(!is.element(i, pos_begin)){
       if(!is.element(m, pos_end)){
         lam2_1[i:m]=((-lambda2[i:m]*sign(beta[i]-beta[i-1])-(lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)))/P[i:m]
         lam2_1[m]=(-lambda2[m]*sign(beta[i]-beta[i-1]) - (lambda2[m]*sign(beta[m+1]-beta[i])*-1))/(P[m])
         }
       if(is.element(m, pos_end)){
         lam2_1[i:m]=((-lambda2[i:m]*sign(beta[i]-beta[i-1])-(lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])*-1)))/P[i:m]
         lam2_1[m]=(-lambda2[m]*sign(beta[i]-beta[i-1]))/(P[m])
         }
         }

     if(!is.element(m, pos_end)){
        lam2_2[i:m]=((-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m])-(lambda2[i:m]*sign(beta[m+1]-beta[i])*-1)))/(P[m]-P[i:m])
        lam2_2[m]=0
        }
     if(is.element(m, pos_end)){
        lam2_2[i:m]=(-lambda2[i:m]*sign(P_grad2[i:m]-P_grad1[i:m]))/(P[m]-P[i:m])
        lam2_2[m]=0
      }
          #####################################################

      lam1_1[i:m] =(lambda1[i:m]*sign(P_grad1[i:m]))
      lam1_2[i:m]=(lambda1[i:m]*sign(P_grad2[i:m]))
      lam1_2[m]=0
      }


#######################################################################

  if(!nzb[i]){

      cu_g1=P_grad
      cu_g1[i:m]=cumsum(grad[i:m])
	  #cu_g1[i:m]=unlist(  tapply(grad[i:m], GRPi[i:m],  cumsum )) #change to within group cumulative length; 
      cu_g2=numeric(length(cu_g1))
   ##################################################
      P=numeric(length(cu_g1))
      P[i:m]=1
      P[i:m]=cumsum(P[i:m])
	 # P[i:m]<- unlist( tapply(P[i:m], GRPi[i:m],  cumsum)  ) #change to within group cumulative length; 
      P_grad1[i:m]=(cu_g1[i:m])/P[i:m]
      P_grad2= numeric(p)
 
	  
    ##########################################

    ###############Calculate lambda1 and lambda2 with signs ################

     if(is.element(i, pos_begin)){
	       if(is.element(m, pos_end)){
                       lam2_1[i:m]= (-lambda2[i:m]*sign(0-P_grad1[i:m])*-1)/P[i:m]
                       lam2_1[m]=0
                       }
           if (!is.element(m, pos_end)){
                         lam2_1[i:m]= (-lambda2[i:m]*sign(0-P_grad1[i:m])*-1)/P[i:m]
                         lam2_1[m]=(-lambda2[m]*sign(beta[m+1]-beta[m])*-1)/(P[m])
                        }
                                }

      if(!is.element(i, pos_begin)){
             if (is.element(m, pos_end)){
                 if(beta[i-1]==0){
                  lam2_1[i:m]=(-lambda2[i:m]*sign(P_grad1[i:m]-0)- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                  lam2_1[m]=(-lambda2[m]*sign(P_grad1[m]-0))/(P[m])
                  }
                 if(beta[i-1]!=0){
                  lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                  lam2_1[m]=(-lambda2[m]*(sign(beta[m]-beta[i-1])))/(P[m])
                  }
                 }
             if(!is.element(m, pos_end)){
                   if(beta[m+1]==0 & beta[i-1]==0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(P_grad1[i:m]-0)))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= lam2_1[m]/P[m]
                   }
                   if(beta[m+1]==0 & beta[i-1]!=0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= lam2_1[m]/P[m]
                   }
                   if(beta[m+1]!=0 & beta[i-1]==0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(P_grad1[i:m]-0)))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m]= ((-lambda2[m]*(sign(P_grad1[m]-0)))- (lambda2[m]*sign(beta[m+1]-beta[m])*-1))/P[m]

                   }
                   if(beta[m+1]!=0 & beta[i-1]!=0){
                   lam2_1[i:m]=((-lambda2[i:m]*(sign(beta[i:m]-beta[i-1])))- (lambda2[i:m]*sign(0-P_grad1[i:m])*-1))/P[i:m]
                   lam2_1[m] = ((-lambda2[m]*(sign(beta[m]-beta[i-1])))- (lambda2[m]*sign(beta[m+1]-beta[m])*-1))/P[m]

                   }
                }

              }
                         lam1_1[i:m]=(lambda1[i:m]*sign(P_grad1[i:m]))

        }


       }

     if(nzb[i] && (length(i:m)!=1)){

      direc1_1=numeric(length(cu_g1))
      direc1_2=numeric(length(cu_g2))

      direc1_1[i:m]=P_grad1[i:m]-lam1_1[i:m]+lam2_1[i:m]
      direc1_1[length(direc1_1)]=P_grad1[length(direc1_1)]-lam1_1[length(direc1_1)]+lam2_1[length(direc1_1)]  #
      direc1_2[i:m]=P_grad2[i:m]-lam1_2[i:m]+lam2_2[i:m]
      test=(direc1_1>direc1_2)
      test1=(direc1_1<direc1_2)
      testd = ((test&((P_grad1)>(P_grad2)))|(test1&((P_grad1)<(P_grad2))))
      testd[m] = F
      if(any(testd[i:m])){
      test2=(i:m)[which(testd[i:m])]

      direc_j=matrix(0,nrow=length(direc1_1),ncol=length(direc1_2))
      if(length(test2)>1){
      for(ik in 1:length(test2)){direc_j[i:test2[ik],test2[ik]]=(direc1_1)[test2[ik]]}
      for(ik in 1:(length(test2)-1)){direc_j[(test2[ik]+1):m,test2[ik]]=(direc1_2)[test2[ik]]}
      if(test2[length(test2)]!=m){ direc_j[(test2[length(test2)]+1):m,test2[length(test2)]]=direc1_2[test2[length(test2)]]   }
       } else
      if(length(test2)==1){
      if (test2==m){
      direc_j[i:test2,test2]=(direc1_1)[test2]
      }
      if (test2!=m){
      direc_j[i:test2,test2]=(direc1_1)[test2]
      direc_j[(test2+1):m,test2]=(direc1_2)[test2]
      }
      }
      norm1 <- apply(direc_j,2,function(il){sqrt(sum(t(il)*il))})
      vmax2=max(norm1)
 
      direc_i[i:m]=direc_j[(i:m),which(norm1==vmax2)]
      if(which(vmax2==norm1)==m){vmax=0}
      if(which(vmax2==norm1)!=m){
      vmax= sqrt(sum(t(direc_i)*direc_i)) }
      }
      }
  if(!nzb[i]){
             if(length(i:m)!=1){
               lam=numeric(length(lam1_1))
               lam=((lam1_1)-(lam2_1))

                if(any((P_grad1[i:m]*sign(P_grad1[i:m]))>(lam[i:m]*sign(P_grad1[i:m])))){
                 d_i = which((P_grad1[i:m]*sign(P_grad1[i:m]))>(lam[i:m]*sign(P_grad1[i:m])))
                 direc1_1= P_grad1 - lam1_1 + lam2_1

                 direc_j=matrix(0,length(i:m),length(i:m))
                 for(im in 1:length(d_i)){
                 direc_j[1:d_i[im],d_i[im]]=((P_grad1[i:m])[d_i[im]]-((lam1_1[i:m])[d_i[im]])+((lam2_1[i:m])[d_i[im]]))    
                        				 }
										 
                 norm_d <- apply(direc_j,2,function(il){sqrt(sum(t(il)*il))})
                 vmax1=max(norm_d)
                 direc1_1[i:m]=direc_j[,which(norm_d==vmax1)]
                 direc_i[i:m]=direc1_1[i:m]

                 vmax = sqrt(sum(t(direc_i)*direc_i))

                 }

              }
             if(length(i:m)==1){

                   if((grad[m]*sign(grad[m]))>(lam1_1[m]-lam2_1[m])*sign(grad[m])){
                      direc_i[m]=grad[m] - lam1_1[m] + lam2_1[m]
                       vmax= sqrt(sum(t(direc_i)*direc_i))
                      }
                   }
              }
# print(i); print(vmax);  print(direc_i); 
                           if (oldvmax>=vmax) {direc_i=olddirec_i
                               vmax=oldvmax}


                      olddirec_i=direc_i
                      oldvmax=vmax

          if(nzb[i]){
                     if(length(i:m)==1){io=i+1}
                     if(length(i:m)!=1){ if(j==m) {io=m+1}
                                         if(j!=m) {io=j+1}  }
                   i=io  } else {i=i+1}
 
         }
		 
} 
    return(list(direc_i = direc_i))

}

 



 
 GroupFusedLasso<-function(GRPlen, lambda1, lambda2)
 {
 source(paste("penalized_0.9-45/penalized/R/checkinput.R", sep=''));
 source(paste("penalized_0.9-45/penalized/R/logit.R", sep=''));
 # source(paste("penalized_0.9-45/penalized/R/core.R", sep=''));
  
myfun1 <- function(response, penalized, unpenalized, lambda1=0, lambda2=0, positive = FALSE, data, fusedl=FALSE,  
  model = c("cox", "logistic", "linear", "poisson"), startbeta, startgamma, steps =1, epsilon = 1e-10, 
  maxiter, standardize = FALSE, trace = TRUE) {
  return(prep<-match.call())
  }
  
 myfun2 <- function(response, penalized, unpenalized, lambda1=0, lambda2=0, positive = FALSE, data, fusedl=FALSE,  
  model = c("cox", "logistic", "linear", "poisson"), startbeta, startgamma, steps =1, epsilon = 1e-10, 
  maxiter, standardize = FALSE, trace = TRUE) {
  return(prep<-parent.frame())
  }
   
 prep1<- myfun1(response=yy, fmla , unpenalized=~0, lambda1=lambda1, lambda2=lambda2, fusedl=TRUE)
 prep2<- myfun2(response=yy, fmla , unpenalized=~0, lambda1=lambda1, lambda2=lambda2, fusedl=TRUE)
  
 prep <- .checkinput(prep1, prep2) 
 fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)$fit
   pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  nr <- nrow(prep$X)
  fusedl <- prep$fusedl
  
  wl1<-1
beta <- prep$beta
lambda1 = lambda1 * wl1 * prep$baselambda1 
lambda2 = lambda2 * prep$baselambda2 
chr = prep$chr 
positive = FALSE
X = prep$X 
trace = FALSE
epsilon = 1e-10 
maxiter = Inf
  
 
out1<-fusedlasso0(beta, GRPlen, chr, lambda1, lambda2,fit, X, positive, trace = FALSE,   epsilon=1e-10, maxiter=Inf) 

 return(out1$beta)
}


DMcmbn<-function(M1, M2)  
{   
   temp<-rep(1:M2, M1)
   temp1<-sort(rep(1:M1, M2),  method='quick')
   pair0<- cbind(temp1, temp)
   pair1 <-pair0[temp>temp1,]
   pair2 <-pair0[temp<temp1,] 
   return(list(pair1=pair1, pair2=pair2))
} 

# pair<-DMcmbn(1:numofvar, 1:numofvar ) 
# pair1<-pair$pair1
# GRPlen<-rep(0, nrow(pair1))
# for(i in 1:nrow(pair1))
# { i1<-pair1[i, 1]
  # i2<-pair1[i, 2]
 # xx1<- paste(dat[, i1], dat[, i2], sep='')
 # GRPlen[i]<-length(table(xx1))
 # } 
 
# fmla<-as.formula(paste("penalized =~ x2_3_ + x1_2_ + x1_3_" ))   
# m1<-  GroupFusedLasso(GRPlen[c(3,1,2)], lambda1=2.5, lambda2=2.5)
# m1
# pen <- penalized(yy, fmla,   unpenalized = ~0, lambda1 =0 , lambda2 = 2.5, fusedl = TRUE) 
# coef(pen,'all')


# length(unique(variab[which(round(coef(pen1,'all'),6)!=0)]))
 
