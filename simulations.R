rm(list=ls(all=TRUE))
library(mvtnorm)
library(mgcov)

##############################################################
# R functions for generating correlation matrix #
##############################################################
### AR model ###
sigAR = function(p, rho, hetero=NULL){
    A = diag(rep(1, p))
    for (i in 1 : (p - 1)){
        for (j in (i + 1) : p){
            A[i, j] = rho^(abs(i - j))
            A[j, i] = A[i, j]
        }
    }
    if (!is.null(hetero)){
        D = runif(p, 0.1, 10)
        A = diag(D) %*% A %*% diag(D)
    }
    return(A)
}
### Block model ###
sigBD = function(p, a, b, k, hetero=NULL){
    # a and b are the limits of the uniform distribution and k is the block size
    A = diag(rep(1, p))
    m = p / k
    for (h in 1 : m){
        for (i in ((h - 1) * k + 1) : (h * k)){
            for (j in ((h - 1) * k + 1) : (h * k)){
                if (i > j){
                    temp = runif(1, a, b)
                    A[i, j] = temp
                    A[j, i] = temp
                }
            }
        }
    }
    if (!is.null(hetero)){
        D = runif(p, 0.1, 10)
        #D = runif(p, 2, 5000)
        A = diag(D) %*% A %*% diag(D)
    }
    return(A)
}
### MA model ###
sigMA = function(p, offdiag, hetero=NULL){
    # a and b are the limits of the uniform distribution and k is the block size
    A = diag(rep(1, p))
    for (i in 1:p){
        for (j in 1:p){
            for (k in 1:length(offdiag)){
                if (abs(i-j)==k){A[i,j] <- offdiag[k]}
            }
            #if (abs(i-j)==1){A[i,j] <- a}
        }
    }
    if (!is.null(hetero)){
        D = runif(p, 0.1, 10)
        A = diag(D) %*% A %*% diag(D)
    }
    return(A)
}


##############################################################
# R function for the closed-form thresholding (Qiu 2019) #
##############################################################
opt_select = function(dat){
    X = dat; n = dim(X)[1]; p = dim(X)[2]
    a=min(sqrt(2 + log(n) / log(p)),2) # Defining a for the intervals S1 and S2.
    indexM = matrix(1, p, p) - diag(rep(1, p))
    Sn = ((n-1)/n)*cov(X)
    Xbar = colSums(X)/n
    # Evaluating theta hat using the estimator proposed by Cai and Liu (2011)
    thetahat = matrix(0, p, p)
    for(i in 1 : p){
        for(j in 1 : p){
            s = c()
            for(k in 1 : n){
                s[k] = (((X[k, i] - Xbar[i]) * (X[k, j] - Xbar[j])) - Sn[i, j])^2
            }
            thetahat[i, j] = sum(s) / n
        }
    }
    thetahatDiag = thetahat * indexM

    a0 =sqrt(1/log(log(p))) # Contamination correction factor
    # Evaluating the proposed estimator
    temp = c1c2Est(n, Sn, thetahat, a1 = 2 - a+a0, a2 = 2 , c0 = sqrt(log(p)))
    deltaProp = temp$deltaStarhat
    print(deltaProp)
    # Evaluating the adaptive thresholding covariance estimator using the proposed delta
    SigmaProp = matrix(0, p, p)
    SigmaProp[which(abs(Sn) > deltaProp * sqrt(thetahatDiag * log(p) / n), arr.ind = TRUE)] = Sn[which(abs(Sn) > deltaProp * sqrt(thetahatDiag * log(p) / n), arr.ind = TRUE)]

    return(SigmaProp)

}
c1c2Est = function(n, Sn, t.hat, a1 , a2, c0 ){
    # Sn is the sample covariance matrix.
    # t.hat is the matrix of theta^_j1j2.
    # c0 prevents getting negative values for N2hat.
    p = dim(Sn)[1]
    ALamhat1 = c()
    kk = 0
    for (i in 1 : (p-1 )){
        for (j in (i+1 ) : p){
            kk = kk + 1
            ALamhat1[kk] = abs(Sn[i, j] / sqrt(t.hat[i, j]) * sqrt(n / log(p)))
        }
    }
    Mk = sum((ALamhat1 > a1) & (ALamhat1 < a2))
    nhat = Mk - 2 * (pnorm((a2) * sqrt(log(p))) - pnorm((a1)* sqrt(log(p)))) * ((p^2 - p) / 2)
    if(nhat > c0)  {on = (log(nhat / sqrt(log(p)))) / log(p); deltaStarhat = sqrt(2 * (2 - on))}
    if(nhat <= c0)  {deltaStarhat = 2}
    list(deltaStarhat = deltaStarhat)
}


#####################################################
#R function for selecting delta using the CV method #
#####################################################
find.del=function(n,p,H,N,X,mul=2,seed=NULL){
    # H is the number of splits of the data
    # 2N is the number of points in the partition of the interval [0,2]
    # X is the data matrix
    set.seed(seed)
    nn1=ceiling(n*(1-(1/log(n))))
    sub=c()
    Ra.ij=array(0,dim=c(H,mul*N))
    for (l in 1:H){
        sub=sample(nrow(X),nn1,replace = FALSE, prob = NULL)
        sam1=X[sub,]
        sam2=X[-sub,]
        n1=nrow(sam1)
        n2=nrow(sam2)
        Sig.hat1=((n1-1)/n1)*cov(sam1)
        Sig.hat2=((n2-1)/n2)*cov(sam2)
        one=matrix(1,nrow=n1,ncol=n1)
        Xbar=one%*%sam1/n1
        t.hat=array(0,dim=c(p,p))
        that1=c()
        for( i in 1:p){
            for(j in 1:p){
                s=c()
                for(k in 1:n1){
                    s[k]=((sam1[k,i]-Xbar[k,i])*(sam1[k,j]-Xbar[k,j])-Sig.hat1[i,j])^2
                }
                that1[j]=sum(s)/n1
            }
            t.hat[i,]=that1
        }
        Raj.i=c()
        to=mul*N
        for (m in 1:to){
            Sig.hat.star1=matrix(0,nrow=p,ncol=p)
            am=m/N
            Lam.am=am*sqrt(t.hat*log(p)/n1)
            Sig.hat.star1[which(abs(Sig.hat1)>Lam.am)]=Sig.hat1[which(abs(Sig.hat1)>Lam.am)]
            Raj.i[m]=sum((Sig.hat.star1-Sig.hat2)^2)
        }
        Ra.ij[l,]=Raj.i
    }

    aj=c()
    Rhat.a=c()
    for (j in 1:to){
        Rhat.a=c(Rhat.a,sum(Ra.ij[,j])/H)
        aj[j]=j/N
    }
    j.hat=which(Rhat.a==min(Rhat.a))

    delta.hat=aj[min(j.hat)]
    return(delta.hat)
}


######################################################
################## Simulation codes ##################
######################################################
p <-10; n<- 100
num_simu = 100
res_mat = matrix(NA, nrow = num_simu, ncol = 24)
for (ii in 1:num_simu){
    print(ii)
    seed <- ii
    set.seed(1000*seed)
    Sigma = sigMA(p,c(0.3), hetero = TRUE)                              # MA model
    #Sigma = sigAR(p,1,0.5, hetero = TRUE)                              # AR model
    #Sigma = sigBD(p,1,0.5,0.5,p/5, hetero = TRUE)                      # Block model
    ## check     sqrt(1/diag(Sigma)) %*% t(sqrt(1/diag(Sigma))) * Sigma
    dat <- rmvnorm(n,mean=rep(0,ncol(Sigma)),sigma=Sigma)
    S = cov(dat); Sig.hat=((n-1)/n)*S

    thres_res = COmet(dat=dat, N=100, mul=3)
    aic_est = thres_res$cov_list[[which.min(thres_res$aic)]]
    bic_est = thres_res$cov_list[[which.min(thres_res$bic)]]

    one=matrix(1,nrow=n,ncol=n)
    Xbar=one%*%dat/n
    t.hat=array(0,dim=c(p,p))
    that1=c()
    for( i in 1:p){
        for(j in 1:p){
            s=c()
            for(k in 1:n){
                s[k]=((dat[k,i]-Xbar[k,i])*(dat[k,j]-Xbar[k,j])-Sig.hat[i,j])^2
            }
            that1[j]=sum(s)/n
        }
        t.hat[i,]=that1
    }
    am <- find.del(n,p,H=100,N=100, X=dat, mul=3, seed = 1000000)
    Lam.am=am*sqrt(t.hat*log(p)/n)
    Sig.hat.star=matrix(0,nrow=p,ncol=p)
    Sig.hat.star[which(abs(Sig.hat)>Lam.am)]=Sig.hat[which(abs(Sig.hat)>Lam.am)]
    aht_est <- Sig.hat.star

    opt_est = opt_select(dat)
    amat <- opt_est!=0
    amat <- 1*amat; diag(amat) <- 0
    colnames(amat) <- 1:p; rownames(amat) <- 1:p
    mlopt_est <- covchaud(amat,dat)$mat

    TPR_aic = (sum((aic_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_bic = (sum((bic_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_aht = (sum((aht_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_opt = (sum((opt_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    FPR_aic = (sum((aic_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_bic = (sum((bic_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_aht = (sum((aht_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_opt = (sum((opt_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    res_mat[ii,1] = TPR_aic
    res_mat[ii,2] = TPR_bic
    res_mat[ii,3] = TPR_aht
    res_mat[ii,4] = TPR_opt
    res_mat[ii,5] = FPR_aic
    res_mat[ii,6] = FPR_bic
    res_mat[ii,7] = FPR_aht
    res_mat[ii,8] = FPR_opt
    res_mat[ii,9] = sqrt(sum((aic_est-Sigma)^2))
    res_mat[ii,10] = sqrt(sum((bic_est-Sigma)^2))
    res_mat[ii,11] = sqrt(sum((aht_est-Sigma)^2))
    res_mat[ii,12] = sqrt(sum((opt_est-Sigma)^2))
    res_mat[ii,13] = sqrt(sum((Sig.hat-Sigma)^2))
    Sig.Inv = solve(Sigma)
    res_mat[ii,14] = -log(det(aic_est%*%Sig.Inv))+sum(diag(aic_est%*%Sig.Inv))-p
    res_mat[ii,15] = -log(det(bic_est%*%Sig.Inv))+sum(diag(bic_est%*%Sig.Inv))-p
    res_mat[ii,16] = -log(det(aht_est%*%Sig.Inv))+sum(diag(aht_est%*%Sig.Inv))-p
    res_mat[ii,17] = -log(det(opt_est%*%Sig.Inv))+sum(diag(opt_est%*%Sig.Inv))-p
    res_mat[ii,18] = -log(det(Sig.hat%*%Sig.Inv))+sum(diag(Sig.hat%*%Sig.Inv))-p

    res_mat[ii,19] = which.min(thres_res$bic)/100 #lambda_seq[which.min(thres_res$bic)]
    res_mat[ii,20] = am
    res_mat[ii,21] = eigen(aht_est)$values[p]
    res_mat[ii,22] = eigen(opt_est)$values[p]

    res_mat[ii,23] = sqrt(sum((mlopt_est-Sigma)^2))
    res_mat[ii,24] = -log(det(mlopt_est%*%Sig.Inv))+sum(diag(mlopt_est%*%Sig.Inv))-p

}

avg = colMeans(res_mat)
# TPR / FPR
tb = rbind(avg[1:4], avg[5:8]); rownames(tb)=c("TPR","FPR"); colnames(tb)=c("AIC","BIC","CV","CF"); tb
# Frobenius loss
frob=res_mat[,c(13,9,10,23,11,12)]; colnames(frob)=c("S","AIC","BIC","CCF","CV","HCF"); boxplot(frob)
# Entropy loss
ent=res_mat[,c(18,14,15,24,17,18)]; colnames(ent)=c("S","AIC","BIC","CCF","CV","HCF"); boxplot(ent)
# Positive-definiteness
sum(res_mat[,21]>0); sum(res_mat[,22]>0)
