rm(list=ls(all=TRUE))
library(mvtnorm)
library(mgcov)

#iseed <- iseed
p <-10; n<- 100
num_simu = 1
res_mat = matrix(NA, nrow = num_simu, ncol = 24)
for (ii in 1:num_simu){
    print(ii)
    #seed <- num_simu*(iseed-1)+ii
    seed <- ii
    set.seed(1000*seed)
    Sigma = sigMA(p,c(0.3), hetero = TRUE)                              # Setting 1
    #Sigma = sigMA(p,seq(1-1/(p/2), 0, -1/(p/2)), hetero = TRUE)         # Setting 2
    #Sigma = sigBD(p,0.5,0.5,p/5, hetero = TRUE)                         # Setting 3
    #Sigma = sigBD(p,0.5,0.5,p/2, hetero = TRUE)                         # Setting 4
    ## check     sqrt(1/diag(Sigma)) %*% t(sqrt(1/diag(Sigma))) * Sigma
    dat <- rmvnorm(n,mean=rep(0,ncol(Sigma)),sigma=Sigma)
    S = cov(dat); Sig.hat=((n-1)/n)*S

    thres_res = COmet(dat=dat, N=100, mul=3)
    aic1_est = thres_res$cov_list1[[which.min(thres_res$aic1)]]
    bic1_est = thres_res$cov_list1[[which.min(thres_res$bic1)]]

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

    TPR_aic1 = (sum((aic1_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_bic1 = (sum((bic1_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_aht = (sum((aht_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    TPR_opt = (sum((opt_est!=0)*(Sigma!=0))-p)/(sum(Sigma!=0)-p)
    FPR_aic1 = (sum((aic1_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_bic1 = (sum((bic1_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_aht = (sum((aht_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    FPR_opt = (sum((opt_est!=0)*(Sigma==0)))/(sum(Sigma==0))
    res_mat[ii,1] = TPR_aic1
    res_mat[ii,2] = TPR_bic1
    res_mat[ii,3] = TPR_aht
    res_mat[ii,4] = TPR_opt
    res_mat[ii,5] = FPR_aic1
    res_mat[ii,6] = FPR_bic1
    res_mat[ii,7] = FPR_aht
    res_mat[ii,8] = FPR_opt
    res_mat[ii,9] = sqrt(sum((aic1_est-Sigma)^2))
    res_mat[ii,10] = sqrt(sum((bic1_est-Sigma)^2))
    res_mat[ii,11] = sqrt(sum((aht_est-Sigma)^2))
    res_mat[ii,12] = sqrt(sum((opt_est-Sigma)^2))
    res_mat[ii,13] = sqrt(sum((Sig.hat-Sigma)^2))
    Sig.Inv = solve(Sigma)
    res_mat[ii,14] = -log(det(aic1_est%*%Sig.Inv))+sum(diag(aic1_est%*%Sig.Inv))-p
    res_mat[ii,15] = -log(det(bic1_est%*%Sig.Inv))+sum(diag(bic1_est%*%Sig.Inv))-p
    res_mat[ii,16] = -log(det(aht_est%*%Sig.Inv))+sum(diag(aht_est%*%Sig.Inv))-p
    res_mat[ii,17] = -log(det(opt_est%*%Sig.Inv))+sum(diag(opt_est%*%Sig.Inv))-p
    res_mat[ii,18] = -log(det(Sig.hat%*%Sig.Inv))+sum(diag(Sig.hat%*%Sig.Inv))-p

    res_mat[ii,19] = which.min(thres_res$bic1)/100 #lambda_seq[which.min(thres_res$bic1)]
    res_mat[ii,20] = am
    res_mat[ii,21] = eigen(aht_est)$values[p]
    res_mat[ii,22] = eigen(opt_est)$values[p]

    res_mat[ii,23] = sqrt(sum((mlopt_est-Sigma)^2))
    res_mat[ii,24] = -log(det(mlopt_est%*%Sig.Inv))+sum(diag(mlopt_est%*%Sig.Inv))-p

}

colMeans(res_mat)
