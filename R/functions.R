#' Estimate the COMET using adaptive thresholding
#'
#' @param dat n by p matrix of multivariate data
#' @param N number of intervals
#' @param mul maximum level of threshold (times of 1 standard error of each sample covariance)
#' @param tol tolerance level of errors
#'
#' @return cov_list1: list of covariance matrix estimators for all thresholds
#' @return aic1: sequence of the values of AIC for all thresholds
#' @return bic1: sequence of the values of BIC for all thresholds
#' @export
#'
#' @examples
#' x = matrix(rnorm(100*10, 0, 1),100,10)
COmet = function(dat, N=100, mul=2, tol=1e-6){
    n = dim(dat)[1]; p = dim(dat)[2]
    S = stats::cov(dat)
    Sig.hat=((n-1)/n)*S
    if (TRUE){
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
    }
    #S1=matrix(diag(Sig.hat), nrow = p, ncol = p)
    #S2=matrix(diag(Sig.hat), nrow = p, ncol = p, byrow = TRUE)
    #var_mat = (Sig.hat^2+S1*S2)

    aictemp_old = 1000000000000; bictemp_old = 1000000000000
    to=mul*N
    cov_list1 = vector("list",to); aicval1 = numeric(to); bicval1 = numeric(to)
    amat_old<-(S!=0)*1; diag(amat_old)=0; A_mit=S
    for (i in 1:to){
        print(i)
        Sig.hat.star=matrix(0,nrow=p,ncol=p)
        am=i/N
        Lam.am=am*sqrt(t.hat*log(p)/n)  #am*se_mat
        #Lam.am=am*sqrt(t.hat)/(n^(0.5-0.01))
        #Lam.am=am*sqrt(var_mat)/(n^(0.5-0.01))
        Sig.hat.star[which(abs(Sig.hat)>Lam.am)]=Sig.hat[which(abs(Sig.hat)>Lam.am)]

        #print(i)
        S_ht <- Sig.hat.star #(hard_thresh(R, lambda_seq[i])!=0) * S
        amat <- S_ht!=0
        amat <- 1*amat; diag(amat) <- 0
        colnames(amat) <- 1:p; rownames(amat) <- 1:p
        if (sum(amat_old)!=sum(amat)){
            A_mit <- covchaud(amat,dat,tol = tol)$mat #covchaud_icf(amat,dat)$mat
            amat_old=amat
        }
        #A_mit <- covchaud_icf(amat,dat)$mat #fitCovGraph(amat, S + cc*diag(p), n=n)$Shat

        cov_list1[[i]] = A_mit

        k=(p+sum(cov_list1[[i]]!=0))/2
        aicval1[i] = n*(log(det(cov_list1[[i]])) + sum(diag(S%*%solve(cov_list1[[i]])))) + 2*k   #sum(cov_list1[[i]]!=0)  # 2*(p+sum(cov_list1[[i]]!=0))/2  # need to minimize + df_mle1[i] #
        bicval1[i] = n*(log(det(cov_list1[[i]])) + sum(diag(S%*%solve(cov_list1[[i]])))) + log(n)*k  #sum(cov_list1[[i]]!=0)  # need to minimize

    }
    return(list("cov_list1"=cov_list1, "aic1"=aicval1, "bic1"=bicval1))
}


covchaud <- function(amat, dat, cc = 0, tol=1e-06){
    n = nrow(dat); p = ncol(dat)
    S = ((n-1)/n)*stats::cov(dat)
    #mat = diag(diag(S))
    #mat = diag(p)
    #mat = S
    mat = S+cc*diag(p)
    matdiff = 1
    while (matdiff > tol){
        mat_old = mat
        for (i in 1:p){
            B = mat[-i,i]*0
            z = dat[,-i] %*% solve(mat[-i,-i])
            #z = demean(dat[,-i]) %*% solve(mat[-i,-i])
            if (sum(amat[i,-i])!=0){
                zsp = z[,amat[i,-i]!=0]
                lsfit2 = solve(crossprod(demean(as.matrix(zsp))) + cc*diag(ncol(as.matrix(zsp)))) %*% crossprod(demean(as.matrix(zsp)), dat[,i])
                #lsfit2 = glmnet(as.matrix(zsp), dat[,i], lambda = cc, alpha = 0, standardize = TRUE, intercept = FALSE)$beta
                B[as.numeric(which(amat[-i,i]!=0))] = lsfit2
            }
            B[as.numeric(which(amat[-i,i]==0))] = 0
            lambda = S[i,i]-t(B)%*%solve(mat[-i,-i])%*%S[-i,i]
            #lambda = S[i,i]+cc-t(B)%*%solve(mat[-i,-i])%*%S[-i,i]
            mat[-i,i] = B
            mat[i,-i] = t(mat[-i,i])
            mat[i,i] = lambda + mat[i,-i] %*% solve(mat[-i,-i]) %*% mat[-i,i]
        }
        matdiff = sum(abs(mat-mat_old))
        #print(matdiff)
    }
    return(list("mat"=mat))
}


demean = function(dat){
    meanmat = matrix(rep(colMeans(dat), nrow(dat)), ncol = ncol(dat), byrow = TRUE)
    datdm = dat - meanmat
    return(datdm)
}
