#' @keywords internal
Estep_expo <- function(tMat,Nmat, Time, G, known) {
  N <- nrow(Time)
  expo <- matrix(NA, ncol = G, nrow = N)
  expov <- matrix(NA, ncol = G, nrow = N)
  for (i in 1:G) {
    qjj <- diag(tMat[, , i])
    tMat_Alt<-tMat[,,i]
    diag(tMat_Alt)<-1
    Nmat_Comp <- rapply(Nmat, function(x) sum(x*log(tMat_Alt)))
    expo[, i] <- apply(Time, 1, function(x) sum(qjj * x))+Nmat_Comp
  }

  if(any(is.infinite(expo))){
    return(list(prob=1))
  } else if (any(is.na(expo))) {
    return(list(prob=1))
  } else{
    #v<--max(expo)
    #v <- -rowMeans(expo)
    #v<-(max(expo)-min(expo))/2
    #v<--median(expo)
    #v <- 0
    v<-rep(0,N)

    for (i in 1:N){

      if (known[i]==0){
        if (max(expo[i,]) < -500){
          diff<--500-max(expo[i,])
          v[i]<-diff
        }
      } else  if (known[i]>0){

        if (expo[i,known[i]]< -500){
          diff<--500-expo[i,known[i]]
          v[i]<-diff
        }
      }
      expo[i,-known[i]]<-1



      expov[i,]<-expo[i,]+v[i]
    }

    return(list(expov = expov, v = v,prob=0))
  }
}
#' @keywords internal
Estep_expo_Dis <- function(tMat,Nmat,G,N,known) {
  expo <- matrix(NA, ncol = G, nrow = N)

  for (i in 1:G) {
    tMat_Alt<-tMat[,,i]
    Nmat_Comp <- rapply(Nmat, function(x) sum(x*log(tMat_Alt)))
    expo[, i] <- Nmat_Comp
  }

  if(any(is.infinite(expo))){
    return(list(prob=1))
  } else if (any(is.na(expo))) {
    return(list(prob=1))
  } else{
    #v<--max(expo)
    #v <- -rowMeans(expo)
    #v<-(max(expo)-min(expo))/2
    #v<--median(expo)
    #v <- 0
    v<-rep(0,N)
    for (i in 1:N){

      if (max(expo[i,]) < -500){
        if (known[i]==0){
          diff<--500-max(expo[i,])
          v[i]<-diff
        } else {
          diff<--500-expo[i,known[i]]
          v[i]<-diff
        }
      }

    }

    expov<-expo
    return(list(expov = expov, v = v,prob=0))
  }
}

#' @keywords internal
E_Step <- function(Nmat, alpha, Pi, tMat, Time, G, Contin, initstate, finalstate, known) {

  N <- nrow(initstate)
  z <- matrix(NA, ncol = G, nrow = N)
  zmat <- matrix(NA, ncol = G, nrow = N)
  #print(known)
  if (Contin) {
    #	for (g in 1:G) {
    #		qjj <- diag(tMat[, , g])
    #		for (i in 1:N) {
    #			z[i, g] <- Pi[g] * sum(alpha[, g] * initstate[i, ]) * exp(sum(qjj * Time[i, ])) * prod(tMat[,
    #				, g]^Nmat[[i]]) * (sum(qjj * finalstate[i, ])) * (-1)
    #		}
    #	}
    expo <- Estep_expo(tMat,Nmat, Time, G,known)

    if(expo$prob==0){

      for (i in 1:G) {
        qjj <- diag(tMat[, , i])

        #for (i in 1:N) {
        #	ztemp <- sum(exp(log(Pi/Pi[g]) + log((alpha[which(initstate[i, ] > 0), ])/(alpha[which(initstate[i,
        #					] > 0), g])) + apply(tMat, 3, function(x) sum(Nmat[[i]] * log(x/tMat[, , g]))) + #apply(tMat,
        #					3, function(x) log(diag(x)[which(finalstate[i, ] == 1)]/diag(tMat[, , g])[which(finalstate[i,
        #						] == 1)])) + apply(tMat, 3, function(x) sum(Time[i, ] * (diag(x) - diag(tMat[, ,
        #					g]))))))
        #			zmat[i, g] <- 1/ztemp
        #		}

        if (G == 1) {
          invzmat <- rowSums(exp(do.call(rbind, replicate(N, log(Pi/Pi[i]), simplify = FALSE)) + (apply(initstate,
                                                                                                    1, function(x) log(alpha[which(x > 0), ]/alpha[which(x > 0), i]))) + do.call(rbind, lapply(Nmat,
                                                                                                                                                                                               function(x) apply(tMat, 3, function(y) sum(x * log(y/tMat[, , i]))))) + (apply(finalstate, 1,
                                                                                                                                                                                                                                                                              function(x) apply(tMat, 3, function(y) log(diag(y)[which(x == 1)]/diag(tMat[, , i])[which(x ==
                                                                                                                                                                                                                                                                                                                                                                          1)])))) + (apply(Time, 1, function(x) apply(tMat, 3, function(y) sum(x * (diag(y) - diag(tMat[,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        , i]))))))))
          zmat[, i] <- 1/invzmat

        } else {
          invzmat <- rowSums(exp(do.call(rbind, replicate(N, log(Pi/Pi[i]), simplify = FALSE)) + t(apply(initstate,
                                                                                                     1, function(x) log(alpha[which(x > 0), ]/alpha[which(x > 0), i]))) + do.call(rbind, lapply(Nmat,
                                                                                                                                                                                                function(x) apply(tMat, 3, function(y) sum(x * log(y/tMat[, , i]))))) + t(apply(finalstate, 1,
                                                                                                                                                                                                                                                                                function(x) apply(tMat, 3, function(y) log(diag(y)[which(x == 1)]/diag(tMat[, , i])[which(x ==
                                                                                                                                                                                                                                                                                                                                                                            1)])))) + t(apply(Time, 1, function(x) apply(tMat, 3, function(y) sum(x * (diag(y) - diag(tMat[,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           , i]))))))))
          zmat[, i] <- 1/invzmat

        }

        alpha_Comp <- apply(initstate, 1, function(x) sum(alpha[, i] * x))
        Final_Comp <- apply(finalstate, 1, function(x) sum(qjj * x))

        z[, i] <- Pi[i] * alpha_Comp * exp(expo$expov[,i]) * (Final_Comp) * (-1)
      }
    }

  } else {
    for (i in 1:G) {

      if(G==1){
        invzmat<-rowSums(exp(do.call(rbind, replicate(N, log(Pi/Pi[i]), simplify = FALSE))+(apply(initstate,
                                                                                              1, function(x) log(alpha[which(x > 0), ]/alpha[which(x > 0), i])))+do.call(rbind, lapply(Nmat,
                                                                                                                                                                                       function(x) apply(tMat, 3, function(y) sum(x * log(y/tMat[, , i])))))))

        zmat[,i]<-1/invzmat
      }
      else{
        invzmat<-rowSums(exp(do.call(rbind, replicate(N, log(Pi/Pi[i]), simplify = FALSE))+t(apply(initstate,
                                                                                               1, function(x) log(alpha[which(x > 0), ]/alpha[which(x > 0), i])))+do.call(rbind, lapply(Nmat,
                                                                                                                                                                                        function(x) apply(tMat, 3, function(y) sum(x * log(y/tMat[, , i])))))))
        zmat[,i]<-1/invzmat


      }



      tMat_Alt<-tMat[,,i]
      diag(tMat_Alt)<-1
      expo <- Estep_expo_Dis(tMat,Nmat, G,N,known)
      z[,i]<-Pi[i]*apply(initstate,1,function(x) sum(alpha[,i]*x))*exp(expo$expov[,i]+expo$v)
    }

  }

  for (n in 1:N){
    tempz<-rep(0,G)
    if (known[n]>0){
      tempz[known[n]]<-1
      zmat[n,]<-tempz
    }
  }




  if(expo$prob==1){
    w.message<-"Inifinite exponential term or NA in exponential"
    return(list(prob=1,w.message))
  }





  if (any(is.na(zmat))) {
    w.message <- "NA in E step"
    #print(w.message)
    return(list(prob = 1, w.message))
  }
  if (any(is.na(z))) {
    w.message <- "NA in z"
    #print(w.message)
    return(list(prob = 1, w.message))
  }

  if (any(is.infinite(z))) {
    w.message <- "Inifinity in z"
    #print(z)
    #print(w.message)
    return(list(prob = 1, w.message))
  }



  return(list(zmat = zmat, z = z, v = expo$v, prob = 0, expo=expo$expov))

}
#' @keywords internal
pig_update <- function(Nmat, Time, zmat, initstate, finalstate, G) {
  #Update pig
  N <- nrow(initstate)
  J <- ncol(initstate)
  Pi <- colSums(zmat)/N
}
#' @keywords internal
alpha_update <- function(Nmat, Time, zmat, initstate, finalstate, G) {
  N <- nrow(initstate)
  J <- ncol(initstate)
  #Update alpha
  alpha <- matrix(NA, nrow = J, ncol = G)
  for (g in 1:G) {
    alpha_temp <- colSums(zmat[, g] * initstate)
    alpha[, g] <- alpha_temp/sum(alpha_temp)
    alpha[which(alpha[,g]<10^-30),g]<-10^-30
  }
  return(alpha)
}
#' @keywords internal
Q_update <- function(Nmat, Time, zmat, initstate, finalstate, G) {
  N <- nrow(initstate)
  J <- ncol(initstate)
  Lambda <- matrix(NA, nrow = J, ncol = G)
  Njmat <- list()
  for (j in 1:J) {
    Njmat[[j]] <- do.call(rbind, lapply(Nmat, function(x) x[j, ]))
  }
  for (g in 1:G) {
    tot_time <- apply(Time, 2, function(x) sum(zmat[, g] * x))
    total_final <- apply(finalstate, 2, function(x) sum(zmat[, g] * x))
    total_trans <- rapply(Njmat, function(x) sum(zmat[, g] * x))
    Lambda[, g] <- (tot_time * total_trans)/(total_final + total_trans)

  }

  Q <- array(NA, dim = c(J, J, G))
  for (g in 1:G) {
    for (j in 1:J) {
      qjjpg <- colSums(zmat[, g] * Njmat[[j]])/Lambda[j, g]
      qjjpg[which(qjjpg < 10^(-30))] <- 10^-30
      qjjpg[j] <- -sum(qjjpg)
      Q[j, , g] <- qjjpg
    }
  }

  return(Q)
}
#' @keywords internal
G_update <- function(Nmat, zmat, initstate, G) {
  N <- nrow(initstate)
  J <- ncol(initstate)
  Gmat <- array(NA, dim = c(J, J, G))
  for (g in 1:G) {
    for (j in 1:J) {
      Njmat <- matrix(NA, nrow = N, ncol = J)
      for (i in 1:N) {
        Njmat[i, ] <- Nmat[[i]][j, ]
      }
      gjjpg <- colSums(zmat[, g] * Njmat)/sum(zmat[, g] * Njmat)
      gjjpg[which(gjjpg<10^-30)] <- 10^-30
      Gmat[j, , g] <- gjjpg

    }
  }
  return(Gmat)
}

#'EM Algorithm for Continuous Time Markov Models
#'
#' This function fits the continuous time first order markov model for a specified set of groups and returns the model chosen by the BIC.
#'
#' @param x A list of states
#' @param t A list of times spent in each state
#' @param J The total number of states
#' @param G A vector containing the number of groups to test
#' @param itemEM The number of emEM iterations for initialization (defaults to 5)
#' @param starts The number of random starting values for the emEM algorithm (defaults to 100)
#' @param maxit The maximum number of iterations after initialization (defaults to 5000)
#' @param tol The tolerance for convergence (defaults to 0.001)
#' @param Contin Fit the continuous time model (defaults to TRUE). If FALSE, fit the discrete model.
#' @param Verbose Display Messages (defaults to TRUE)
#' @param seed Sets the seed for the emEM algorithm (defaults to 1)
#' @param known A vector of labels for semi-supervised classification. 0 indicates unknown observations. The known labels are denoted by their group number (1,2,3, etc.).
#' @param crit The model selection criterion to use ("BIC" or "ICL"). Defaults to "BIC".
#' @param returnall If true, returns the results for all groups considered. Defaults to FALSE.
#'
#' @return Returns a list with parameter and classification estimates for the best model chosen by the selection criterion.
#' @export
#' @examples
#' library(gtools)
#' data(SimData)
#' x<-SimData[[1]]
#' t<-SimData[[2]]
#' Click_2G<-ClickClust_EM(x,t,5,2)

ClickClust_EM<-function(x,t,J,G,itemEM=5,starts=100,maxit=5000,tol=0.001,Contin=TRUE,Verbose=TRUE,seed=1,known=NULL,crit="BIC",returnall=FALSE){
  res<-list()
  it<-0
  for (g in G){
    # print(g)
    it<-it+1
    temp<-tryCatch(ClickClust_Cont_Int(x,t,maxit_em=itemEM,starts=starts,J=J,Contin=Contin,g=g,tol=tol,maxit=maxit,seed=seed,known=known,crit=crit),error=function(e) NULL)
   # print(temp)

    if (is.null(temp)){
      if (crit=="BIC"){
        res[[it]]<-list(BIC=-Inf)
      } else if (crit=="ICL"){
        res[[it]]<-list(ICL=-Inf)
      }
    } else if (temp$prob==1){
      if (crit=="BIC"){
        res[[it]]<-list(BIC=-Inf)
      } else if (crit=="ICL"){
        res[[it]]<-list(ICL=-Inf)
      }
    } else{
      res[[it]]<-temp
    }
  }
  #return(temp)
  #print(res)
  if (crit=="BIC"){
    BIC<-unlist(lapply(res, function(x) x$BIC))
  } else if (crit=="ICL"){
    BIC<-unlist(lapply(res, function(x) x$ICL))
  }
  if (max(BIC)==-Inf){
    if (Verbose==TRUE){
      message("None of the models could be fitted")
    }
    w.message<-"None of the models could be fitted"
    return(list(w.message=w.message,Flag=TRUE))
  } else {
    indexbest<-which.max(BIC)
    if (Verbose==TRUE){
      if (crit=="BIC"){
        message(paste("The best model chosen by the BIC has ",G[indexbest], " groups",sep=""))
      } else if (crit=="ICL"){
        message(paste("The best model chosen by the ICL has ",G[indexbest], " groups",sep=""))
      }
    }
    if (returnall==TRUE){
      Resfin<-res
    } else {
      Resfin<-res[[indexbest]]
    }
    return(Resfin)
  }
}

#' @keywords internal
ClickClust_Cont_Int <- function(x=NULL, t=NULL, init = NULL,
                                maxit_em = 5, starts = 100, J, Contin = TRUE,g=NULL, tol = 0.001, maxit = 5000,seed=1,known=NULL,crit) {

  #Put each sequence of x and t into a dataframe for ease of manipulation
  data <- list()
  N <- length(x)
  if(Contin){
    for (i in 1:N) {
      data[[i]] <- data.frame(x = x[[i]], t = t[[i]])
    }
  } else{
    for (i in 1:N){
      data[[i]]<-data.frame(x=x[[i]])
    }
  }

  if (is.null(known)){
    known<-rep(0,N)
  }

  #Create matrix with number of transistions from one state to the other for each individual
  #Get all permutations of length 2 of the state space
  Perm <- gtools::permutations(J, 2, repeats.allowed = TRUE)
  #Create strings for each permutation
  Perms <- paste(Perm[, 1], Perm[, 2], sep = " ")


  Nmat <- list()

  for (i in 1:N) {
    #Get transisiton counts for individual i
    tabi <- table(factor(paste(head(data[[i]]$x, -1), tail(data[[i]]$x, -1),sep=" "), levels = Perms))
    Nmat[[i]] <- matrix(as.vector(tabi), ncol = J, byrow = TRUE)
  }
  #print(Nmat)


  if (Contin){
    #Create matrix with time spent in each state for each observation
    Time <- matrix(NA, nrow = N, ncol = J)

    for (i in 1:N) {
      temp <- as.vector(tapply(data[[i]]$t, factor(data[[i]]$x, levels = 1:J), sum))
      temp[is.na(temp)] <- 0
      Time[i, ] <- temp
    }
  }

  #Set up initial state indicators
  initial <- matrix(rapply(x, function(x) head(x, 1)), ncol = 1)
  initstate <- t(apply(initial, 1, function(x) My_Indicator(x, J)))

  if(Contin){
    #Set up final state indicators
    final <- matrix(rapply(x, function(x) tail(x, 1)), ncol = 1)
    finalstate <- t(apply(final, 1, function(x) My_Indicator(x, J)))
  }


  #Set up BIC vector
  if (crit=="BIC"){
    BIC <- rep(-Inf, length(g))
    BIC_best <- -Inf
    iter <- 1
  } else if (crit=="ICL"){
    ICL <- rep(-Inf, length(g))
    ICL_best <- -Inf
    iter <- 1
  }

  #Run the EM algorithm for each g


    #Initialize the parameters
    #pig <- rep(1/g, g)



      init <- Init_emEM(Nmat, Time, starts, Contin, maxit_em, initstate, finalstate, g,seed,N,known)
      #print(init)

      zmat<-init[[1]]

    if (init$prob == 1) {
      w.message <- paste("Could not initialize with ", g, "groups", sep = "")
      #print(w.message)
      Initialized <- FALSE
      return(list(w.message,prob=1))
    } else {
      Initialized <- TRUE
    }

    if (Initialized) {

      #print(paste("Initialized for ",g," groups",sep=""))

      #Initialize the group indicators


      Res <- EM_Main(Nmat, Time, zmat, initstate, finalstate, g, Contin, maxit, tol,known)
      if (crit=="BIC"){
        if (Res$prob == 1) {
         # print(paste("Model could not be fitted with ", g, " groups"))
          return(NULL)
        } else {
          ll <- max(Res$likelihood)

          BIC_new <- BIC_Calc(ll, J, g, N)
          BIC[iter] <- BIC_new
        }
        if (BIC_new > BIC_best) {
          BIC_best <- BIC_new
          Final_Res <- Res
          G_best <- g
        }
      } else if (crit=="ICL"){
        if (Res$prob == 1) {
          #print(paste("Model could not be fitted with ", g, " groups"))
          return(NULL)
        } else {
          ll <- max(Res$likelihood)
          ztemp<-Res$zmat
          #print(ztemp)
          BICtemp <- BIC_Calc(ll, J, g, N)
          ICL[iter] <- ICL_Calc(BICtemp,ztemp)
          ICL_new<-ICL[iter]
        }
        if (ICL_new > ICL_best) {
          ICL_best <- ICL_new
          Final_Res <- Res
          G_best <- g
        }
      }



      iter <- iter + 1
    }

  Final_Res$Classification <- imap(Final_Res$zmat)
  Final_Res$G <- G_best
  if (crit=="BIC"){
    Final_Res$BIC <- BIC
  } else if (crit=="ICL"){
    Final_Res$ICL<-ICL
  }
  return(Final_Res)

}
#' @keywords internal
EM_Main <- function(Nmat, Time, zmat, initstate, finalstate, g, Contin, maxit, tol, known) {
  conv <- 0
  wknown<-which(known>0)
  ll <- NULL
  it <- 1
  zmat <- zmat
  while (conv == 0) {

    #Perform the M steps
    #print(zmat)
    #Update Pig
    pig <- pig_update(Nmat, Time, zmat, initstate, finalstate, g)

    #Update Alpha
    alphag <- alpha_update(Nmat, Time, zmat, initstate, finalstate, g)

    #Update Matrix
    if (Contin) {
      Mat <- Q_update(Nmat, Time, zmat, initstate, finalstate, g)
    } else {
      Mat <- G_update(Nmat, zmat, initstate, g)
    }

    #Perform the E step
    zmat_up <- E_Step(Nmat, alphag, pig, Mat, Time, g, Contin, initstate, finalstate, known)
    if (zmat_up$prob == 1) {
      w.message <- paste("Problem with E Step on iteration ", it, sep = "")
      #print(w.message)
      return(list(w.message,prob=1))
    }
    zmat <- zmat_up$zmat
    #Calculate the log likelihood
    #Find max of the row sums of zmat_up$z, call it v
    # Add v to each row of (the row sums of) zmat_up$z
    lltempClust <- log(rowSums(as.matrix(zmat_up$z[which(known==0),])))
    lltempDA <- zmat[which(known>0),]*log(zmat_up$z[which(known>0),])

    ll[it] <- sum(lltempClust)+sum(lltempDA)-sum(zmat_up$v)



    if (is.na(ll[it])) {
      w.message <- "Problem with likelihood calculation"
      #print(w.message)
      return(list(prob = 1, w.message))
    }

    if (it > 3) {
      #Check for decreasing likelihood
      if (ll[it] - ll[it - 1] < 0) {
        w.message <- "Warning! Decreasing Likelihood"
       # print(w.message)
        #print(ll)
        #print(zmat_up$z)
        #print(zmat_up$v)
        return(list(w.message, ll, pig, alphag, Mat, prob = 1,zmat_up$v))
      }

      denomi <- ll[it - 1] - ll[it - 2]

      if (denomi == 0) {
        conv <- 1
      } else {
        ak <- (ll[it] - ll[it - 1])/(ll[it - 1] - ll[it - 2])
        linf <- ll[it - 1] + (ll[it] - ll[it - 1])/(1 - ak)

        if (abs(linf - ll[it - 1]) < tol) {
          conv <- 1
        }

      }

    }
    if (it == maxit) {
      return(list(prob=1))
    }
    it <- it + 1
  }
  #print(ll)

  return(list(pig = pig, alpha = alphag, tMatrix = Mat, zmat = zmat, likelihood = ll, it = it, prob = 0))

}
#' @keywords internal
Init_emEM <- function(Nmat, Time, starts, Contin, maxit_em, initstate, finalstate, g,seed,N,known) {
  set.seed(seed)
  J <- ncol(initstate)
  ll_em <- rep(-Inf, starts)
  ll_best <- -Inf
  Init_Result <- list(prob = 1)

  for (i in 1:starts) {
    zmat <-matrix(sample(1:1000,g*N,replace=TRUE),N,g)
    zmat<-zmat/rowSums(zmat)
    for (n in 1:N){
      tempz<-rep(0,g)
      if (known[n]>0){
        tempz[known[n]]<-1
        zmat[n,]<-tempz
      }
    }
    #print(zmat)


    Res_emEM <- emEM(Nmat, Time, Contin, maxit_em, zmat, initstate, finalstate, g,known)

    if (is.numeric(Res_emEM)) {

      ll_em[i] <- Res_emEM

      if (Res_emEM > ll_best) {
        ll_best <- Res_emEM
        Init_Result <- list(zmat, prob = 0)

      }
    }
  }

  return(Init_Result)
}
#' @keywords internal
emEM <- function(Nmat, Time, Contin, maxit, zmat, initstate, finalstate, g, known) {
  it <- 1
  ll <- NULL
  while (it < (maxit + 1)) {
    #Perform M step
    #Update Pig
    pig <- pig_update(Nmat, Time, zmat, initstate, finalstate, g)


    #Update Alpha
    alphag <- alpha_update(Nmat, Time, zmat, initstate, finalstate, g)

    #Update Q or G

    if (Contin) {
      Mat <- Q_update(Nmat, Time, zmat, initstate, finalstate, g)
    } else {
      Mat <- G_update(Nmat, zmat, initstate, g)
    }

    #Perform E step
    zmat_up <- E_Step(Nmat, alphag, pig, Mat, Time, g, Contin, initstate, finalstate, known)
    zmat <- zmat_up$zmat
    if (zmat_up$prob == 1) {
      w.message <- paste("Problem with E Step on iteration ", it, " in emEM", sep = "")
      #print(w.message)
      return(-Inf)
    }
    #Calculate the log likelihood

    lltempClust <- log(rowSums(as.matrix(zmat_up$z[which(known==0),])))
    #print(sum(lltempClust))
    lltempDA <- rowSums(as.matrix(zmat[which(known>0),])*log(zmat_up$z[which(known>0),]))
    #print(sum(lltempDA))
    #print(sum(zmat_up$v))
    ll[it] <- sum(lltempClust)+sum(lltempDA)-sum(zmat_up$v)
    if (is.na(ll[it])){
      #print(list(Clust=sum(lltempClust),DA=sum(lltempDA),Eup=zmat_up))
      return(-Inf)
    }

    if (is.infinite(ll[it])) {
      w.message <- "Problem calculating likelihood in emEM"
      #print(w.message)
      return(-Inf)
    }
    if (it > 1) {
      if (ll[it] - ll[it - 1] < 0) {
        w.message <- "Decreasing Likelihood in emEM"
        #print(zmat)
        #print(ll)
        #print(w.message)

        return(list(w.message))
      }
    }

    it <- it + 1

  }
  #if(g>1){
  #print(zmat_up$z)
  #}
  return(max(ll))
}
#' @keywords internal
My_Indicator <- function(x, J) {
  ind <- rep(0, J)
  ind[x] <- 1
  ind <- matrix(ind, ncol = J)
  return(ind)
}
#' @keywords internal
BIC_Calc <- function(ll, J, g, N) {
  freepar<-g-1+g*(J-1)+g*J*(J-1)
  BIC <- 2 * ll - freepar* log(N)
}
#' @keywords internal
ICL_Calc<-function(BIC,z){
  ent <- apply(z, 1, max)
  #print(ent)
  #print(BIC)
  ICL_temp <- BIC+2*sum(log(ent))
}
#' @keywords internal
hardzcalc <- function(z) {
  N <- nrow(z)
  tg <- ncol(z)
  hardz <- matrix(NA, ncol = tg, nrow = N)
  for (i in 1:N) {
    temp <- rep(0, tg)
    temp[which.max(z[i, ])] <- 1
    hardz[i, ] <- temp
  }
  return(hardz)
}
#' @keywords internal
imap <- function(z, warn = TRUE) {
  nrowz <- nrow(z)
  zout <- numeric(nrowz)
  J <- 1:ncol(z)
  for (i in 1:nrowz) {
    if (length(unique(z[i, ])) == 1) {
      zout[i] <- sample(1:ncol(z), 1)
    } else {
      zout[i] <- (J[z[i, ] == max(z[i, ])])[1]
    }
  }
  return(zout) #}
}

#'Simulated Data
#'
#' This is a simulated dataset with 2 groups. It is in the form of a list with the first element being the list of states and the second element being the list of time stamps.
#'
#' @docType data
#' @usage data(SimData)
#' @keywords data
"SimData"



