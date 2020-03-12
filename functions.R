spectralClustering <- function (affinity, K, type = 4, 
                                 fast = T, 
                                 maxdim = 50, 
                                 t = 0, neigen = NULL,
                                 kernel = FALSE) 
{
    if(kernel) {
        eps.val <- diffusionMap::epsilonCompute(1-affinity, 0.1)
        L <- exp(-(1-affinity)^2/(eps.val))
        d <- rowSums(L)
        D <- diag(d)
        d[d == 0] <- .Machine$double.eps
    } else {
        d <- rowSums(affinity)
        d[d == 0] <- .Machine$double.eps
        D <- diag(d)
        L <- affinity
    }
    # add kernel???
    neff <- K + 1
    if (type == 1) {
        NL <- L
    }
    else if (type == 2) {
        Di <- diag(1/d)
        NL <- Di %*% L
    }
    else if (type == 3) {
        Di <- diag(1/sqrt(d))
        NL <- Di %*% L %*% Di
    } else if (type == 4) {
        v <- sqrt(d)
        NL <- L/(v%*%t(v))
    }
    if(!fast) {
        eig <- eigen(NL)
    }else {
        f = function(x, A = NULL){ # matrix multiplication for ARPACK
            as.matrix(A %*% x)
        }
        n <- nrow(affinity)
        NL <- as(NL, "dgCMatrix")
        eig <- arpack(f, extra = NL, sym = TRUE,
                      options = list(which = 'LA', nev = neff, n = n, ncv = max(min(c(n,4*neff)))))
    }
    psi = eig$vectors / (eig$vectors[,1] %*% matrix(1, 1, neff))#right ev
    eigenvals <- eig$values
    cat('Computing Spectral Clustering \n')
    res <- sort(abs(eigenvals), index.return = TRUE, decreasing = T)
    U <- eig$vectors[, res$ix[1:K]]
    normalize <- function(x) x/sqrt(sum(x^2))
    if (type == 3 | type == 4) {
        U <- t(apply(U, 1, normalize))
    }
    # This part is equal to performing kmeans
    # labels <- kmeans(U, centers = K, nstart = 1000)$cluster
    eigDiscrete <- .discretisation(U)
    eigDiscrete <- eigDiscrete$discrete
    labels <- apply(eigDiscrete, 1, which.max)
    cat('Computing Diffusion Coordinates\n')
    if (t <= 0) {# use multi-scale geometry
        lambda = eigenvals[-1]/(1 - eigenvals[-1])
        lambda = rep(1,n) %*% t(lambda)
        if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
            lam = lambda[1,]/lambda[1,1]
            # neigen = min(which(lam < .05)) # default number of eigenvalues
            neigen = min(neigen, maxdim, K)
            eigenvals = eigenvals[1:(neigen+1)]  
            cat('Used default value:',neigen,'dimensions\n')
        }
        X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
    }
    else{# use fixed scale t
        lambda = eigenvals[-1]^t
        lambda = rep(1, n) %*% t(lambda)
        if (is.null(neigen)) {#use no. of dimensions corresponding to 95% dropoff
            lam = lambda[1, ]/lambda[1, 1]
            neigen = min(which(lam < .05)) # default number of eigenvalues
            neigen = min(neigen, maxdim)
            eigenvals = eigenvals[1: (neigen + 1)]  
            cat('Used default value:', neigen, 'dimensions\n')
        }
        X = psi[, 2:(neigen + 1)] * lambda[, 1:neigen] #diffusion coords. X
    }
    return(list(labels = labels, 
                eigen_values = eig$values, 
                eigen_vectors = eig$vectors,
                X = X))
}

.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}


.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:1000) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}


getThreshold <- function(mixmdl, verbose = FALSE){
  
  membership <- apply(mixmdl$posterior, 1, which.max)
  m_list <- sort(unique(membership))
  
  mu_list <- mixmdl$mu
  names(mu_list) <- c(1:length(mu_list))
  mu_list <- mu_list[m_list]
  
  if (length(mu_list) > 1) {
    idx1 <- as.numeric(names(mu_list)[order(mu_list)][1])
    idx2 <- as.numeric(names(mu_list)[order(mu_list)][2])
    
    root <- try(uniroot(funMixModel, interval = c(mixmdl$mu[idx1] - mixmdl$sigma[idx1], mixmdl$mu[idx2] + mixmdl$sigma[idx2]),
                        mu1 = mixmdl$mu[idx1], mu2 = mixmdl$mu[idx2],
                        sd1 = mixmdl$sigma[idx1], sd2 = mixmdl$sigma[idx2],
                        rho1 = mixmdl$lambda[idx1], rho2 = mixmdl$lambda[idx2]),
                silent = T)
    
    
    if (class(root) != "try-error") {
      
      t <- root$root
    }else{
      t <- 0
    }
    
  }else{
    t <- 0
  }
  
  return(t)
}


funMixModel <- function(x, mu1, mu2, sd1, sd2, rho1, rho2) {
  
  dnorm(x, mean = mu1, sd = sd1) * rho1 - dnorm(x, mean = mu2, sd = sd2) * rho2
  
  
}
