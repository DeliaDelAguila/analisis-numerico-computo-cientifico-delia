
# SVD function
r_svd <- function(A) {
  
  A <- as.matrix(A)
  svd <- list()
  n <- nrow(A)
  m <- ncol(A)
  
  # Valores singulares de A
  AtA <- t(A)%*%A
  sv <- as.matrix(sqrt(abs(eigen(AtA)$values)))
  
  # Calculo de S
  S <- diag(n)*array(sv,c(n,n))
  svd$S <- S
  
  # Eigenvectores y V
  V <- eigen(AtA)$vectors
  svd$V <- V
  
  # Calculo de U
  U <- A %*% V %*% solve(S)
  svd$U <- U
  
  svd
}

# Prueba SVD

A <- cbind(c(1,2,3),c(8,9,10),c(17,38,29)); A
A_svd <- r_svd(A); A_svd$U %*% A_svd$S %*% t(A_svd$V)






library(Matrix)
A <- Hilbert(8); A
A_svd <- r_svd(A); A_svd$U %*% A_svd$S %*% t(A_svd$V)

# Minimos cuadrados

r_mc <- function(A,b) {
  A_svd <- r_svd(A)
  A_svd$V %*% solve(A_svd$S) %*% t(A_svd$U) %*% b
}

# Prueba Minimos cuadrados

A <- cbind(c(1,0,2),c(1,2,5),c(1,5,-1))
x <- c(5,3,-2)
b <- A%*%x
xm <- r_mc(A,b)

A <- Hilbert(8)
x <- c(2, 3, 5, 7, 11, 13, 17, 19)
b <- A%*%x
xm <- r_mc(A,b)

A%*%x
A%*%xm


  
  