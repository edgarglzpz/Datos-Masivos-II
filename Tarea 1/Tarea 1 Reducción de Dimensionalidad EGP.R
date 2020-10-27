#1.Elije una imagen de tu elección y aplica SVD para reducir dimensionalidad
library(imager)
photo_path <- system.file('extdata/SE2.jpeg',package='imager') 
photo <- load.image(photo_path)
grey_photo <- grayscale(photo)
par(mfrow=c(2,3),mar=c(1,1,1,1))
plot(photo)
plot(grey_photo)

svd1<-svd(scale(grey_photo))
svd1
str(svd1)

#matrix multiplication of U,D and V*
#including only first 5 singluar vectors
approx5<-svd1$u[,1:5] %*% diag(svd1$d[1:5]) %*% t(svd1$v[,1:5])
#including only first 10 singular vectors
approx10<-svd1$u[,1:10] %*% diag(svd1$d[1:10]) %*% t(svd1$v[,1:10])
#including only first 20 singular vectors
approx25<-svd1$u[,1:25] %*% diag(svd1$d[1:25]) %*% t(svd1$v[,1:25])
#including only first 20 singular vectors
approx35<-svd1$u[,1:35] %*% diag(svd1$d[1:35]) %*% t(svd1$v[,1:35])
#including only first 20 singular vectors
approx50<-svd1$u[,1:50] %*% diag(svd1$d[1:50]) %*% t(svd1$v[,1:50])
approx88<-svd1$u[,1:88] %*% diag(svd1$d[1:88]) %*% t(svd1$v[,1:88])

par(mfrow=c(2,3),mar=c(1,1,1,1))
#plotting for reduced images
plot(as.cimg(approx5),main="(a) 5 singular Vectors",axes=FALSE)
plot(as.cimg(approx10),main="(b) 10 singular Vectors ",axes=FALSE)
plot(as.cimg(approx35),main="(c) 35 singular Vectors",axes=FALSE)
plot(as.cimg(approx50),main="(c) 50 singular Vectors",axes=FALSE)
plot(as.cimg(approx88),main="(d) 88 singular vectors",axes=FALSE)
plot(as.cimg(grey_photo),main="(d) Full image",axes=FALSE)

#Matriz equivalente

d <- svd1$d
u <- svd1$u
v <- svd1$v

Matrix <- u %*% diag(d) %*% t(v)
#Matrix <- grayscale(photo) %>% as.matrix

library(matrixcalc)
##Calidad de la aproximación
Calidad <- frobenius.norm(Matrix - approx88)/frobenius.norm(Matrix)
Calidad

# PCA

# 3. Calcular la pseudo inversa de las siguientes matrices

# A
X1 <- rbind(c(1,0),c(1,2))

pseudoinversa <- function(M){
  s<-svd(M)
  D <- diag(s$d)
  Dinv <- diag(1/s$d)
  U <- s$u; V <- s$v
  B <- s$v %*% diag(1/s$d) %*% t(s$u)
  return(B)
}

pseudoinversa(X1)

# B

X2 <- rbind(c(1,0,8),c(1,2,-1))

pseudoinversa(X2)


# 4. Resolver los siguientes sistemas de ecuaciones Ax = b.
# A
A <- rbind(c(1,0,8),c(1,2,-1))
b <- cbind(c(2,1))

A_inv=pseudoinversa(A)
A_inv

X1 <- A_inv %*% b
X1

# Comprobacion

A %*% X1

# B
A2 <- cbind(c(1,0),c(1,0))
b2 <- cbind(c(1,1))

A2_inv <- pseudoinversa(A2)
A2_inv

X2<-A2_inv%*%b2
X2
#Buscamos una forma que no devuelva NaN

pseudoinversa2 <- function(M){
  s<-svd(M)
  D <- diag(s$d)
  Dinv <- diag(1/s$d)
  U <- s$u; V <- s$v
  B <- s$v %*% diag(c(1/s$d[1],0)) %*% t(s$u)
  return(B)
}

A2_inv <- pseudoinversa2(A2)
A2_inv
X2_aprox<-A2_inv %*% b2
X2_aprox

# Comprobación
A2 %*% X2_aprox