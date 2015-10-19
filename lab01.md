# Многомерное нормальное распределение

## Многомерное нормальное распределение.

### 1. Изучить функцию, позволяющую получать многомерные нормально распределенные выборки $Х = (х_1, \dots, х_N), N >= 3$

Define function for generating a random positive-definite matrix with user-specified positive eigenvalues. If eigenvalues are not specified, they are generated from a uniform distribution.

```r
PDmatrix <- function (n, ev = runif(n, 0, 10)) {
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}
```

Initialize variables

```r
N <- 10^3
mu <- runif(3, 0, 100)
sigma <- PDmatrix(3)

print(mu)
```

```
## [1] 79.12593 31.61143 67.17831
```

```r
print(sigma)
```

```
##            [,1]      [,2]      [,3]
## [1,]  6.4903437 -2.630640 0.9761571
## [2,] -2.6306397  7.609814 1.1830541
## [3,]  0.9761571  1.183054 7.1047546
```

Generate multidimensional sample with init values

```r
library(MASS)
mnorm <- mvrnorm(N, mu, sigma)

head(mnorm)
```

```
##          [,1]     [,2]     [,3]
## [1,] 79.48959 31.02184 66.86238
## [2,] 76.54124 30.59157 64.78002
## [3,] 85.26558 32.75451 69.13405
## [4,] 80.79882 33.06798 66.18413
## [5,] 76.85633 36.52537 66.94344
## [6,] 76.80828 33.95826 66.51868
```

### 2. Выделить две компоненты и построить график соответствующей плотности (маргинальной);
2d:

```r
library(ggplot2)

twocomp <- data.frame(x = mnorm[, 1], y = mnorm[, 2])

ggplot(twocomp, aes(x = twocomp$x, y = twocomp$y)) + geom_point() + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
```

![](lab01_files/figure-html/unnamed-chunk-1-1.png) 

3d:

```r
library(sm)
```

```
## Package 'sm', version 2.2-5.4: type help(sm) for summary information
## 
## Attaching package: 'sm'
## 
## The following object is masked from 'package:MASS':
## 
##     muscle
```

```r
# sm.density(twocomp, display="rgl")
sm.density(twocomp)
```

```
## Warning: weights overwritten by binning
```

```
## Loading required package: rgl
## Loading required package: rpanel
## Loading required package: tcltk
## Package `rpanel', version 1.1-3: type help(rpanel) for summary information
```

![](lab01_files/figure-html/unnamed-chunk-2-1.png) 

Histogram:

```r
xc <- cut(twocomp$x, 20)
yc <- cut(twocomp$y, 20)

z <- table(xc, yc)

library(plot3D)
hist3D(z=z, border = "black")
```

![](lab01_files/figure-html/unnamed-chunk-3-1.png) 

### 3. Выделить три компоненты и построить диаграмму рассеивания scatter diagram (оси Ох1, Ох2, Ох3);

```r
library(scatterplot3d)

# plot3d(mnorm)
scatterplot3d(mnorm, highlight.3d = TRUE, pch = 20)
```

![](lab01_files/figure-html/unnamed-chunk-4-1.png) 

### 4. Оценить параметры норм распределения, например, по формулам из лекции;
Expeсted value

```r
apply(mnorm, 2, mean)
```

```
## [1] 79.03212 31.85334 67.27706
```

Covariance matrix

```r
var(mnorm)
```

```
##            [,1]      [,2]      [,3]
## [1,]  6.1112065 -2.451546 0.9032096
## [2,] -2.4515464  7.313995 1.1753664
## [3,]  0.9032096  1.175366 7.1128551
```

### 5. Сгенерировать выборку из многомерного РВ, отличного от нормального;

```r
mnotnorm <- log(mnorm^c(1, 2, 3))
head(mnotnorm)
```

```
##           [,1]      [,2]      [,3]
## [1,]  4.375626  6.869383 12.607909
## [2,]  8.675659 10.262173  4.170997
## [3,] 13.337313  3.489041  8.472095
## [4,]  4.391962  6.997131 12.577322
## [5,]  8.683876 10.794021  4.203848
## [6,] 13.023937  3.525132  8.394966
```

### 6. Выделить две компоненты и построить соответсвующую гистограмму (убедиться в ненормальности);

2d:

```r
twocomp <- data.frame(x = sample(mnotnorm[, 1]), y = sample(mnotnorm[, 2]))

ggplot(twocomp, aes(x = twocomp$x, y = twocomp$y)) + geom_point() + 
  geom_rug(col=rgb(.5,0,0,alpha=.2))
```

![](lab01_files/figure-html/unnamed-chunk-8-1.png) 

3d:

```r
# sm.density(twocomp, display="rgl")
sm.density(twocomp)
```

![](lab01_files/figure-html/unnamed-chunk-9-1.png) 

Histogram:

```r
xc <- cut(twocomp$x, 20)
yc <- cut(twocomp$y, 20)

z <- table(xc, yc)

library(plot3D)
hist3D(z=z, border = "black")
```

![](lab01_files/figure-html/unnamed-chunk-10-1.png) 

### 7. Для нормальной выборки выбрать 2 компоненты, отличные от первой и вычислить множественный коэф. корреляциии 1-й комп. и выбранных
By formula $$R_{z, xy} = \sqrt{\frac{r_{xz}^2 + r_{yz}^2 - 2 r_{yz} r_{xz} r_{xy}}{1 - r_{xy}^2}}$$ where `z` -- first component.

```r
mcor <- function(d, x, y, z = 1) {
  r_xy <- cor(d[, x], d[, y])
  r_xz <- cor(d[, x], d[, z])
  r_yz <- cor(d[, y], d[, z])
  
  sqrt((r_xz^2 + r_yz^2 - 2 * r_yz * r_xz * r_xy) / (1 - r_xy^2))
}
```

Multiple correlation coefficent:

```r
mcor(mnorm, 2, 3)
```

```
## [1] 0.4174063
```
