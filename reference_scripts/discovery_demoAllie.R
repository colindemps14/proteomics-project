library(BiocManager)
BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
library(pcalg) ## install pcalg and related dependencies

set.seed(123)
n <- 10000
X1 <- rnorm(n, 0 , .5)
X4 <- rnorm(n, 0, 1.4)
X2 <- rnorm(n, -1 + 0.8*X1 + 1.2*X4, 1)
X3 <- rnorm(n, -1 + 2.2*X4, 1) ## 1-->2<--4-->3

data <- cbind(X1,X2,X3,X4) ## L is not observed in data

indepTest <- gaussCItest ## specify the independence test
suffStat <- list(C = cor(data), n = n) ## using the correlation matrix

pc.est <- pc(suffStat, indepTest, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est, main="PC Graph")


##changing 
n <- 10000
X1_1 <- rnorm(n, 0 , .9)
X4_1 <- rnorm(n, 0, .8)
X2_1 <- rnorm(n, -1 + 3*X1_1 + .2*X4_1, 1)
X3_1 <- rnorm(n, -1 + 2.2*X2_1, 1) 

data1 <- cbind(X1_1,X2_1,X3_1,X4_1) ## L is not observed in data

indepTest1 <- gaussCItest ## specify the independence test
suffStat1 <- list(C = cor(data1), n = n) ## using the correlation matrix
suffStat1

pc.est1 <- pc(suffStat1, indepTest1, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est1, main="PC Graph 1")

#only show relationships if similar enough in distribution? wasnt showing with crazy difference

n <- 10000
X1_2 <- rnorm(n, 0 , 1)
X4_2 <- rnorm(n, .31*X1_2, 2)
X2_2 <- rnorm(n, -1 + 9*X1_2,1)
X3_2 <- rnorm(n, -2 + 2.2*X1_2,0.5) 

data2 <- cbind(X1_2,X2_2,X3_2,X4_2) ## L is not observed in data

indepTest2 <- gaussCItest ## specify the independence test
suffStat2 <- list(C = cor(data2), n = n) ## using the correlation matrix
suffStat2

pc.est2 <- pc(suffStat2, indepTest2, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est2, main="PC Graph 2")

#why no arrows? not a CPDAG thats why?

n <- 10000
X1_3 <- rnorm(n, 3, 0.8)
X4_3 <- rnorm(n, .31, 2*X1_3)
X2_3 <- rnorm(n, -1 + 2*X4_3,1)
X3_3 <- rnorm(n, -2 + 2.2*X1_3*4*X4_3*2*X2_3,0.5*X4_3) #too many relationships = cancel out?

data3 <- cbind(X1_3,X2_3,X3_3,X4_3) ## L is not observed in data

indepTest3 <- gaussCItest ## specify the independence test
suffStat3 <- list(C = cor(data3), n = n) ## using the correlation matrix
suffStat3

pc.est3 <- pc(suffStat3, indepTest3, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est3, main="PC Graph 3")



n <- 10000
X1_4 <- rnorm(n, 3, 0.8)
X4_4 <- rnorm(n, .31, 2)
X2_4 <- rnorm(n, -1, 1)
X3_4 <- rnorm(n, 2.2*X1_3,0.5) #too many relationships = cancel out?

data4 <- cbind(X1_4,X2_4,X3_4,X4_4) ## L is not observed in data

indepTest4 <- gaussCItest ## specify the independence test
suffStat4 <- list(C = cor(data4), n = n) ## using the correlation matrix
suffStat4

pc.est4 <- pc(suffStat4, indepTest4, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
plot(pc.est4, main="PC Graph")



## WARNING: the default plotting method for CPDAGs in package 'pcalg' will display undirected (--) edges
## as if they are bidirected (<->). This is just a somewhat annoying feature of the plot function, they
## are not really bidirected edges, they are undirected. CPDAGs only contain directed and undirected edges!

