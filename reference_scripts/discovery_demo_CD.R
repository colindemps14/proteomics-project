library(pcalg) ## install pcalg and related dependencies
library(Rgraphviz)

#DAG 0
set.seed(123)
n <- 10000
X1 <- rnorm(n, 0 , .5)
X4 <- rnorm(n, 0, 1.4)
X2 <- rnorm(n, -1 + 0.8*X1 + 1.2*X4, 1)
X3 <- rnorm(n, -1 + 2.2*X4, 1) ## 1-->2<--4-->3

#DAG 1
#chain 
set.seed(123)

n <- 10000

X1 <- rnorm(n, 0, 1)

X2 <- rnorm(n, 1 + 1.5*X1, 1)

X3 <- rnorm(n, -1 + 0.8*X2, 1)

X4 <- rnorm(n, 0.5 + 1.2*X3, 1)

#DAG 2
#fork
set.seed(123)
n <- 10000

X1 <- rnorm(n)

X2 <- rnorm(n, 0.8*X1, 1)

X3 <- rnorm(n, 1.5*X1, 1)

X4 <- rnorm(n, -1 + 2*X1, 1)

#DAG 3

#Diamond
set.seed(123)
n <- 10000

X1 <- rnorm(n)

X2 <- rnorm(n, 1 + X1, 1)

X3 <- rnorm(n, -1 + 1.5*X1, 1)

X4 <- rnorm(n, 0.8*X2 + 1.2*X3, 1)

#DAG 4
set.seed(123)
n <- 10000

# Root nodes
X1 <- rnorm(n, 0, 1)
X2 <- rnorm(n, 0, 1)

# Collider
X3 <- rnorm(n, 1 + 1.2*X1 + 0.8*X2, 1)

# Child
X4 <- rnorm(n, -1 + 1.5*X3, 1)

data <- cbind(X1,X2,X3,X4) ## L is not observed in data

indepTest <- gaussCItest ## specify the independence test
suffStat <- list(C = cor(data), n = n) ## using the correlation matrix

pc.est <- pc(suffStat, indepTest, alpha = 0.05, p = 4, verbose=TRUE) ## estimate at CPDAG
Title = "PC Graph 4"
plot(pc.est, main = Title)

## WARNING: the default plotting method for CPDAGs in package 'pcalg' will display undirected (--) edges
## as if they are bidirected (<->). This is just a somewhat annoying feature of the plot function, they
## are not really bidirected edges, they are undirected. CPDAGs only contain directed and undirected edges!

