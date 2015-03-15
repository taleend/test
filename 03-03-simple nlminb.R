x <- rnbinom(100, mu = 10, size = 10)
hdev <- function(par) 
{
    -sum(dnbinom(x, mu = par[1], size = par[2], log = TRUE))
}
e1 <- nlminb(c(9, 12), hdev)$par
e2 <- nlminb(c(20, 20), hdev, lower = 0, upper = Inf)$par
e3 <-nlminb(c(20, 20), hdev, lower = 0.001, upper = Inf)$par