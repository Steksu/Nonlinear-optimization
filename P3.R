options(digits = 10)
setwd("D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt")



f3 <- function(x) {(1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2}

library(rgl) # wczytanie bliblioteki
f3w <-function(x,y) {(1 - x)^2 + 100 * (y - x^2)^2} # podanie funkcji
x <- seq(-1.5, 1.5, length.out = 30) # podanie zakresu
y <- seq(-1, 2, length.out = 30)
z <- outer(x,y, f3w) # obliczenie warto?ci
persp3d(x,y,z,col="lightblue",xlab="x",ylab="y",zlab="z",shade=0.8) # tworzenie wykresu


hookejeeves <- function (f,x,step,alpha,tol) {
  it_p <- 0
  it_r <- 0
  while (step >= tol) {
    xb <- x
    
    x <- trial_stage(f,xb,step)
    it_p = it_p + 1
    #cat("Etap probny wykona? si? poraz: ", it_p, "\n")
    if (f(x) < f(xb)) {
      
      while (f(x) < f(xb)) {
        it_r = it_r + 1
        #cat("Etap roboczy wykona? si? poraz: ", it_r, "\n")
        old.xb <- xb
        xb <- x
        x <- 2*xb-old.xb
        x <- trial_stage(f,x,step)
        it_p = it_p + 1
        #cat("Etap probny wykona? si? poraz: ", it_p, "\n")
      }
      x <- xb
    } else {
      step <- alpha * step
    }
  }
  return(xb)
}
trial_stage <- function (f,x,step) {
  n <- length(x)
  versor <- diag(n)
  i <- 1
  
  while (i <= n) {
    
    if (f(x + step*versor[i,]) < f(x)) {
      x <- x + step*versor[i,]
    } else if (f(x - step*versor[i,]) < f(x)) {
      x <- x - step*versor[i,]
    }
    i <- i + 1
  }
  return(x)
}


external.penalty <- function (f, penalty,x, alphax, tol) {
  i <- 0
  #cat ("Iteracja: ",i,"x",x," f(x)=", f(x),"\n")
  fcelu <- 0
  repeat {
    objective <- function (x) {
      
      fcelu <<- fcelu + 1
      f(x) + alphax * penalty(x)
      
    } 
    
    newx <- hookejeeves(objective,x,0.5,0.5, tol)
    alphax <- 2*alphax
    if (dist(rbind(newx, x)) < tol) {
    
      #cat ("Minimum x = (", newx[1], ",", newx[2], ") o wartości f(x) = ", f(newx),"\n")
      return( c(newx[1], newx[2], f(newx), fcelu))
    }
    i <- i+1
    x <- newx
    #cat("Iteracja", i, "x=", x, "f(x)", f(x), "F(x)=", objective(x), "alpha=", alphax, "\n")
  }
}


internal.penalty <- function (f, penalty,x, alphax, tol) {
  i <- 0
  fcelu <- 0
  alfa <- alphax
  #cat ("Iteracja: ",i,"x",x," f(x)=", f(x),"\n")
  repeat {
    objective <- function (x) {
      
      fcelu <<- fcelu + 1
      f(x) + alphax * penalty(x)
      
      } 
    newx <- hookejeeves(objective,x,0.5,0.5, tol)
    alphax <- alfa*alphax
    if (dist(rbind(newx, x)) < tol) {
      
      #cat ("Minimum x = (", newx[1], ",", newx[2], ") o wartości f(x) = ", f(newx),"\n")
      return( c(newx[1], newx[2], f(newx), fcelu))
    }
    i <- i+1
    x <- newx
    #cat("Iteracja", i, "x=", x, "f(x)", f(x), "F(x)=", objective(x), "alpha=", alphax, "\n")
  }
}


#### Zewnetrzna funkcja kary

a <- 1.9

g1 <- function(x){x[1] + x[2] - a}
g2 <- function(x){(x[1]-1)^3 - x[2] + 1/8}
g3 <- function(x){ -x[1] - 1/2}

# Funkcja zewnętrznej funkcji kary
S_zew <- function(x) {
  g_values <- c(g1(x), g2(x), g3(x))
  return(sum(pmax(0, g_values)^2))
}
    
    
# Funkcja wewnętrznej funkcji kary
S_wew <- function(x) {
  g_values <- c(g1(x), g2(x), g3(x))
  return(-sum(1/g_values))
}



library(readxl)
library(R.utils)
points = read_excel("pointsp3.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), fcelu = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  cat(x_value, y_value)
  result <- withTimeout({
    tryCatch(
      expr = {
        result <- external.penalty(f3, S_zew, c(x_value, y_value), 0.5,  1e-6)
        result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 4))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 25, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 4)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip3.xlsx")





#### Wewnetrzna funkcja kary


a <- 1.9

g1 <- function(x){x[1] + x[2] - a}
g2 <- function(x){(x[1]-1)^3 - x[2] + 1/8}
g3 <- function(x){ -x[1] - 1/2}

# Funkcja zewnętrznej funkcji kary
S_zew <- function(x) {
  g_values <- c(g1(x), g2(x), g3(x))
  return(sum(pmax(0, g_values)^2))
}


# Funkcja wewnętrznej funkcji kary
S_wew <- function(x) {
  g_values <- c(g1(x), g2(x), g3(x))
  return(-sum(1/g_values))
}



library(readxl)
points = read_excel("pointsp3.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), fcelu = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        result <- internal.penalty(f3, S_wew, c(x_value, y_value), 0.5,  1e-6)
        result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 4))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 5, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 4)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip3.xlsx")




