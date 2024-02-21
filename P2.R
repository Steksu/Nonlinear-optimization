install.packages('writexl')
install.packages('R.utils')
library(R.utils)

setwd("D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt")

f2 <- function(x) {
  x[1]^2 - x[2]^2 - cos(1.8 * pi * x[1]) - cos(1.6 * pi * (x[2] - 0.5)) + 3
}

result <- optim(c(0, 0), f2, control = list(fnscale = -1))

cat("Global Minimum Value:", result$value, "\n")
cat("Optimal Parameters:", result$par, "\n")

mf2 <- function(x) {
  (-1)*(x[1]^2 - x[2]^2 - cos(1.8 * pi * x[1]) - cos(1.6 * pi * (x[2] - 0.5)) + 3)
}

library(rgl) # wczytanie bliblioteki
f1w <-function(x,y) {x^2 - y^2 - cos(1.8 * pi * x) - cos(1.6 * pi * (y - 0.5)) + 3} # podanie funkcji
x <- seq(-10, 10, length.out = 300) # podanie zakresu
y <- seq(-10, 10, length.out = 300)
z <- outer(x,y, f1w) # obliczenie warto?ci
persp3d(x,y,z,col="lightblue",xlab="x",ylab="y",zlab="z",shade=0.8) # tworzenie wykresu


golden <- function(f,lower, upper, tol){
  
  alpha <- (sqrt(5)-1)/2
  n <- ceiling(log(2*tol/(upper-lower),alpha))
  iter = 0
  x1 <- alpha * lower + (1-alpha) * upper
  f.x1 <- f(x1)
  x2 <- (1-alpha) * lower + alpha * upper
  f.x2 <- f(x2)
  while(abs(upper - lower) > 2 * tol){
    if (f.x1 < f.x2){
      upper <- x2
      x2 <- x1
      f.x2 <- f.x1
      x1 <- alpha * lower + (1-alpha) * upper
      f.x1 <- f(x1)
    } else {
      lower <- x1
      x1 <- x2
      f.x1 <- f.x2
      x2 <- (1-alpha) * lower + alpha * upper
      f.x2 <- f(x2)
    }
    iter = iter + 1
    #cat(" iteracja: ", iter, "a= ",lower, "b= ", upper, "\n")
  }
  #cat("Liczba iteracji", n, "\n")
  return((lower + upper) / 2)
}

library(numDeriv)
cauchy.golden <- function(f, x, tol) {
  iteracja <- 1
  cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
  repeat {
    g <- function(a) {f(x - a * grad(f,x))}
    step <- golden(g,0,5,tol)
    new.x <- x - step * grad(f, x)
    iteracja = iteracja + 1
    cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
    if (dist(rbind(new.x,x)) < tol) {
      return(new.x)
    }
    x <- new.x
  }
}

cauchy.golden(f2,c(1,1),0.001)


golden_max <- function(f,lower, upper, tol){
  
  alpha <- (sqrt(5)-1)/2
  n <- ceiling(log(2*tol/(upper-lower),alpha))
  iter = 0
  x1 <- alpha * lower + (1-alpha) * upper
  f.x1 <- f(x1)
  x2 <- (1-alpha) * lower + alpha * upper
  f.x2 <- f(x2)
  while(abs(upper - lower) > 2 * tol){
    if (f.x1 > f.x2){
      upper <- x2
      x2 <- x1
      f.x2 <- f.x1
      x1 <- alpha * lower + (1-alpha) * upper
      f.x1 <- f(x1)
    } else {
      lower <- x1
      x1 <- x2
      f.x1 <- f.x2
      x2 <- (1-alpha) * lower + alpha * upper
      f.x2 <- f(x2)
    }
    iter = iter + 1
    #cat(" iteracja: ", iter, "a= ",lower, "b= ", upper, "\n")
  }
  #cat("Liczba iteracji", n, "\n")
  return((lower + upper) / 2)
}


library(numDeriv)
cauchy.golden.max <- function(f, x, tol) {
  iteracja <- 1
  cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
  repeat {
    g <- function(a) {f(x + a * grad(f,x))}
    step <- golden_max(g,0,5,tol)
    new.x <- x + step * grad(f, x)
    iteracja = iteracja + 1
    cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
    if (dist(rbind(new.x,x)) < tol) {
      return(new.x)
    }
    x <- new.x
  }
}

cauchy.golden(f2,c(1,1),0.001)

hestenes.stiefel <- function(f, x, tol) {
  beta <- 1
  d <- -grad(f,x)
  iteracja <- 1
  cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x),"\n")
  repeat {
    g <- function(a) {f(x + a * d)}
    step <- 0.15
    new.x <- x + step * d
    iteracja = iteracja + 1
    cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
    if (dist(rbind(new.x,x)) < tol) {
      return(new.x)
    }
    beta <- t(grad(f,new.x)) %*% (grad(f,new.x) - grad(f,x)) / (t(d) %*% (grad(f,new.x) - grad(f,x)))
    d <- -grad(f,new.x) + as.vector(beta) * d
    x <- new.x
  }
}


f2(c(-0.5929943,-1))

library(matlib)

bfgs <- function(f, x, tol) {
  
  n <- length(x)
  v <- diag(n)
  g <- function(a) {f(x-a*grad(f,x))}
  
  step <- golden(g,0,5,tol)
  new.x <- x - step * grad(f, x)
  
  iteracja <- 1
  repeat {
    
    cat(" Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
    s <- grad(f,new.x) - grad(f,x)
    r <- new.x - x
    ma <- t(r) %*% s
    a <- (r %*% t(r))/ma[1,1]
    mc <- ma
    c <- (t(s) %*% v %*% s)/mc[1,1]
    md <- ma
    d <- -(r %*% t(s) %*% v + v %*% s %*% t(r))/md[1,1]
    v <- v + (1+as.vector(c)) * a + d
    x <- new.x
    new.x <- x - v %*% grad(f,x)
    new.x <- new.x[,1]
    dist <- dist(rbind(new.x,x))
    if (dist(rbind(new.x, x)) < tol) {
      
      cat("Iteracja", iteracja+1, "x = ", new.x, "f(x)= ", f(new.x), "\n")
      cat("Macierz V= ", v, "\n")
      cat("Hesjan H^-1= ",inv(hessian(f,new.x)), "\n" )
      return(new.x)
    }
    iteracja <- iteracja + 1
  }
}


bfgs(mf2, c(-0.093331688,	-0.3975422), 0.001)
f2(c(-1.787710, -1.269643))


#### BFGS stałokroowka
bfgs <- function(f, x, tol) {
  
  result <- tryCatch(
    expr = {
      n <- length(x)
      v <- diag(n)
      g <- function(a) {f(x-a*grad(f,x))}
      step <- 0.02
      new.x <- x - step * grad(f, x)
      grad_count <- 1
      iteracja <- 1
      repeat {
        
        #cat(" Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
        s <- grad(f,new.x) - grad(f,x)
        grad_count <- grad_count + 2
        r <- new.x - x
        ma <- t(r) %*% s
        a <- (r %*% t(r))/ma[1,1]
        mc <- ma
        c <- (t(s) %*% v %*% s)/mc[1,1]
        md <- ma
        d <- -(r %*% t(s) %*% v + v %*% s %*% t(r))/md[1,1]
        v <- v + (1+as.vector(c)) * a + d
        x <- new.x
        new.x <- x - v %*% grad(f,x)
        grad_count <- grad_count + 1
        new.x <- new.x[,1]
        dist <- dist(rbind(new.x,x))
        if (dist(rbind(new.x, x)) < tol) {
          
          #cat("Iteracja", iteracja+1, "x = ", new.x, "f(x)= ", f(new.x), "\n")
          #cat("Macierz V= ", v, "\n")
          #cat("Hesjan H^-1= ",inv(hessian(f,new.x)), "\n" )
          return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
        }
        iteracja <- iteracja + 1
      }
    },
    error = function(e) {
      # Error handling code
      cat("An error occurred: ", conditionMessage(e), "\n")
      return(c(NA,NA,NA,iteracja, grad_count))  # Value to be returned in case of an error
    }
  )
  
  return(result)
  
  
 
}

# Wczytanie losowych punktów do ramki danych
library(readxl)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), iteracja = numeric(), grad_count = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        bfgs_result <- bfgs(f2, c(x_value, y_value), 1e-9)
        bfgs_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)


library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")


#### BFGS zmiennokrokowa

bfgs <- function(f, x, tol) {
  n <- length(x)
  v <- diag(n)
  g <- function(a) {f(x-a*grad(f,x))}
  step <- golden(g,0,0.1,tol)
  new.x <- x - step * grad(f, x)
  grad_count <- 1
  iteracja <- 1
  repeat {
    
    #cat(" Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
    s <- grad(f,new.x) - grad(f,x)
    grad_count <- grad_count + 2
    r <- new.x - x
    ma <- t(r) %*% s
    a <- (r %*% t(r))/ma[1,1]
    mc <- ma
    c <- (t(s) %*% v %*% s)/mc[1,1]
    md <- ma
    d <- -(r %*% t(s) %*% v + v %*% s %*% t(r))/md[1,1]
    v <- v + (1+as.vector(c)) * a + d
    x <- new.x
    new.x <- x - v %*% grad(f,x)
    grad_count <- grad_count + 1
    new.x <- new.x[,1]
    dist <- dist(rbind(new.x,x))
    if (dist(rbind(new.x, x)) < tol) {
      
      #cat("Iteracja", iteracja+1, "x = ", new.x, "f(x)= ", f(new.x), "\n")
      #cat("Macierz V= ", v, "\n")
      #cat("Hesjan H^-1= ",inv(hessian(f,new.x)), "\n" )
      return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
    }
    iteracja <- iteracja + 1
  }
}

# Wczytanie losowych punktów do ramki danych
library(readxl)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), iteracja = numeric(), grad_count = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        bfgs_result <- bfgs(f2, c(x_value,y_value), 1e-9)
        bfgs_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)


library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")



#### Najszybszego spadku Cauchy stałokrokowa

cauchy.const <- function(f, x, tol) {
  
  result <- tryCatch(
    expr = {
      iteracja <- 1
      grad_count <- 1
      #cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
      repeat {
        step <- 0.15
        new.x <- x - step * grad(f, x)
        iteracja = iteracja + 1
        grad_count = grad_count + 1
        
        #cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
        if (dist(rbind(new.x,x)) < tol) {
          cat(iteracja, grad_count)
          return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
        }
        x <- new.x
      }
    },
    error = function(e) {
      # Error handling code
      cat("An error occurred: ", conditionMessage(e), "\n")
      return(c(NA,NA,NA,iteracja, grad_count))  # Value to be returned in case of an error
    }
  )
  
  return(result)
  
  
}



# Wczytanie losowych punktów do ramki danych
library(readxl)

library(R.utils)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), iteracja = numeric(), grad_count = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        cauchy.const_result <- cauchy.const(f2, c(x_value, y_value), 1e-9)
        cauchy.const_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")


#### Najszybszego spadku Cauchy zmiennokrokowa

cauchy.golden <- function(f, x, tol) {
  
  result <- tryCatch(
    expr = {
      iteracja <- 1
      grad_count <- 1
      #cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
      repeat {
        g <- function(a) {f(x - a * grad(f,x))}
        step <- golden(g,0,0.1,tol)
        new.x <- x - step * grad(f, x)
        grad_count = grad_count + 2
        iteracja = iteracja + 1
        #cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
        
        if (dist(rbind(new.x,x)) < tol) {
          return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
        }
        x <- new.x
      }
    },
    error = function(e) {
      # Error handling code
      cat("An error occurred: ", conditionMessage(e), "\n")
      return(c(NA,NA,NA,iteracja, grad_count))  # Value to be returned in case of an error
    }
  )
  
  return(result)
  
  
}




# Wczytanie losowych punktów do ramki danych
library(readxl)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), f = numeric(), iteracja = numeric(), grad_count = numeric() )


for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        cauchy.golden_result <- cauchy.golden(f2, c(x_value, y_value), 1e-9)
        cauchy.golden_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}





colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")



#### Metoda gradientów sprzęzonych stałokrokowa

hestenes.stiefel <- function(f, x, tol) {
  
  result <- tryCatch(
    expr = {
      beta <- 1
      d <- -grad(f,x)
      iteracja <- 1
      grad_count <- 1
      #cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x),"\n")
      repeat {
        #g <- function(a) {f(x + a * d)}
        step <- 0.15
        new.x <- x + step * d
        iteracja = iteracja + 1
        #cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
        if (dist(rbind(new.x,x)) < tol) {
          return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
        }
        beta <- t(grad(f,new.x)) %*% (grad(f,new.x) - grad(f,x)) / (t(d) %*% (grad(f,new.x) - grad(f,x)))
        d <- -grad(f,new.x) + as.vector(beta) * d
        grad_count = grad_count + 5
        x <- new.x
      }
    },
    error = function(e) {
      # Error handling code
      cat("An error occurred: ", conditionMessage(e), "\n")
      return(c(NA,NA,NA,iteracja, grad_count))  # Value to be returned in case of an error
    }
  )
  
  return(result)
  
  
}




library(readxl)

library(R.utils)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), iteracja = numeric(), grad_count = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  result <- withTimeout({
    tryCatch(
      expr = {
        hestenes.stiefel_result <- hestenes.stiefel(f2, c(x_value, y_value), 1e-9)
        hestenes.stiefel_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")


#### Metoda gradientów sprzęzonych zmiennokrokowa

hestenes.stiefel <- function(f, x, tol) {
  
  result <- tryCatch(
    expr = {
      beta <- 1
      d <- -grad(f,x)
      iteracja <- 1
      grad_count <- 1
      #cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x),"\n")
      repeat {
        g <- function(a) {f(x + a * d)}
        step <- golden(g,0,0.1,tol)
        new.x <- x + step * d
        iteracja = iteracja + 1
        #cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
        if (dist(rbind(new.x,x)) < tol) {
          return( c(new.x[1], new.x[2], f(new.x), iteracja, grad_count))
        }
        beta <- t(grad(f,new.x)) %*% (grad(f,new.x) - grad(f,x)) / (t(d) %*% (grad(f,new.x) - grad(f,x)))
        d <- -grad(f,new.x) + as.vector(beta) * d
        grad_count = grad_count + 5
        x <- new.x
      }
    },
    error = function(e) {
      # Error handling code
      cat("An error occurred: ", conditionMessage(e), "\n")
      return(c(NA,NA,NA,iteracja, grad_count))  # Value to be returned in case of an error
    }
  )
  
  return(result)
  
  
}

library(readxl)

library(R.utils)
points = read_excel("pointsp2.xlsx", sheet=1,range="A1:B101", col_types = "numeric")

results <- data.frame(x = numeric(), y = numeric(), iteracja = numeric(), grad_count = numeric() )

for (i in 1:nrow(points)) {
  x_value <- points$x[i]
  y_value <- points$y[i]
  
  cat(x_value)
  cat(y_value)
  
  result <- withTimeout({
    tryCatch(
      expr = {
        hestenes.stiefel_result <- hestenes.stiefel(f2, c(x_value, y_value), 1e-9)
        hestenes.stiefel_result
      },
      error = function(e) {
        cat("Error: ", conditionMessage(e), "\n")
        return(rep(NA, 5))
      },
      finally = {
        cat("Execution completed.\n")
      }
    )
  }, timeout = 10, onTimeout = "silent")
  
  if (inherits(result, "timeout")) {
    result <- rep(NA, 5)
  }
  results <- rbind(results, result)
}

colnames(results) <- c("x1*", "x2*", "y*", "Liczba wywołań funkcji celu", "Liczba wywołań gradientu")
print(results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wynikip2.xlsx")




#### Wykresy 


### Metoda najszybszego spadku		# golden(g,0,0.1,tol)			

cauchy <- function(f, x, tol, results) {
  iteracja <- 1
  cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
  results <- rbind(results, c(x[1], x[2]))
  repeat {
    g <- function(a) {f(x - a * grad(f,x))}
    step <- 0.15
    new.x <- x - step * grad(f, x)
    iteracja = iteracja + 1
    cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
    results <- rbind(results, c(new.x[1], new.x[2]))
    if (dist(rbind(new.x,x)) < tol) {
      return(results)
    }
    x <- new.x
  }
}

results <- data.frame(x = numeric(), y = numeric())

results <- cauchy(f2,c(-0.53674843,-0.121136928),1e-9, results)



library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wykresy.xlsx")

### Metoda gradientów sprzężonych					

hestenes.stiefel <- function(f, x, tol, results) {
  beta <- 1
  d <- -grad(f,x)
  iteracja <- 1
  cat("Iteracja", iteracja, "x = ", x, "f(x)= ", f(x),"\n")
  results <- rbind(results, c(x[1], x[2]))
  repeat {
    g <- function(a) {f(x + a * d)}
    step <- golden(g,0,0.1,tol)	
    new.x <- x + step * d
    iteracja = iteracja + 1
    cat("Iteracja", iteracja, "x = ", new.x, "f(x)= ", f(new.x), "Krok = ", step, "\n")
    results <- rbind(results, c(new.x[1], new.x[2]))
    if (dist(rbind(new.x,x)) < tol) {
      return(results)
    }
    beta <- t(grad(f,new.x)) %*% (grad(f,new.x) - grad(f,x)) / (t(d) %*% (grad(f,new.x) - grad(f,x)))
    d <- -grad(f,new.x) + as.vector(beta) * d
    x <- new.x
  }
}
				


results <- data.frame(x = numeric(), y = numeric())

results <- hestenes.stiefel(f2,c(-0.53674843,-0.121136928),1e-9, results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wykresy.xlsx")


### Metoda quasi-Newtona	


bfgs <- function(f, x, tol, results) {
  
  n <- length(x)
  v <- diag(n)
  g <- function(a) {f(x-a*grad(f,x))}
  
  step <- 0.15
  new.x <- x - step * grad(f, x)
  
  iteracja <- 1
  repeat {
    
    cat(" Iteracja", iteracja, "x = ", x, "f(x)= ", f(x), "\n")
    results <- rbind(results, c(x[1], x[2]))
    s <- grad(f,new.x) - grad(f,x)
    r <- new.x - x
    ma <- t(r) %*% s
    a <- (r %*% t(r))/ma[1,1]
    mc <- ma
    c <- (t(s) %*% v %*% s)/mc[1,1]
    md <- ma
    d <- -(r %*% t(s) %*% v + v %*% s %*% t(r))/md[1,1]
    v <- v + (1+as.vector(c)) * a + d
    x <- new.x
    new.x <- x - v %*% grad(f,x)
    new.x <- new.x[,1]
    dist <- dist(rbind(new.x,x))
    
    if (dist(rbind(new.x, x)) < tol) {
      
      cat("Iteracja", iteracja+1, "x = ", new.x, "f(x)= ", f(new.x), "\n")
      results <- rbind(results, c(new.x[1], new.x[2]))
      cat("Macierz V= ", v, "\n")
      cat("Hesjan H^-1= ",inv(hessian(f,new.x)), "\n" )
      return(results)
    }
    iteracja <- iteracja + 1
  }
}


results <- data.frame(x = numeric(), y = numeric())

results <- bfgs(f2,c(-0.53674843,-0.121136928),1e-9, results)

library(writexl)
write_xlsx(results, "D:/materialystudia/SemestrV/Optymalizacja nieliniowa/Projekt/wykresy.xlsx")
