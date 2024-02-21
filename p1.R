f1 <- function(x){sin(x/10)*exp(-((x/10)+pi)^2) - cos(x/10)*exp(-((x/10)-2*pi)^2) + 0.003*(x/10)^2}
plot(f1, -90,100)



expandsion_bds_2pkt <- function(f, x0, delta) {
  a <- b <- x0  # Inicjalizacja przedziału poszukiwań
  i <- 0
  j <- 0
  
  
  
  while (TRUE) {
    # Poszukiwanie górnej granicy
    x_right <- x0 + 2^i * delta
    if (f(x_right) > f(x0)) {
      b <- x_right
      break
    }
    i <- i + 1
  }
    
  while (TRUE) {
    # Poszukiwanie dolnej granicy
    x_left <- x0 - 2^j * delta
    if (f(x_left) > f(x0)) {
      a <- x_left
      break
    }
    j <- j + 1
  }
    
  
  return(cat(a,b,i,j))
}

expandsion_bds_2pkt(f1, -30, 0.00000001)

interpolacja_Lagrange <- function(f, a, c, b, epsilon, gamma) {
  i <- 0
  d <- (1/2) * (f(a) * (c^2 - b^2) + f(c) * (b^2 - a^2) + f(b) * (a^2 - c^2)) /
    (f(a) * (c - b) + f(c) * (b - a) + f(b) * (a - c))
  cat(d)
  if (d<a || d>b || d==c){
    return(cat(" Algorytm nie jest zbieżny"))
  }
  while ((b - a) > epsilon) {
    i <- i + 1
    
    if (is.na(d) || is.na(c) || abs(d - c) < gamma) {
      break
    }
    
    if (a < d && d < c) {
      if (f(d) < f(c)) {
        a <- a
        c <- d
        b <- c
      } else {
        a <- d
        c <- c
        b <- b
      }
    } else if (c < d && d < b) {
      if (f(d) < f(c)) {
        a <- c
        c <- d
        b <- b
      } else {
        a <- a
        c <- c
        b <- d
      }
    }
    
    d <- (1/2) * (f(a) * (c^2 - b^2) + f(c) * (b^2 - a^2) + f(b) * (a^2 - c^2)) /
      (f(a) * (c - b) + f(c) * (b - a) + f(b) * (a - c))
  }
  
  cat("Znalezione minimum:", d, "\n")
  return(d)
}
res <- interpolacja_Lagrange(f1, -30,-19,-8,0.001, 0.00001)
res

newton_armijo <-  function(f, df, d2f, x0, alfa, p, tol) {
  iter = 1
  cat("Iteracja: ", iter, "x = ", x0, "\n")
  
  repeat {
    if (d2f(x0) > 0) {
      # Wyznaczanie długości kroku zgodnie z algorytmem Newtona
      d = -df(x0) / d2f(x0)
    } else {
      # Jeśli druga pochodna nie jest dodatnia, przyjmujemy długość kroku d = 1
      d = 1
    }
    
    # Sprawdzenie warunku Armijo
    if (f(x0 + d) <= f(x0) + alfa * d * df(x0)) {
      x_new = x0 + d
      iter = iter + 1
      cat("Iteracja: ", iter, "x = ", x_new, "\n")
      
      # Sprawdzenie warunku stopu
      if (abs(x_new - x0) < tol) {
        return(x_new)
      }
      
      x0 = x_new
    } else {
      # Zmniejszenie długości kroku
      d = p * d
    }
  }
}


