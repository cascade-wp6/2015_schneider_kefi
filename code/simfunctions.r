########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
########################################################


######################
## mapping function ##
######################

# required to vectorise the counting and plotting. 
# returns a map of the landscape to translate it into a vector with boundaries and another one to back-translate it to a vector without boundaries into the global environment. Needs to be called only once for the dimensions of the lattice. 


mapping <- function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
    
  # derive helper vectors for counting: 
  # transformation vector for evaluation at the border of the grid
  # set evaluation matrix 
  X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  # setting the border of the evaluation matrix X
  X <- cbind(X[,width], X, X[,1] )  
  X <- rbind(X[height,], X, X[1,] ) 
  # transformation vector which adds the border to the lattice:
  x_with_border <- as.integer(t(X))
  
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  # from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
  #x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )  	
  
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  # defining the neighborhood which is to be evaluated	
  # set interaction matrix
  I <- i_matrix	
  # coordinates of neighbours in Interaction matrix I: 
  neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
  # coordinates relative to the evaluated cell (=  which(is.na(I) ) 
  relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  # relative position of the four direct neighbours of a cell
  #interact <- (relrow * dim(X)[2] + relcol)
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}



####################
## count function ##
####################

count  <- function(x, neighbor) {
  
  neighbors <- numeric(length = prod(x$dim))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}


#################
## patch count ##
#################

# identify and count the size of individual patches

patches <- function(x, state) {
  pattern <- x$cells
  pattern <- pattern %in% state
  map <- rep(NA, times = prod(x$dim))
  old <- rep(99, times = prod(x$dim)) 
  
  while(!identical(old[pattern], map[pattern])) {
    old <- map
    count = as.integer(1)
    for(i in which(pattern)) {
      neighbors <- map[x_with_border][x_to_evaluate[i]+interact]
      if(all(is.na(neighbors)) ) { 
        map[i] <- count
      } else {
        map[i] <- min(neighbors, na.rm = TRUE)
      }
      count <- count +1
    }
    
  }
  
  map <- as.factor(map)
  patchvec <- as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
  
  out <- vector()
  if(length(patchvec) > 0) out <- sort(patchvec) else out <- NA
  #out <- data.frame(size = unique(out), freq = sapply(unique(out), function(i) length(which(out >= i)) ))
  return(out)
  
} 



fitPL <- function(psd, p_spanning, n = NULL) {
  
  # code of fitted classes
  
  n_plants <- sum(psd$size * psd$n)/n
  
  out <- list()
  out$best <- NA
  out$AIC <- vector("numeric", length = 3)
  out$dAIC <- vector("numeric", length = 3)
  
  # criteria for vegetated state & desert state
  
  ##### linear power law model for parameter estimation
  PLlm <- lm(I(log(p)) ~  1 - I(log(size)) , data = psd) 
  
  ###########
  
  try( {out$TPLdown <- nls(I(log(p)) ~ alpha * log(size) + Sx * (1 - size) , 
                           data = psd,
                           start = list(alpha =  PLlm$coefficients, Sx = 1/1000),
                           #algorithm = "port",
                           trace = FALSE
  )}, silent = TRUE
  )    
  
  if(!is.null(out$TPLdown) & !coefficients(out$TPLdown)["Sx"] <= 0) {
    out$AIC[1] <- AIC(out$TPLdown) 
  } else {
    out$TPLdown <- list(NA)
    out$AIC[1] <- NA
  }
  
  #####
  
  try({out$PL <- nls(I(log(p)) ~ alpha * log(size), 
                     data = psd,
                     start = list( alpha =  PLlm$coefficients ),
                     trace = FALSE,
                     nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  if(!is.null(out$PL)) {
    out$AIC[2] <- AIC(out$PL)
  } else {
    out$PL  <- list(NA)
    out$AIC[2] <- NA
  }
  
  ###########
  
  try({out$TPLup <- nls(I(log(p)) ~  log(b) + log(1+(size^(alpha))/b ) , 
                        data = psd,
                        start = list( alpha =  PLlm$coefficients, b = p_spanning ) , 
                        nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  
  if(!is.null(out$TPLup)) {
    out$AIC[3] <- AIC(out$TPLup) 
  } else { 
    #result$fit$summary$TPLup  <- list(NA)
    out$TPLup  <- list(NA)
    out$AIC[3] <- NA
  }
  
  ###########
  
  out$dAIC <-   out$AIC -min(out$AIC, na.rm = TRUE)
  
  out$best <- which.min(out$AIC)+1
  
  return(out)
} 


