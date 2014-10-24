########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider & Sonia Kéfi
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

rm(list=ls())

########################################################################################
source("code/simfunctions.r")

################ parameter settings
first_ID = 1  #
set.seed(855235)

global <- list(  
  m0 = 0.05, #intrinsic mortality
  d = 0.1,  	# degradation
  c_ = 0.2, 		# beta*g  
  del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
  r = 0.01, 	# regeneration rate
  f = 0.9		# local fascilitation
)

parameters <- list(
  b = seq(0,1,.02), #seq(0,1,.002)+0.001
  g = c(0.1,0.4), #
  protect = c(0,0.5),
  starting = NA, # by default starting from low cover
  seed = NA
)

iterations <- expand.grid(parameters)

iterations$seed <- sample(10^6:10^7, length(iterations$b)  )
iterations$starting <- runif(length(iterations$b) , 0.8,0.9)
iterations <- cbind(ID = 1:dim(iterations)[1],iterations, global)
str(iterations)

# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
color <- c("black","grey80", "white") # define colors for the cell state levels


# time and resolution of simulation
#timesteps = 1000
delta = 1/5
t_min = 500
t_max = 2500
t_eval <- 200
each = 100

n_snaps <- length(seq(0,2*t_eval,each))
t_snaps <-  ceiling(t_max/(each*n_snaps))*(n_snaps*each)
snapshots <- data.frame(time = seq(100, t_snaps, each), i= seq(100, t_snaps, each)/delta+1, pos = c(1:n_snaps) )

################ map objects
mapping(100,100)


################ output saving

sim <- "21_static"



################ starting parallel backend

# see documentation !!!

library(foreach)
library(doSNOW)

cl <- makeSOCKcluster(rep("localhost", times = 2), outfile='out_messages.txt')
registerDoSNOW(cl)

#################################
## desert state simulation
#################################

if(TRUE){
replicates = 100


################ starting foreach loop
foreach(iteration = iterations$ID, .combine = rbind) %dopar% {
  
  #foreach(iterations, .combine = "c") %do% {
  #subset(iterations, g == 0.3 & b == 0.8 & starting == 0.0001)
  #iteration =    4172   #2051
  
  set.seed(iterations$seed[iteration])
  
  parms <- as.list(iterations[iterations$ID == iteration,])
  
  #foreach(repl = 1:replicates, .combine = rbind) %do% {
  temp <- data.frame(
    replicate = integer(replicates),
    rho_plus = numeric(replicates), 
    q_plus_plus = numeric(replicates),
    recover = logical(replicates),
    runtime = numeric(replicates)
  )
  
  for(repl in 1:replicates) {
    
    parms_temp <- parms #as.list(iterations[iterations$ID == iteration,])
    
    # sampling the initial random grid into a list object
    parms_temp$rho_plus <- 0.001
    
    # how many initial plants
    init_plant <- as.integer(width*height*parms_temp$rho_plus)
    
    # vector of empty cells in state "0"
    cells <- factor(rep("-", times = width*height), levels = states) # 1:length(states))
    # replace init_plant cells, randomly drawn, with "+"
    cells[sample(1:(width*height), init_plant, replace = FALSE)] <- "+"
    
    # which cells are still empty?
    empty <- which(cells != "+")
    # select which will be degraded? fixed to 50% of non occupied cells. 
    
    initial <- list(  
      dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
      cells = cells#contains a random row-wise, factorial vector to fill the grid 
    )
    levels(initial$cells) <- states  #assign cell states 
    class(initial) <- c("list","landscape") # set class of object (required for plotting)
    
    parms_temp$Q_plus <- count(initial, "+")/4
    rnum <- runif(width*height) # one random number between 0 and 1 for each cell
    
    regeneration <- with(parms_temp, (r + f*Q_plus))
    initial$cells[which(initial$cells == "-"   & rnum <= regeneration)] <- "0"
    
    #plot(initial)
    
    ##### set result obj
    
    result <- list()
    
    result$time <-as.integer(0)#seq(0, delta, delta) # write simulated timesteps
    
    #result$mortality <- numeric(length = timesteps/delta+1)
    #result$mortality_border <- numeric(length = timesteps/delta+1)
    
    #save all level global densities: rho_*
    result$rho_plus  <- vector("numeric", length = length(result$time))
    
    result$rho_plus[1] <- parms_temp$rho_plus
    
    
    result$q_plus_plus  <- vector("numeric", length = length(result$time))  		
    result$q_plus_plus[1]  <- mean(subset(count(initial, "+")/4, initial$cells == "+") ) 		
    
    # result$snapshots <- integer(length = n_snaps)
    
    # result$timeseries <- list()
    # result$timeseries[1:n_snaps] <- rep(list(initial), times = n_snaps) 
    
    # ## the following are only stored for the last five snapshots (overridden until equ. is reached)
    # result$patches <- list()
    # result$cumpatch <- list()	
    # result$models <- list()
    # result$snapstats <- data.frame(time = vector("numeric", length = n_snaps), max_patch = vector("numeric", length = n_snaps), rho_plus = vector("numeric", length = n_snaps) )
    # #	result$snapshots <- data.frame(time = vector("numeric", length = 5), max_patch = vector("numeric", length = 5), bestmodel = factor(rep(NA,length=5), levels = c("DES","PL", "TPLup", "TPL", "FULL")), alpha = vector("numeric", length = 5))
    
    x_old <- initial
    
    #stability <- 1
    
    i = 1
    #for(i in seq_along(result$time)[-1]) {    #calculation loop #2:(timesteps/delta+1)
    while(result$rho_plus[i] > 0 & result$rho_plus[i] <= 0.01 & i <= t_max/delta) {
      i <- i +1
      x_new <- x_old 		# copy x_old into an object x_new to allocate memory	
      
      # model specific part:
      # 1 - setting time-step parameters
      
      if(parms_temp$rho_plus == 0) {flag <- FALSE} else {flag <- TRUE}
      
      # 2 - drawing random numbers
      rnum <- runif(width*height) # one random number between 0 and 1 for each cell
      
      # 4 - applying the rules to fill the cells in x_new
      
      if(flag) { # if density is unequal 0, then
        
        # count local density of occupied fields for each cell: 
        parms_temp$Q_plus <- count(x_old, "+")/4
        
        # calculate recolonisation rates of all cells
        recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta)
        
        # calculate death rates
        
        #if(parms_temp$assoc == TRUE) { # for differentiated grazing, i.e. associative protection,
        
        # set prob of death for each cell
        death <- with(parms_temp, (m0+g*(1-protect))*delta)
        #}
        # correct for overshooting death prob
        death[death > 1] <- 1
        
        regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
        
      } else {
        # if vegetation cover (rho) is 0
        recolonisation <-  with(parms_temp, r*delta)
        death <- 0
        Q_plus <- 0
      }
      
      degradation <- with(parms_temp, (d *delta))
      
      
      
      # check for sum of probabilities to be inferior 1 and superior 0
      if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
      if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
      
      # apply rules 
      
      
      if(flag) {	
        x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
        x_new$cells[which(x_old$cells == "+"  & rnum <= death)] <- "0"
      }
      
      
      x_new$cells[which(x_old$cells == "0"  & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
      x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"
      
      
      # 5 saving state of the new grid		
      
      parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) 
      result$rho_plus[i] <- parms_temp$rho_plus
      
      result$q_plus_plus[i]  <- mean(subset(count(x_new, "+")/4, x_new$cells == "+") ) 		
      
      x_old <- x_new
      result$time[i] <- i*delta 
      
      
    } # end of simulation.
    
    ####### processing saved snapshots of the lattice
    
    
    t_fin = result$time[length(result$time)]
    
    result$out <- data.frame(
      replicate = repl,
      rho_plus = result$rho_plus[i], 
      q_plus_plus = result$q_plus_plus[i],
      recover = result$rho_plus[i] > 0.01,
      runtime = t_fin
    )
    
    temp[repl,] <- result$out
    
    #return(result$out)
  } #-> temp
  
  
  out <- data.frame(
    ID = iteration,
    starting = 0.001, 
    #assoc = parms$assoc,
    g = parms$g, 
    b = parms$b, 
    m0 = parms$m0,
    rho_plus = mean(temp$rho_plus[temp$recover]), 
    rho_plus_sd = sd(temp$rho_plus[temp$recover]), 
    q_plus_plus = mean(temp$q_plus_plus[temp$recover]),
    q_plus_plus_sd = sd(temp$q_plus_plus[temp$recover]),
    growth = mean(temp$rho_plus/temp$runtime),
    n_recover = sum(temp$recover),
    prob_recover = mean(temp$recover), 
    runtime = mean(temp$runtime)
  )
  
  #save(result, file = paste(path, "result_", iteration, sep = ""))
  
  gc() 
  return(out)
  
} -> output


write.csv(output, file = paste("data/sim", sim,"/", "output_deg.csv", sep = "") )
}


#################################
## vegetated state simulation
#################################
replicates = 10

min_replicates = 10
max_tries = 50


# run parameters (in parallel)

foreach(iteration = iterations$ID, .combine = rbind) %dopar% { 
  #subset(iterations, g == 0.10 & b == 0.50)
  
  if(Sys.info()[['nodename']] == "kefi118") {
    setwd("/home/schneider/herbivory/")
    #dir.create("data")
    #dir.create(paste("data/sim", sim,"/", sep = ""))
    #dir.create(paste("data/sim", sim,"/results/", sep = ""))
  } #for Linux Workstation 
  
  #iteration = 63 #9 
  
  if(iteration %in% done) {
    
    load(paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))  
    
  } else {
    
    set.seed(iterations$seed[iteration])
    
    # collecting final lattices for each run in j
    collect <- list()
    collect$lattice <- list()
    collect$out <- data.frame(
      replicate = NA,
      rho_plus_ini =NA,
      
      rho_plus = NA, 
      rho_plus_sd = NA,
      rho_plus_fin = NA,
      
      q_plus_plus = NA,
      
      clus_coef = NA,
      clus_coef_sd = NA,
      
      stability = NA,
      runtime = NA
    )[-1,]
    
    
    # run replicates
    
    j = 1
    n = 1
    
    while(n < min_replicates & j < max_tries) {
      
      #for(j in 1:10) {
      
      parms_temp <- as.list(subset(iterations, ID == iteration))
      parms_temp$rho_plus <- runif(1, 0.8,0.9)
      
      ###### initialise grid
      # how many initial plants
      init_plant <- as.integer(width*height*parms_temp$rho_plus)
      
      # vector of empty cells in state "0"
      cells <- factor(rep("0", times = width*height), levels = states) # 1:length(states))
      # replace init_plant cells, randomly drawn, with "+"
      cells[sample(1:(width*height), init_plant, replace = FALSE)] <- "+"
      
      # which cells are still empty?
      empty <- which(cells != "+")
      # select which will be degraded? fixed to 50% of non occupied cells. 
      init_degraded <- sample(empty, length(empty) * (parms_temp$r + parms_temp$f/2) )  
      # replace cell state. 
      cells[init_degraded] <- "-"
      
      initial <- list(  
        dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
        cells = cells#contains a random row-wise, factorial vector to fill the grid 
      )
      levels(initial$cells) <- states  #assign cell states 
      class(initial) <- c("list","landscape") # set class of object (required for plotting)
      
      
      #### initialise result object
      result <- list()
      result$time <- seq(0, t_min, delta)
      
      result$rho_plus <- vector("numeric", length = length(result$time))
      result$rho_plus[1] <- parms_temp$rho_plus
      
      result$q_plus_plus  <- vector("numeric", length = length(result$time))			
      result$q_plus_plus[1]  <- mean(subset(count(initial, "+")/4, initial$cells == "+") ) 		
      
      
      x_old <- initial
      
      stability <- 1
      #stability <- FALSE
      i = 1
      #t_max = 50
      # starting iterations
      while(stability > 0.000001 & i <= t_max/delta & parms_temp$rho_plus > 0 ) {
        
        i <- i +1
        x_new <- x_old 		# copy x_old into an object x_new to allocate memory	
        
        # model specific part:
        # 1 - setting time-step parameters
        
        #if(parms_temp$rho_plus == 0) {flag <- FALSE} else {flag <- TRUE}
        #if(flag) {
        # count local density of occupied fields for each cell: 
        parms_temp$Q_plus <- count(x_old, "+")/4
        #  } else {
        #  Q_plus <- 0
        #}
        
        # 2 - drawing random numbers
        rnum <- runif(width*height) # one random number between 0 and 1 for each cell
        
        # 4 - applying the rules to fill the cells in x_new
        
        #if(flag) { # if density is unequal 0, then
        
        # calculate recolonisation rates of all cells
        recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta)
        
        # calculate death rates
        
          # set prob of death for each cell
        death <- with(parms_temp, (m0+g*(1-protect))*delta)
          
        # correct for overshooting death prob
        death[death > 1] <- 1
        
        regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
        
        
        #} else {
        # if vegetation cover (rho) is 0
        #  regeneration <-  with(parms_temp, r*delta)
        #}
        
        degradation <- with(parms_temp, d *delta)
        
        
        
        # check for sum of probabilities to be inferior 1 and superior 0
        if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in run", iteration, "time step", i, "! decrease delta!!!")) 
        if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in run", iteration, "in time step", i, "! balance parameters!!!")) 
        
        # apply rules 
        
        
        #if(flag) {
        x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
        x_new$cells[which(x_old$cells == "+"  & rnum <= death)] <- "0"
        #}
        
        x_new$cells[which(x_old$cells == "0"  & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
        x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"
        
        
        # 5 saving state of the new grid		
        
        parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) 
        result$rho_plus[i] <- parms_temp$rho_plus
        
        result$q_plus_plus[i]  <- mean(subset(count(x_new, "+")/4, x_new$cells == "+") ) 		
        
        x_old <- x_new
        
        
       
        if(parms_temp$rho_plus > 0) {
            
          if(i > t_min/delta+1 | stability <= 0.000001) {
            
            t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1
            t_2 <- (i-t_eval/delta):(i)
            stability <- (abs(mean(result$rho_plus[t_1]) - mean(result$rho_plus[t_2])))/(mean(result$rho_plus[t_1]))
            result$time[i] <- i*delta
            
          }
          
        } else {
          stability <- 0
          if(i > t_min/delta+1) result$time[i] <- i*delta
        }
          
          
          
       
        
        
        #  if(i %in% snapshots$i) {
        #   
        #   pos <-  snapshots$pos[match(i, snapshots$i)]
        #    result$timeseries[[pos]] <- x_new 
        #    result$snapstats$time[pos] <- snapshots$time[match(i, snapshots$i)]
        #    result$snapstats$rho_plus[pos] <- parms_temp$rho_plus
        #    
        #    
        #  }
        
       # if(i < t_min/delta+1 ){ xlim <- c(0,t_min) } else { xlim <- c( result$time[i] - t_min,  result$time[i])}
       # plot(result$rho_plus[0:i] ~ result$time[0:i], ylim = c(0,1), main = round(stability, digits =4 ), xlim = xlim, type = "l")
        # plot(x_new)
        
      } # end of simulation run (over i)
      
      (j = j + 1)
      
      if(stability < 0.000001 ) { #parms_temp$rho_plus == 0 | 
              
        collect$lattice[[n]] <- x_new
        
        t_fin = result$time[length(result$time)]
        i_stable <- (length(result$time)-(2*t_eval/delta)):length(result$time)
        
        # switch for output at extinction
        if(parms_temp$rho_plus > 0) {
          result$out <- data.frame(
            
            replicate = j,
            rho_plus_ini = result$rho_plus[1],
            
            rho_plus = mean(result$rho_plus[i_stable]), 
            rho_plus_sd = sd(result$rho_plus[i_stable]),
            rho_plus_fin = result$rho_plus[length(result$time)],
            
            q_plus_plus = mean(result$q_plus_plus[i_stable]),
            
            clus_coef = mean(result$q_plus_plus[i_stable]/result$rho_plus[i_stable]),
            clus_coef_sd = mean(result$q_plus_plus[i_stable]/result$rho_plus[i_stable]),
            
            stability = stability,
            runtime = t_fin
          )
        } else {
          result$out <- data.frame(
            
            replicate = j,
            rho_plus_ini = result$rho_plus[1],
            
            rho_plus = 0, 
            rho_plus_sd = NA,
            rho_plus_fin = result$rho_plus[length(result$time)],
            
            q_plus_plus = 0,
            
            clus_coef = NA,
            clus_coef_sd = NA,
            
            stability = 0,
            runtime = t_fin
          ) 
        }
        
        collect$out[n,] <- result$out
        
        (n = n + 1) 
      }
      
    }#-> collect$out
    
    
    result <- list()
    
    # pool replicates to one row (means, sd)
    stable <- collect$out$stability < 0.000001
    
    result$runs <- collect$out 
    result$grids <- collect$lattice
    
    
    result$out <- data.frame(
      ID = iteration,
      protect = parms_temp$protect,
      g = parms_temp$g, 
      b = parms_temp$b, 
      m0 = parms_temp$m0,
      
      n = length(which(stable)),
      runtime = mean(result$runs$runtime[stable]), 
      runtime_sd = sd(result$runs$runtime[stable]), 
      
      rho_plus_ini = mean(result$runs$rho_plus_ini[stable]),
      rho_plus = mean(result$runs$rho_plus[stable]),
      rho_plus_sd = mean(result$runs$rho_plus_sd[stable]),
      
      q_plus_plus = mean(result$runs$q_plus_plus[stable]),
      
      clus_coef = mean(result$runs$clus_coef[stable]),
      largest_patch = NA
    )
    
    
    # calculate cpsd for pooled lattices in collect
    
    result$p <- list()
    result$cumpatch <- list()
    result$runs$largest_patch <-NA
    
    
    for(j in (1:length(result$runs$replicate))[result$runs$rho_plus_fin != 0 & stable]) {
      result$p[[j]] <- patches(result$grids[[j]],"+")
      
      if( !is.na(result$p[j] )) {
        cumbins <- sort(unique(unlist(result$p[j]))) 
        #bins <- seq(1,10001, 10)
        
        result$cumpatch[[j]] <- data.frame(size = cumbins)
        result$cumpatch[[j]]$n <- sapply(cumbins, function(k) length(which(result$p[[j]] >= k)) )
        result$cumpatch[[j]]$p <- sapply(cumbins, function(k) length(which(result$p[[j]] >= k))/length(result$p[[j]]) ) 
        
      } else {
        result$cumpatch[[j]] <- NA
      }
      
      result$runs$largest_patch[j] <- max(result$p[[j]])
      
    }
    
    
    
    
    result$out$largest_patch <- mean(result$runs$largest_patch, na.rm = TRUE)
    
    
    
    result$fit <- list()
    result$fit$best <- NA
    
    

  
  
  stable <- result$runs$stability < 0.000001
  
  ## stage 1: test for desert or vegetated, else 
  result$out$best_model = NA
  
  ### check if desert
  if(result$out$rho_plus < 0.01) { result$out$best_model = 1 } 
  
  ### check if vegetated 
  if(is.na(result$out$best_model) & result$out$rho_plus >= 0.80 ) {
    #dd4 <- do.call("rbind", result$cumpatch)
    result$out$best_model = 5
    
  } 
  ## stage 2: fit power law models and compare via AIC via function
  
  if(is.na(result$out$best_model)) {
    dd4 <- do.call("rbind", result$cumpatch)
    p_spanning <- mean(unlist(sapply(1:length(result$cumpatch), function(x) tail(result$cumpatch[[x]]$p, 1)) )) 
    
    result$fit <- fitPL(dd4, p_spanning, n = result$out$n)
    
    result$out$best_model = result$fit$best 
  }
  
  
  result$out$class = c("DEG", "DOWN","PL", "UP", "COV")[result$out$best]
  
  
  
  
  t_fin = result$time[length(result$time)]
  #i_stable <- (length(result$time)-(2*t_eval/delta)):length(result$time)
  result$out$largest_patch = mean(result$runs$largest_patch, na.rm = TRUE) 
  result$out$largest_patch_sd = sd(result$runs$largest_patch, na.rm = TRUE)
  
  
  #result$out <- cbind(result$out, data.frame(alpha = NA, alpha_se = NA, alpha_p = NA, trunc = NA, trunc_se = NA, trunc_p = NA, lower = NA, lower_sd = NA, lower_p = NA)) 
  
  result$out$alpha     <- NA
  result$out$alpha_se  <- NA 
  result$out$alpha_p  <- NA
  
  result$out$trunc  <- NA
  result$out$trunc_se  <- NA
  result$out$trunc_p  <- NA
  
  result$out$lower  <- NA
  result$out$lower_se  <- NA
  result$out$lower_p  <- NA
  
  if(! result$out$class %in% c("DEG","COV")) {
    result$out$alpha <- coefficients(result$fit[[result$out$best+2]])["alpha"]
    result$out$alpha_se <-  summary(result$fit[[result$out$best+2]])$coefficients["alpha",2]
    result$out$alpha_p <-  summary(result$fit[[result$out$best+2]])$coefficients["alpha",4]
    
  }
  
  if(result$out$class %in% c("UP")) {
    result$out$lower <- coefficients(result$fit$TPLup)["b"]
    result$out$lower_se <-  summary(result$fit$TPLup)$coefficients["b",2]
    result$out$lower_p <-  summary(result$fit$TPLup)$coefficients["b",4]
  }
  
  if(result$out$class %in% c("DOWN")) {
    result$out$trunc <- coefficients(result$fit$TPLdown)["Sx"]
    result$out$trunc_se <- summary(result$fit$TPLdown)$coefficients["Sx",2]
    result$out$trunc_p <- summary(result$fit$TPLdown)$coefficients["Sx",4]
  }
  
  #result$out$updated <- TRUE
  
  save(result, file = paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))
  
  }  # end if(!iteration %in% done) 
  
  gc() 
  return(result$out)
  
}#-> output

#write.csv(output, file = paste("data/sim", sim,"/", "output.csv", sep = "") )



#collect output files to server

system(paste0("scp schneider@162.38.184.118:~/herbivory/data/sim", sim, "/results/*.* data/sim", sim, "/results/"))

foreach(iteration = iterations$ID, .combine = rbind) %do% {     
    load(paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))  
    return(result$out)
  }-> output

write.csv(output, file = paste("data/sim", sim,"/", "output.csv", sep = "") )
  



stopCluster(cl)
