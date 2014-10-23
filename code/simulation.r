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

parameters <- list(  
  m0 = 0.05, #intrinsic mortality
  d = 0.1,  	# degradation
  c_ = 0.2, 		# beta*g  
  del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
  r = 0.01, 	# regeneration rate
  f = 0.9,		# local fascilitation
  b = 0.8, #seq(0,1,.002)+0.001
  g = 0.4, #
  assoc = TRUE,
  seed = 589232
)


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

mapping(width, height)



parms_temp <- parameters
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
while(stability > 0.000001 & i <= t_max/delta & parms_temp$rho_plus > 0) {
  
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
  
  if(parms_temp$assoc == TRUE) { # for differentiated grazing, i.e. associative resistance,
    
    # set prob of death for each cell
    death <- with(parms_temp, (m0+g*(1-Q_plus))*delta)
    
  } else { # for undifferentiated grazing
    # set prob of death for each cell
    death <- with(parms_temp, (m0+g*(1-rho_plus))*delta)				
  }
  
  # correct for overshooting death prob
  death[death > 1] <- 1
  
  regeneration <- with(parms_temp, (r + f*Q_plus)*delta)
  
  
  #} else {
  # if vegetation cover (rho) is 0
  #  regeneration <-  with(parms_temp, r*delta)
  #}
  
  degradation <- with(parms_temp, d *delta)
  
  
  
  # check for sum of probabilities to be inferior 1 and superior 0
  if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
  if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
  
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
  
  
  if(i > t_min/delta+1) {
    t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1
    t_2 <- (i-t_eval/delta):(i)
    

    if(parms_temp$rho_plus > 0) {
      stability <- (abs(mean(result$rho_plus[t_1]) - mean(result$rho_plus[t_2])))/(mean(result$rho_plus[t_1]))
    } else {
      stability <- 0
    }
    
    result$time[i] <- i*delta
    #mu_1 <- mean(result$rho_plus[t_1])
    #mu_2 <- mean(result$rho_plus[t_2])
    
    
  }
    
} # end of simulation run (over i)


par(mfrow = c(2,1))
plot(result$rho_plus[0:i] ~ result$time[0:i], 
     ylim = c(0,1), ylab = expression( rho["+"]), 
     xlab = "time [yr]", type = "l")
     

plot(result$q_plus_plus[0:i] ~ result$time[0:i], 
     ylim = c(0,1), ylab = expression(bar(q["+|+"] )), 
     xlab = "time [yr]", type = "l")
