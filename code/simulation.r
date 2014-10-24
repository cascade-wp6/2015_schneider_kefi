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

set.seed(589232) # setting seed for random number generation

# parameter settings:
parameters <- list(  
  m0 = 0.05, #intrinsic mortality
  d = 0.1,  	# degradation
  c_ = 0.2, 		# beta*g  
  del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
  r = 0.01, 	# regeneration rate
  f = 0.9,		# local fascilitation
  b = 0.78, # environmental quality
  g = 0.4, # grazing intensity
  assoc = TRUE # TRUE: associational resistance (spatially-explicit); FALSE: global resistance (mean field) 
)

# specify lattice:
width = 100
height = 100

# cell states:
states = c("+","0","-")  # vegetated, empty and degraded, respectively

# time and resolution of simulation:
delta = 1/5     # resolution of simulation, each year is broken into 5 updates
t_min = 500     # minimal number of years to simulate
t_max = 2500    # maximal number of years to simulate
t_eval <- 200   # timespan used to evaluate stability
each = 100      # save snapshots each 100 years

# calculate timesteps for snapshots:
n_snaps <- length(seq(0,2*t_eval,each))
t_snaps <-  ceiling(t_max/(each*n_snaps))*(n_snaps*each)
snapshots <- data.frame(time = seq(100, t_snaps, each), i= seq(100, t_snaps, each)/delta+1, pos = c(1:n_snaps) )

# mapping vectors: 
mapping(width, height)   # create mapping vectors (see documentation)

# initialize landscape: 
parms_temp <- parameters   # set temporary parameter object
parms_temp$rho_plus <- runif(1, 0.8,0.9)   # draw initial vegetation cover

init_plant <- as.integer(width*height*parms_temp$rho_plus) # how many initial plants
cells <- factor(rep("0", times = width*height), levels = states) # vector of empty cells in state "0"
cells[sample(1:(width*height), init_plant, replace = FALSE)] <- "+"  # replace init_plant cells, randomly drawn, with "+"
empty <- which(cells != "+") # which cells are still empty?
init_degraded <- sample(empty, length(empty) * (parms_temp$r + parms_temp$f/2) )   # select which will be degraded? fixed to 50% of non occupied cells
cells[init_degraded] <- "-"  # replace cell state of degraded cells

# wrap landscape object:
initial <- list(  
  dim = c(as.integer(width), as.integer(height)),  # first element contains the dimensions of the landscape 
  cells = cells  #contains a random row-wise, factorial vector to fill the grid 
)
levels(initial$cells) <- states  #assign cell states 
class(initial) <- c("list","landscape") # set class of object (required for plotting)

# initialise result object:
result <- list()  # generate list object
result$time <- seq(0, t_min, delta) # add vector of realized timesteps
result$rho_plus <- vector("numeric", length = length(result$time))  # allocate vector for global vegetation cover 
result$rho_plus[1] <- parms_temp$rho_plus   # fill in first value of vegetation cover
result$q_plus_plus  <- vector("numeric", length = length(result$time))  # allocate vector for average local cover
result$q_plus_plus[1]  <- mean(subset(count(initial, "+")/4, initial$cells == "+") )  # fill in first vector of average local cover

# --------------- simulation -----------------

# initialise simulation variables: 
x_old <- initial  # ghost matrix at t_i
stability <- 1  # check value for stability
i = 1  # iterator for simulation timesteps

# starting iterations:
while(stability > 0.000001 & i <= t_max/delta & parms_temp$rho_plus > 0) {
  
  i <- i +1  # increase iterator
  x_new <- x_old 		# copy x_old into an object x_new to allocate memory	
  
  # model specific part:
  # 1 - setting time-step parameters
  parms_temp$Q_plus <- count(x_old, "+")/4  # count local density of occupied fields for each cell:
   
  # 2 - drawing random numbers
  rnum <- runif(width*height) # one random number between 0 and 1 for each cell
  
  # 3 - setting transition probabilities
  recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*delta) # recolonisation rates of all cells 
  if(parms_temp$assoc == TRUE) { # for differentiated grazing, i.e. associative resistance,
    death <- with(parms_temp, (m0+g*(1-Q_plus))*delta)   # set probability of death for each cell
  } else { # for undifferentiated grazing
    death <- with(parms_temp, (m0+g*(1-rho_plus))*delta) # set probability of death for each cell	
  }
  death[death > 1] <- 1   # correct for overshooting death probability
  regeneration <- with(parms_temp, (r + f*Q_plus)*delta) # set probability of regeneration for each cell 
  degradation <- with(parms_temp, d *delta) # set probability of degradation for each cell
  
  # check for sum of probabilities to be inferior 1 and superior 0
  if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in time step", i, "! decrease delta!!!")) 
  if(any(c(recolonisation, degradation, death, regeneration) < 0)) warning(paste("a set probability falls below 0 in time step", i, "! balance parameters!!!")) 
  
  # 4 - apply transition probabilities  

  x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
  x_new$cells[which(x_old$cells == "+" & rnum <= death)] <- "0"
  x_new$cells[which(x_old$cells == "0" & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
  x_new$cells[which(x_old$cells == "-" & rnum <= regeneration)] <- "0"
  
  # 5 - saving state of the new grid		
  
  parms_temp$rho_plus <- sum(x_new$cells == "+")/(width*height) # store current global vegetation cover
  result$rho_plus[i] <- parms_temp$rho_plus # write current vegetation cover to result
  result$q_plus_plus[i]  <- mean(subset(count(x_new, "+")/4, x_new$cells == "+") ) # write current average local cover to result
  x_old <- x_new # replace ghost matrix for next iteration
  
  if(i > t_min/delta+1) { # if we are over the minimal timespan 
    t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1 # vector of t_eval timesteps previous to the last t_eval timesteps
    t_2 <- (i-t_eval/delta):(i) # vector of the last t_eval timesteps 
    
    if(parms_temp$rho_plus > 0) { 
      stability <- (abs(mean(result$rho_plus[t_1]) - mean(result$rho_plus[t_2])))/(mean(result$rho_plus[t_1])) # calculate stability, i.e. difference between the mean cover in the two evaluation periods t_1 and t_2 
    } else {
      stability <- 0 # set stability to 0 if cover is 0, immediate stop of simulation
    }
    
    result$time[i] <- i*delta # save timestep to results
    
  }
    
} # end of simulation run (over i)


par(mfrow = c(2,1))
plot(result$rho_plus[0:i] ~ result$time[0:i], 
     ylim = c(0,1), ylab = expression( rho["+"]), 
     xlab = "time [yr]", type = "l")
     

plot(result$q_plus_plus[0:i] ~ result$time[0:i], 
     ylim = c(0,1), ylab = expression(bar(q["+|+"] )), 
     xlab = "time [yr]", type = "l")
