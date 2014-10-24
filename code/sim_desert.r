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
set.seed(213823)

global <- list(	
	m0 = 0.05, #intrinsic mortality
	d = 0.1,		# degradation
	c_ = 0.2, 		# beta*g  
	del = 0.1,  	# seeds dispersed; (1-del) seeds on nearest neighbourhood
	r = 0.01, 	# regeneration rate
	f = 0.9		# local fascilitation
	)
	
parameters <- list(
	b = seq(0,1,.02), #seq(0,1,.002)+0.001
	g = seq(0,.5,.01), #c(0,0.025, 0.05, 0.075, 1), #.025 #
	assoc = c(FALSE, TRUE),
	starting = c(0.001) # c(0.1, 0.2, 0.5, .99),
)

iterations <- expand.grid(parameters)

iterations$seed <- sample(10^6:10^7, length(iterations$b)  )

iterations <- cbind(ID = 1:dim(iterations)[1],iterations, global)
str(iterations)



replicates = 100

# specify lattice
width = 100
height = 100

# initial cell states
states = c("+","0","-")
#prob = c(9/10,.9/10,0.1/10)
color <- c("black","grey80", "white") # define colors for the cell state levels

# time and resolution of simulation
#timesteps = 1000
delta = 1/5
t_min = 1
t_max = 100
# t_eval <- 500
	# n_snaps <- length(seq(0,2*t_eval,100))
	# t_snaps <-  ceiling(t_max/(100*n_snaps))*(n_snaps*100)
	# snapshots <- data.frame(time = seq(100, t_snaps, 100), i= seq(100, t_snaps, 100)/delta+1, pos = c(1:n_snaps) )

################ helper objects

# derive helper vectors for counting: 
# transformation vector for evaluation at the border of the grid
	# set evaluation matrix 
	X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
	# setting the border of the evaluation matrix X
	X <- cbind(X[,width], X, X[,1] )  
	X <- rbind(X[height,], X, X[1,] ) 
	# transformation vector which adds the border to the lattice:
	x_with_border <- as.integer(t(X))
	
	# from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
	x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]	)		

# defining the neighborhood which is to be evaluated	
	# set interaction matrix
	I <- matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)	
	# coordinates of neighbours in Interaction matrix I: 
	neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
	# coordinates relative to the evaluated cell (=  which(is.na(I) ) 
	relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
	relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
	
	# relative position of the four direct neighbours of a cell
	interact <- (relrow * dim(X)[2] + relcol)

	
	
################ output saving

sim <- "18_deg"
dir.create( paste("data/sim", sim,"/", sep = ""))
path <- paste("data/sim", sim,"/results/", sep = "")
dir.create(path)

################ starting parallel backend

# see documentation !!!

library(foreach)
library(doSNOW)

cl <- makeSOCKcluster(rep("localhost", times = 2), outfile='out_messages.txt')
registerDoSNOW(cl)

################ preallocate memory
temp_num <- numeric(length(iterations$ID))
temp_char <- character(length(iterations$ID))
temp_logi <- logical(length(iterations$ID))
temp_int <- integer(length(iterations$ID))

output <- data.frame(
			ID = temp_int,
			starting = temp_num, 
			assoc = temp_logi,
			g = temp_num, 
			b = temp_num, 
			m0 = temp_num,
			rho_plus = temp_num, 
			rho_plus_sd = temp_num,
			q_plus_plus = temp_num,
			q_plus_plus_sd = temp_num,
			speed = temp_num,
			runtime = temp_num
	)
	
	

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
	parms_temp$rho_plus <- parms_temp$starting
	
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
	
		if(parms_temp$assoc == TRUE) { # for differentiated grazing, i.e. associative protection,
				
			# set prob of death for each cell
			death <- with(parms_temp, (m0+g*(1-Q_plus))*delta)
			
		} else { # for undifferentiated grazing
			# set prob of death for each cell
			death <- with(parms_temp, (m0+g*(1-rho_plus))*delta)				
		}
		
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
		
		# if(i > t_min/delta+1) {
		# t_1 <- (i-2*t_eval/delta):(i-t_eval/delta)-1
		# t_2 <- (i-t_eval/delta):(i)
		
		# stability <- abs(mean(result$rho_plus[t_1]) - mean(result$rho_plus[t_2]))
		# result$time[i] <- i*delta
		# }
		
		# if(i %in% snapshots$i) {

			# pos <-  snapshots$pos[match(i, snapshots$i)]
			# result$timeseries[[pos]] <- x_new 
			# result$snapstats$time[pos] <- snapshots$time[match(i, snapshots$i)]
			# result$snapstats$rho_plus[pos] <- parms_temp$rho_plus

		# }
		
		### check output, replace snapshot data.frame ??? 
		
		#if(i %in% seq(0,t_max,10)/delta) {
		#if(i < t_min/delta+1 ){ xlim <- c(0,t_min) } else { xlim <- c( result$time[i] - t_min,  result$time[i])}
		#plot(result$rho_plus[0:i] ~ result$time[0:i], ylim = c(0,1), xlim = xlim, type = "l")
		#try(axis(4, at = mean(result$rho_plus[t_2]), labels = FALSE, col = "grey90" ))
		#mtext(round(stability, digits = 5), side = 3)
		#plot(x_new, add = TRUE)
		#}
		
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
			starting = parms$starting, 
			assoc = parms$assoc,
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


write.csv(output, file = paste("data/sim", sim,"/", "output.csv", sep = "") )



stopCluster(cl)







