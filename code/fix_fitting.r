source("code/simfunctions.r")


######## test cases ######################

iteration = 2818
#subset(output, best_model == 2)$ID
i = c(3283, 3284, 3285, 3286, 3335, 3336, 3337, 3387, 3388)


par(mfrow = c(3,3))
for(iteration in i){
load(paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))

dd4 <- do.call("rbind", result$cumpatch)

result$fit <- fitPL(dd4, p_spanning, n = result$out$n)


plot(p ~ size, data = dd4, log = "xy", pch = 20, xlim = c(1,10000), ylim = c(0.001,1) )

x_size <- exp(seq(log(1),log(10000), length = 100))


try(lines( x_size,
       exp(predict(result$fit[[result$out$best+2]], list(size = x_size))),
       col = "black", 
       lwd = 2
))


try(lines( x_size,
           exp(predict(result$fit[[5]], list(size = x_size))),
           col = "red", 
           lwd = 1
))


}


################ starting parallel backend

library(foreach)
library(doSNOW)

workstation <-  list(host = "162.38.184.118", user = "schneider",
                     rscript = "/usr/lib/R/bin/Rscript",
                     snowlib = "/usr/lib/R/library")


#workerlist <- c(rep(list(workstation), times = 22), rep("localhost", times = 10)) 

cl <- makeSOCKcluster(rep("localhost", times = 10), master="162.38.184.88", outfile='out_messages.txt')

#clusterCall(cl,function() Sys.info()[c("nodename","machine")])

registerDoSNOW(cl)



foreach(iteration = iterations$ID, .combine = rbind) %dopar% { 
  #subse

  load(paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))
  
  
    
  ## fit power law models and compare via AIC via function
  
  if(result$out$best_model %in% c(2,3,4)) {
    dd4 <- do.call("rbind", result$cumpatch)
    p_spanning <- mean(unlist(sapply(1:length(result$cumpatch), function(x) tail(result$cumpatch[[x]]$p, 1)) )) 
    
    result$fit <- fitPL(dd4, p_spanning, n = result$out$n)
    
    result$out$best_model = result$fit$best 
  }
  
  
  result$out$class = c("DEG", "DOWN","PL", "UP", "COV")[result$out$best]
  

result$out$predict_largest <- result$out$largest_patch^result$out$alpha
result$out$mismatch_PL <- log(result$out$mean_p_largest / result$out$predict_largest)



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

result$out$updated3 <- TRUE

save(result, file = paste("data/sim", sim,"/results/", "result_", iteration, sep = ""))


order <- c(
  "ID", "assoc", "g", "b", "m0" ,  "n", "runtime", "runtime_sd", "rho_plus_ini",   
  "rho_plus", "rho_plus_sd", "q_plus_plus", "clus_coef", "largest_patch", "best_model", 
  "class",
  "largest_patch_sd", "alpha", "alpha_se", "alpha_p", "trunc", "trunc_se", "trunc_p",
  "lower", "lower_se", "lower_p", "mean_p_largest", "predict_largest", "mismatch_PL")
return(result$out[order])

}-> output

write.csv(output, file = paste("data/sim", sim,"/", "output.csv", sep = "") )

stopCluster(cl)