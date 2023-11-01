#' @title Simulation of patch-based movement trajectory
#' @description Simulate a movement trajectory with a user defined number of patches and interpatch movement.
#' @details See moveNT package 
#' 
#' @param type whether movement within patches should be based on a 2states process (from package moveHMM) or a Bivariate Ornstein-Uhlenbeck process (OU) (from package adehabitatLT)
#' @param npatches Number of patches, default=5
#' @param ratio Ratio (in percent) of locations associated to interpatch movement, default=5
#' @param nswitch Number of switch/depart from patches, default=150
#' @param ncore Number of locations within a patch per visit, default=200
#' @param spacecore Minimum distance between center of patches, default=200
#' @param seq_visit Specify the sequence of visit among patches, default is random sequence
#' @param stepDist Distribution for step length if 2states specified in type, see simData of moveHMM package
#' @param angleDist Distribution for turn angle if 2states specified in type, see simData of moveHMM package
#' @param stepPar Parameters for step length distribution if 2states specified in type, see simData of moveHMM package
#' @param anglePar Parameters for turn angle distribution if 2states specified in type, see simData of moveHMM package
#' @param s Parameters for the OU process, see simm.mou of adehabitatLT package
#' @param grph Whether a graph of the trajectory should be produced, default=F
#' @keywords traj2adj adj2stack #check with GBR
#' @return A ltraj (adehabitatLT) object
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' traj2<-sim_mov(type="2states", npatches=2, grph=T)
#' 
#' @export 

sim_mov<-function(type=c("2states", "OU"), npatches=5, ratio=5, nswitch=150, ncore=200,spacecore=200, seq_visit=sample(1:npatches, nswitch, replace=T),
                  stepDist= "gamma", angleDist = "vm",  stepPar = c(0.5,3,1,5), anglePar = c(pi,0,0.5,2), s=diag(40,2), grph=F) {
  
  coordx<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  coordy<-sample(seq(0,20,2), npatches, replace=F)*spacecore
  nmig=ncore/ratio
  out<-data.frame()
  for (i in 1:(nswitch-1)){
    
    if(type=="2states") {
      core<-moveHMM::simData(nbAnimals=1,nbStates=2,stepDist=stepDist,angleDist=angleDist,stepPar=stepPar, anglePar=anglePar,zeroInflation=F,obsPerAnimal=ncore)
      corex<-core$x+coordx[seq_visit[i]]
      corey<-core$y+coordy[seq_visit[i]]
      Corri1<-rep(2, ncore)
    }
    
    if(type=="OU") {
      core<-adehabitatLT::simm.mou(date=1:ncore, b=c(coordx[seq_visit[i]],coordy[seq_visit[i]]), s=s)
      corex<-ld(core)$x
      corey<-ld(core)$y
      Corri1<-rep(2, ncore)
    }
    
    if(seq_visit[i] != seq_visit[i+1]) {
      mig<-adehabitatLT::simm.bb(date=1:nmig, begin=c(tail(corex,1), tail(corey,1)), end=rnorm(2, c(coordx[seq_visit[i+1]],coordy[seq_visit[i+1]]), sd=25))
      Corri2<-rep(1, nmig)
      sub<-cbind(c(corex, ld(mig)$x), c(corey, ld(mig)$y), c(Corri1, Corri2))
      
    }
    if(seq_visit[i] == seq_visit[i+1]) {
      sub<-cbind(corex, corey, Corri1)
      colnames(sub)<-c("V1", "V2", "V3")
    }
    out<-rbind(out, sub)
  }
  names(out)<-c("x", "y", "Corri")
  out<-adehabitatLT::as.ltraj(out[,1:2], as.POSIXct(1:nrow(out), origin = "1960-01-01", tz="GMT"), id="id", infolocs=data.frame(out$Corri))
  if(grph==T) {plot(out)}
  return(out)
}


#' @title Calculation of timing history from movement data
#' @description Transform an ltraj object to a timing history (entry and exit time in each pixel) using a user-specified grid size.
#' @details Using a specific grid cell size, this function return the timing history (entry and exit time in each pixels) based on an individual trajectory.
#' 
#' @param mov Movement trajectory, need to be a ltraj object
#' @param res Grid size (based on coordinate system of movement trajectory)
#' @param grid User specified grid (a raster), needs to have a larger extent than the movement trajectory
#' @keywords timing
#' @return A list of objects containing the timing history and the grid use.
#'
#' @examples
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' timing_ls<-traj2timing(traj1, res=50, grid=NULL)
#'
#' @export

#First main function, calculate time in and out of each pixels
traj2timing<-function(mov, res=100, grid=NULL) {
   mov<-adehabitatLT::ld(mov)
  mov[,13]<-1:nrow(mov)
  tt<-sp::SpatialPoints(mov[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  if(is.null(grid)){ras<-raster::raster(xmn=floor(tt1[1])-2*res, ymn=floor(tt1[2])-2*res,xmx=ceiling(tt2[1])+2*res, ymx=ceiling(tt2[2])+2*res, res=res)}
  if(!is.null(grid)){ras<-grid}
  values(ras)<-1:ncell(ras)
  mov$pix_start<-raster::extract(ras,tt)
  mov<-mov[!is.na(mov$dist),]
  timing<-list()
  timing[[ncell(ras)]]<-NA
  mm<-max(table(mov$pix_start))+1
  timing<-lapply(timing, function(x) data.frame(time_in=as.POSIXct(rep(NA,mm)), time_out=as.POSIXct(rep(NA,mm))))
  
  timing[[mov$pix_start[1]]]$time_in[2]<-mov$date[1] 
  
 for (i in 2:nrow(mov)) {
 
  if(mov$pix_start[(i-1)]!=mov$pix_start[i]) {
     
  a<-sum(!is.na(timing[[mov$pix_start[i-1]]]$time_out))+2
  b<-sum(!is.na(timing[[mov$pix_start[i]]]$time_in))+2
  
 
  if(mov$dist[i-1]<=res) {fract<-0.5} 
  if(mov$dist[i-1]>res) {fract<-res/mov$dist[i-1]/2}
  
   timing[[mov$pix_start[i-1]]]$time_out[a]<-mov$date[i-1]+mov$dt[i-1]*fract
   timing[[mov$pix_start[i]]]$time_in[b]<-mov$date[i]-mov$dt[i-1]*fract 
  }
  }
  timing<-lapply(timing, na.omit)
  return(list(timing, ras))
}


### Intermediary functions 
      #add the time columns
      time_diff<-function(x, unit="secs") {
        x$total_time<-as.numeric(difftime(x$time_out, x$time_in, units=unit))
        return(x)}
      
      
      #add an interval 4th columns 
      time_interval<-function(x, unit="secs") {
        time2<-c(x$time_in[-1], NA)
        x$time_interval<-as.numeric(difftime(time2,x$time_out, units=unit))
        return(x)
      }
      
      #frequency (nrow)
      freq_visit<-function(x) {return(nrow(x))}
      
      #total time (sum 3rd column)
      total_duration<-function(x) {return(ifelse(nrow(x)>0, sum(x$total_time, na.rm=T), 0))}
      
      #avg duration (mean 3rd column)
      mean_duration<-function(x) {return(ifelse(nrow(x)>0, mean(x$total_time, na.rm=T), 0))}
      
      #cv duration (cv 3rd column)
      cv_duration<-function(x) {return(ifelse(nrow(x)>0, cv(x$total_time, na.rm=T), 0))}
      sd_duration<-function(x) {return(ifelse(nrow(x)>0, sd(x$total_time, na.rm=T), 0))}
      
      
      #avg interval
      mean_interval<-function(x) {return(ifelse(nrow(x)>0, mean(x$time_interval, na.rm=T), 0))}
      #cv interval 
      cv_interval<-function(x) {return(ifelse(nrow(x)>0, cv(x$time_interval, na.rm=T), 0))}
      
      #sd interval 
      sd_interval<-function(x) {return(ifelse(nrow(x)>0, sd(x$time_interval, na.rm=T), 0))}


      
#' @title Timing history to raster stack conversion 
#' @description Extract use metrics from timing history
#' @details This function extract the intensity of use metrics from the timing history and return a stack object of a raster layer for each metric. 
#' 
#' @param timing_ls timing history for each pixel
#' @param unit_time time unit for timing_ls parameter
#' @keywords 
#' @returns a raster stack object 
#' @examples  
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' timing_ls<-traj2timing(traj1, res=50, grid=NULL)
#' stck<-timing2stack(timing_ls) 
#' plot(stck)
#' 
#' @export 

timing2stack<-function(timing_ls, unit_time="secs") {
  timing<-timing_ls[[1]]
  ras<-timing_ls[[2]]
  timing<-lapply(timing, function(x) time_diff(x, unit=unit_time))
  timing<-lapply(timing, function(x) time_interval(x, unit=unit_time))
  grid<-stack(ras,ras,ras,ras,ras,ras,ras,ras)
  values(grid[[1]])<-unlist(lapply(timing, freq_visit))
  values(grid[[2]])<-unlist(lapply(timing, total_duration))
  values(grid[[3]])<-unlist(lapply(timing, mean_duration))
  values(grid[[4]])<-unlist(lapply(timing, cv_duration))
  values(grid[[5]])<-unlist(lapply(timing, sd_duration))
  values(grid[[6]])<-unlist(lapply(timing, mean_interval))
  values(grid[[7]])<-unlist(lapply(timing, cv_interval))
  values(grid[[8]])<-unlist(lapply(timing, sd_interval))
    names(grid)<- c("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval")
  return(grid)
  }


#' @title Test several pixel sizes (resolution) and estimate the coefficient of variation in residency time
#' @description Test several pixel sizes (resolution) and estimate the coefficient of variation in residency time
#' @details The function estimates the coefficient of variation in residency time calculated over grid of varying pixel size and produce a graph showing their values   
#' 
#' @param mov simulated movement, replace with trajectory object of interest
#' @param res_seq a vector of pixel size. 
#' @param unit_time time unit for parameters
#' @keywords 
#' @returns a graph showing the coefficient of variation in residency time as a function of pixel size.  
#' @examples code examples using function - these are vital ###########################
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' res_test(traj1, res_seq=c(50,100,150))
#' @export 


res_test<-function(mov, res_seq=c(50,100,150), unit_time="secs") {
  ls<-pbapply::pblapply(res_seq, function (x) traj2timing(mov, x))
  timing<-lapply(ls, function(x) x[[1]])
  timing<-lapply(1:length(timing), function(y) lapply(timing[[y]], function(x) time_diff(x, unit=unit_time)))
  rt<-lapply(1:length(timing), function(y) unlist(lapply(timing[[y]], total_duration)))
  rt_var<-lapply(rt, function(x) log(cv(x[x>0])))
  #rt_var<-lapply(rt, function(x) log(var(x)))
  plot(res_seq, rt_var, xlab="Resolution", ylab="CV RT", type="l")
}
  
#' @title Function for clustering intensity of use metrics 
#' @description Apply mixture-model clustering to intensity of use metrics("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval")
#' @details The function applies mixture-model clustering using the mclust package to a stack of intensity of use metrics. 
#' 
#' @param stck  stack that is produced from timing2stack
#' @param col specific parameters to be used ("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval"). By default, the function exclude SD Duration and SD interval. 
#' @param nb_clust select the range of possible clusters that may be identified 
#' @param min_fix set the minimum number of fixes within a cell for it to be included in the clustering (3 is the minimum when considering cv_duration and cv_interval)
#' @keywords 
#' @returns return a list with the clustering outputs, a raster of the classification ,and a raster of associated uncertainty. 
#' @examples 
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' timing_ls<-traj2timing(traj1, res=50, grid=NULL)
#' stck<-timing2stack(timing_ls, col=1,2,3,4,6,7) 
#' test<-clust_use(stck)
#' plot(test[[2]], col=c("grey", "blue", "green", "red", "purple", "orange"))
#' plot(test[[3]]) #Uncertainty 
#' mask<-mask_uncertain(test, p=0.01)
#' plot(mask, col=c("grey", "blue", "green", "red", "purple", "orange"))
#' 
#' @export

clust_use<-function(stck,col=c(1,2,3,4,6,7),nb_clust=1:5, min_fix=3 ) {
  if(require("mclust")){
    print("mclust is loaded correctly")
  } else {
    print("trying to install mclust")
    install.packages("mclust")
    if(require(mclust)){
      print("mclust installed and loaded")
    } else {
      stop("could not install mclust")
    }
  }
data<-values(stck)  
data<-data[,col]
data<-scale(data[data[,1]>=min_fix,]) ## Need to be visited at least twice to be considered or 3 times for CV/SD interval included 
clust<-Mclust(data, G=nb_clust)
print(clust$parameters$mean)
print(clust$parameters$pro)
ras1<-stck[[1]]
values(ras1)[values(ras1)<min_fix]<-0
values(ras1)[values(ras1)>=min_fix]<-clust$classification
ras2<-stck[[1]]
values(ras2)[values(ras2)<min_fix]<-0
values(ras2)[values(ras2)>=min_fix]<-clust$uncertainty
out<-list(clust, ras1, ras2)
return(out)
}


#' @title Function for scaling back centers (means) of clusters to actual values 
#' @description Convert the centers (means) of clusters to original scale
#' @details The function reverse-standardize (*sd+mean) the values of the centers of each cluster. 
#' 
#' @param stck  stack that is produced from timing2stack
#' @param clust output produced by clust_use 
#' @param col specific parameters to be used ("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval"). By default, the function exclude SD Duration and SD interval. 
#' @param nb_clust select the range of possible clusters that may be identified 
#' @param min_fix set the minimum number of fixes within a cell for it to be included in the clustering (3 is the minimum when considering cv_duration and cv_interval)
#' @keywords 
#' @returns return matrix of the cluster centers on their original sclaes 
#' @examples 
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' timing_ls<-traj2timing(traj1, res=50, grid=NULL)
#' stck<-timing2stack(timing_ls, col=c(1,2,3,4,6,7), col=c(1,2,3,4,6,7), min_fix=3) 
#' test<-clust_use(stck)
#' backscaling_clust_use(stck, test, col=c(1,2,3,4,6,7), min_fix=3)
#' 
#' @export
backscaling_clust_use<-function(stck, clust, col=c(1,2,3,4,6,7), min_fix=3) {
  data<-values(stck)  
  data<-data[,col]
  avg<-colMeans(data[data[,1]>=min_fix,])
  sd<-apply(data[data[,1]>=min_fix,], 2, sd)
  tt<-clust[[1]]$parameters$mean*sd+avg
  print(tt)
  return(tt)
}
  


#' @title Mask uncertain pixel from classification
#' @description Remove pixel with uncertainty higher than a specific value from a plot
#' @details Remove pixel with uncertainty higher than a specific value from a plot
#' 
#' @param out  output of clust_use
#' @param p uncertainty threshold (values with uncertainty higher than p will be masked). Default = 0.05 
#' @keywords 
#' @returns return a raster object
#' @examples 
#' traj1<-sim_mov(type="OU", npatches=3, grph=T)
#' timing_ls<-traj2timing(traj1, res=50, grid=NULL)
#' stck<-timing2stack(timing_ls, col=1,2,3,4,6,7) 
#' test<-clust_use(stck)
#' plot(test[[2]], col=c("grey", "blue", "green", "red", "purple", "orange"))
#' plot(test[[3]]) #Uncertainty 
#' mask<-mask_uncertain(test, p=0.01)
#' plot(mask, col=c("grey", "blue", "green", "red", "purple", "orange"))
#' 
#' @export

mask_uncertain<-function(out, p=0.05) {
  ras<-out[[2]]
  values(ras)[values(out[[3]])>p]<-0
  return(ras)
}



#' Looping over all individuals
#'
#' Extract the timing history and intensity of use metrics for all individuals in a traj object
#' @param traj An object produce by the function adehabitatLT with multiple individuals
#' @param res Grid size, will be apply to all individuals
#' @keywords timing2stack traj2timing
#' @return A list object containing a raster stack object for each individual
#' @export
#' @examples
#' data(puechcirc)
#' traj<-na.omit(puechcirc)
#' ls1<-loop_id(traj, res=300)
#' table<-table_cluster(traj, ls1)
#' ind<-ind_clust(table)
#' pop<-pop_clust(traj, ind)
#' stack<-clust_stack(ls1, pop, ind, table, min_fix = 3)
#' plot(stack[[1]][[1]]) #Plot first individuals 
#' plot(stack[[2]][[1]]) #Plot second individuals 
loop_id<-function(traj, res=100){
  tt<-SpatialPoints(ld(traj)[,1:2])
  tt1<-apply(coordinates(tt), 2, min)
  tt2<-apply(coordinates(tt), 2, max)
  ras<-raster(xmn=floor(tt1[1])-2*res, ymn=floor(tt1[2])-2*res,xmx=ceiling(tt2[1])+2*res, ymx=ceiling(tt2[2])+2*res, res=res)
  id<-unique(adehabitatLT::id(traj))
  id2<-adehabitatLT::id(traj)
  out<-list()
  for (i in 1:length(id)) {
    try(out[[i]]<-timing2stack(traj2timing(traj[which(id2==id[i])], res=res, grid=ras)))
    cat(id[i], '\n')
  }
  
  return(out)
}

#' Convert a list of timing2stack object to a data.frame for clustering
#'
#' Convert output of loop_id function to a data.frame.
#' @param traj The trajectory used in loop (a traj object)
#' @param ls The output of the loop_id function
#' @keywords timing2stack traj2timing loop_id
#' @return A data.frame object.
#' @export
#' @examples
#' data(puechcirc)
#' traj<-na.omit(puechcirc)
#' ls1<-loop_id(traj, res=300)
#' table<-table_cluster(traj, ls1)
#' ind<-ind_clust(table)
#' pop<-pop_clust(traj, ind)
#' stack<-clust_stack(ls1, pop, ind, table, min_fix = 3)
#' plot(stack[[1]][[1]]) #Plot first individuals 
#' plot(stack[[2]][[1]]) #Plot second individuals  
table_cluster<-function(traj, ls) {
  id<-unique(adehabitatLT::id(traj))
  if(length(id)!=length(ls)) {stop("traj and grid don't have the same number of individuals")}
  out<-data.frame()
  for (i in 1:length(id)) {
    tt1<-data.frame(values(ls[[i]]))
    #tt3<-data.frame(na.omit(tt1))
    tt1$ID<-id[i]
    out<-rbind(out, tt1)
  }
  return(out)
}

#' Individual-level clustering of intensity of use metrics
#'
#' Perform individual-level clustering (first step) of intensity of use metrics. This function uses the output of table_cluster and perform a mixture-model. Users can select which variables will be used and the maximum number of clusters. See also mclust
#' @param table An output from the table_cluster function
#' @param nb_clust The number of clusters to be tested, see the documentation for mclust for more information. Default = 1:5.
#' @param min_fix The minimum of locations in a pixel needed to be included in the analysis. Default=3. 
#' @param vars The variable to be included. Default = c("Number.visits", "Total.duration", "Mean.duration", "CV.duration", "Mean.interval", "SD.interval")
#' @keywords timing2stack traj2timing loop_id table_cluster
#' @return A list object with each element representing an individual.
#' @export
#' @examples
#' data(puechcirc)
#' traj<-na.omit(puechcirc)
#' ls1<-loop_id(traj, res=300)
#' table<-table_cluster(traj, ls1)
#' ind<-ind_clust(table)
#' pop<-pop_clust(traj, ind)
#' stack<-clust_stack(ls1, pop, ind, table, min_fix = 3)
#' plot(stack[[1]][[1]]) #Plot first individuals 
#' plot(stack[[2]][[1]]) #Plot second individuals 
ind_clust<-function(table, nb_clust=1:5, min_fix=3, vars=c("Number.visits", "Total.duration", "Mean.duration", "CV.duration", "Mean.interval", "SD.interval")) {
  if(require("mclust")){
    print("mclust is loaded correctly")
  } else {
    print("trying to install mclust")
    install.packages("mclust")
    if(require(mclust)){
      print("mclust installed and loaded")
    } else {
      stop("could not install mclust")
    }
  }
  id<-unique(table$ID)
  ls<-list()
  
  for (i in 1:length(id)) {
    tt<-table[table$ID==id[i],vars]
    tt<-scale(tt[tt[,1]>=min_fix,])
    try(ls[[i]]<-Mclust(tt, G=nb_clust))
    print(id[i])
  }
  return(ls)
}

#' Population-level clustering of intensity of use metrics
#'
#' Combine individual-level clustering of movement metrics into a population-level clustering (second step). Users can  define the number of clusters. See also mclust
#' @param traj The trajectory object
#' @param ls Individual-level clustering object, the output of ind_clust.
#' @param n_clust The number of clusters to be tested, see the documentation for mclust for more information. Default = 1:5.
#' @keywords atiming2stack traj2timing loop_id table_cluster ind_clust
#' @return A list object with each element representing an individual.
#' @export
#' @examples
#' data(puechcirc)
#' traj<-na.omit(puechcirc)
#' ls1<-loop_id(traj, res=300)
#' table<-table_cluster(traj, ls1)
#' ind<-ind_clust(table)
#' pop<-pop_clust(traj, ind)
#' stack<-clust_stack(ls1, pop, ind, table, min_fix = 3)
#' plot(stack[[1]][[1]]) #Plot first individuals 
#' plot(stack[[2]][[1]]) #Plot second individuals 
pop_clust<-function(traj, ls, nb_clust=1:5) {
  id<-unique(adehabitatLT::id(traj))
  nvar<-nrow(ls[[1]]$parameters$mean)
  coef<-data.frame()
  for (i in 1:length(id)) {
    G<-ls[[i]]$G
    gg<-data.frame(cbind(t(ls[[i]]$parameters$mean), ls[[i]]$parameters$pro, rep(id[i], G))   )
    coef<-rbind(coef, gg)
  }
  coef[,-ncol(coef)] <- sapply(coef[-ncol(coef)],function(x) as.numeric(as.character(x)))
  names(coef)[nvar+1]<-"Prop"
  names(coef)[nvar+2]<-"ID"
  clust<-Mclust(coef[,1:nvar], G=nb_clust)
  coef$clust<-clust$classification
  print(clust$parameters$mean)
  print(clust$parameters$pro)
  out<-list(clust, coef)
  return(out)
}

#' Back-association of population-level clustering to individual clusters
#'
#' Generate individual-level rasters of population-level clustering. For each individual, the function generates a raster stack containing a raster of the most likely cluster, and several rasters giving the probability of observing each cluster.
#' @param grid The output of the loop function
#' @param pop_clust The output of the pop_clust function
#' @param ind_clust The output of the ind_clust function
#' @param table The output of table_cluster
#' @keywords timing2stack traj2timing loop_id table_cluster ind_clust pop_clust
#' @return A list of raster stack object.
#' @export
#' @examples
#' data(puechcirc)
#' traj<-na.omit(puechcirc)
#' ls1<-loop_id(traj, res=300)
#' table<-table_cluster(traj, ls1)
#' ind<-ind_clust(table)
#' pop<-pop_clust(traj, ind)
#' stack<-clust_stack(ls1, pop, ind, table, min_fix = 3)
#' plot(stack[[1]][[1]]) #Plot first individuals 
#' plot(stack[[2]][[1]]) #Plot second individuals 
clust_stack<-function(ls1, pop, ind, table, min_fix=3) {
  id<-unique(table$ID)
  coef2<-cbind(pop[[2]], pop[[1]]$z)
  out_ls<-list()
  n.clust1<-length(unique(coef2$clust))
  
  for (i in 1:length(ls1)) {
    coef3<-coef2[coef2$ID==id[i],]
    n.clust<-length(coef3$clust)
    cl<-((coef3$clust))
    class<-data.frame(cbind(1:length(cl), cl))
    class<-rbind(class, c(0,0))
    prop<-coef3[,as.character(cl)]
    
    out2<-table[table$ID==id[i],]
    out2$clust<-0
    out2$clust[out2$Number.visits>=min_fix]<-ind[[i]]$classification
    out2$id<-1:nrow(out2)
    out2<-merge(out2, class, by.x="clust", by.y="V1", sort=F)
    out2<-out2[order(out2$id),]
    
    nc<-ncol(ind[[i]]$z)
    mat<-data.frame(matrix(1, nrow=nrow(out2), ncol=nc))
    mat[out2$Number.visits>=min_fix,]<-ind[[i]]$z
    out2<-cbind(out2,mat)
    
    for (j in 1:nrow(out2)) {
      if(out2$clust[j]!=0) {out2[j,(ncol(out2)-n.clust+1):ncol(out2)]<- out2[j,(ncol(out2)-n.clust+1):ncol(out2)]*prop[out2$clust[j],]} #Multiply the two probability
    }
    

    r0<-ls1[[i]][[1]]
    values(r0)<-0
    tt <- values(r0)
    gr<-stack(r0)
    tt<-out2$cl
    gr[[1]]<-setValues(gr[[1]],tt)
    
    for (z in 2:(n.clust1+1)) { gr[[z]]<-r0 }
    for (k in 1:n.clust) {
      tt<-out2[,(ncol(out2)-n.clust)+k]
      gg<-setValues(gr[[1]],tt)
      gr[[cl[k]+1]]<-mosaic(gr[[cl[k]+1]], gg, fun=max)
    }
    names(gr)<-c("Clust", paste("Prop", paste("Clust", 1:n.clust1, sep="")))
    out_ls[[i]]<-gr
    print(id[i])
  }
  return(out_ls)
}



