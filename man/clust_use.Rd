% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/usescape_PACKAGE.R
\name{clust_use}
\alias{clust_use}
\title{Function for clustering intensity of use metrics}
\usage{
clust_use(stck, col = c(1, 2, 3, 4, 6, 7), nb_clust = 1:5, min_fix = 3)
}
\arguments{
\item{stck}{stack that is produced from timing2stack}

\item{col}{specific parameters to be used ("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval"). By default, the function exclude SD Duration and SD interval.}

\item{nb_clust}{select the range of possible clusters that may be identified}

\item{min_fix}{set the minimum number of fixes within a cell for it to be included in the clustering (3 is the minimum when considering cv_duration and cv_interval)}
}
\value{
return a list with the clustering outputs, a raster of the classification ,and a raster of associated uncertainty.
}
\description{
Apply mixture-model clustering to intensity of use metrics("Number visits", "Total duration", "Mean duration", "CV duration", "SD duration", "Mean interval", "CV interval", "SD interval")
}
\details{
The function applies mixture-model clustering using the mclust package to a stack of intensity of use metrics.
}
\examples{
traj1<-sim_mov(type="OU", npatches=3, grph=T)
timing_ls<-traj2timing(traj1, res=50, grid=NULL)
stck<-timing2stack(timing_ls, col=1,2,3,4,6,7) 
test<-clust_use(stck)
plot(test[[2]], col=c("grey", "blue", "green", "red", "purple", "orange"))
plot(test[[3]]) #Uncertainty 
mask<-mask_uncertain(test, p=0.01)
plot(mask, col=c("grey", "blue", "green", "red", "purple", "orange"))

}
