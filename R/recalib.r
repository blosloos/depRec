#' @title Function for m/z and RT recalibration.
#'
#' @export
#'
#' @description Function for m/z and RT recalibration.
#'
#' @param peaklist Matrix or data.frame of sample peaks. 3 columns: m/z, intensity, RT.
#' @param mz Numeric vector. Masses to be matched within large tolerances of \code{tolmz}.
#' @param tolmz Numeric. +/- m/z tolerance (precision)
#' @param ppm Logical. \code{tolmz} given in ppm?
#' @param ret Numeric vector. Retention times to be matched; units equal to those of the input files
#' @param tolret Numeric. RT tolerance for matches with \code{peaklist}; units equal to those of the input files
#' @param what Character strings "mass" or "ret".
#' @param one Logical. Only use recalibration \code{mz, ret} that can be matched unambiguously.
#' @param knot Integer. Number of spline knots.
#' @param plot_it Logical. Produce recalibration plot?
#' @param path_1 Logical \code{FALSE} or character string. If not \code{FALSE}, filepath to output the plot
#' @param path_2 Logical \code{FALSE} or character string. Filepath for saving GAM model.
#' @param stopit Logical. Triggers a full R error (\code{TRUE}) or returns an error message string (\code{FALSE}) if failing. 
#' @param intermediate_results Logical. Call by reference to have intermediates in the enclosing environment?
#' @param plot_ppm FALSE or vector of numerics. If a vector, plots the ppm mass deviations from the vector.
#' @param max_recal FALSE or numeric. If the recalibration proposes mass corrections larger than that value (either ppm or absolute, see parameter \code{ppm}), dismiss the recalibration.
#' @param up_bound FALSE or numeric. Maximum mass differences to be considered.
#' @param low_bound FALSE or numeric. Minimum mass differences to be considered.
#'
#' @return Recalibrated \code{peaklist}.
#' 
#' @details enviMass workflow function. \code{mz} and \code{ret} must be of same length and specify the points to recalibrate with.
#' A minimum of 15 matches must be found to conduct the recalibration.
#'


recalib <- function(
  peaklist,
  mz,
  tolmz = 10,
  ppm = TRUE,
  ret,
  tolret,
  what = "mass",
  one = TRUE,
  knot = 5,
  plotit = FALSE,
  path_1 = FALSE,
  path_2 = FALSE,  
  stopit = FALSE,
  intermediate_results = FALSE,
  plot_ppm = FALSE,
  max_recal = FALSE,
  up_bound = FALSE,
  low_bound = FALSE
){


  ##############################################################################
  # search for concerned peaks #################################################
  if(what != "mass" && what != "ret") stop("what what?")
  peaks <- search_peak(peaklist, mz, tolmz, ppm, ret, tolret, onlymax = one);
  # collect concurring peaks ###################################################
  getit1 <- c();  # expected mz
  getit2 <- c();  # observed peaklist
  for(i in 1:length(peaks)){
    if(peaks[i]!="FALSE"){
      put <- as.numeric(strsplit(peaks[i],"/")[[1]]);
      if(what=="mass"){
		getit1<-c(getit1,rep(mz[i],length(put)));  # expected
		getit2<-c(getit2,peaklist[put,1]);         # observed
      }else{ # what=="ret"
		getit1<-c(getit1,rep(ret[i],length(put))); # expected
		getit2<-c(getit2,peaklist[put,3]);         # observed
      }
    }
  }
  getit3 <- c(getit1 - getit2); # observed-expected
  if(!is.logical(up_bound)){getit2 <- getit2[getit3 <= (up_bound[1] / 1000)]; getit3 <- getit3[getit3 <= (up_bound[1] / 1000)]}
  if(!is.logical(low_bound)){getit2 <- getit2[getit3 >= (low_bound[1] / 1000)]; getit3 <- getit3[getit3 >= (low_bound[1] / 1000)]}
  ##############################################################################
  # stop if too few data found! ################################################
  if(length(getit3) < 15){
		if(path_1!="FALSE"){png(filename = path_1, bg = "white")}    
		plot.new();
		plot.window(xlim=c(1,1),ylim=c(1,1));
		box();
		text(1,1,label="not available",cex=1.5,col="darkred") 
		if(path_1!="FALSE"){dev.off()}   
		if(!stopit){
			return("Too few data points for fit!\n");
		}else{
			stop("Too few data points for fit!\n")
		}
  }
  ##############################################################################
  # train gam model ############################################################
  that<-data.frame(getit2, getit3)
  names(that)<-c("obs","delta")
  #attach(that)
  model <- mgcv::gam(delta ~ s(obs, bs = "ts", k = knot), data = that);
  if(plotit==TRUE){
    if(what=="mass"){    
		if(path_1!="FALSE") png(filename = path_1, bg = "white")
		ylim <- c(min(getit3) * 1000, max(getit3) * 1000)
		if(ylim[1] > 0) ylim[1] <- 0
		if(ylim[2] < 0) ylim[2] <- 0
		plot(getit2, getit3 * 1000, pch = 19, cex = 0.5, xlab = "m/z", ylab = "Expected m/z - observed m/z [mmu]",
			main = "Recalibration results", ylim = ylim);
		abline(h = 0, col = "darkgreen");
		points(getit2[order(getit2)], predict(model)[order(getit2)] * 1000, col = "red", type = "l", lwd = 2);
		if(plot_ppm[1]!="FALSE"){
			ppm_mass<-seq(0,max(getit2),10)
			for(k in 1:length(plot_ppm)){
				ppm_ppm<-(ppm_mass*plot_ppm[k]/1E6*1000)
				lines(ppm_mass,ppm_ppm,lty=2,col="gray")
				lines(ppm_mass,-ppm_ppm,lty=2,col="gray")
				plotmass<-median(ppm_mass)
				text(
					plotmass,ppm_ppm[ppm_mass==plotmass],
					labels=paste(as.character(plot_ppm[k]),"ppm"),
					col="gray",cex=1.2
				)
				text(
					plotmass,-ppm_ppm[ppm_mass==plotmass],
					labels=paste("-",as.character(plot_ppm[k])," ppm",sep=""),
					col="gray",cex=1.2
				)
			}
		}
		if(max_recal != "FALSE"){
			if(ppm){
				ppm_mass<-seq(0,max(getit2),10)
				ppm_ppm<-(ppm_mass*max_recal/1E6*1000)
				lines(ppm_mass,ppm_ppm,lty=2,lwd=1.5,col="red")
				lines(ppm_mass,-ppm_ppm,lty=2,lwd=1.5,col="red")			
			}else{
				abline(h=max_recal,lty=2,lwd=1.5,col="red")
				abline(h=-max_recal,lty=2,lwd=1.5,col="red")				
			}
		}
		if(up_bound != "FALSE"){
			if(ppm){
				ppm_mass <- seq(0, max(getit2), 10)
				ppm_ppm <- (ppm_mass * up_bound / 1E6 * 1000)
				lines(ppm_mass, ppm_ppm, lty = 2, lwd = 1.5, col = "red")		
			}else{
				abline(h = up_bound, lty = 2, lwd = 1.5, col = "red")			
			}
		}		
		if(low_bound != "FALSE"){
			if(ppm){
				ppm_mass <- seq(0, max(getit2), 10)
				ppm_ppm <- (ppm_mass * low_bound / 1E6 * 1000)
				lines(ppm_mass, ppm_ppm, lty = 2, lwd = 1.5, col = "red")		
			}else{
				abline(h = low_bound, lty = 2, lwd = 1.5, col = "red")			
			}
		}		
		if(path_1!="FALSE") dev.off()   
    }else{
		if(path_1!="FALSE") png(filename = path_1, bg = "white")  
		plot(getit2,getit3,pch=19,cex=0.5,xlab="Retention time",ylab="Expected RT - observed RT",
			main="Recalibration results");
		abline(h=0,col="red");
		points(getit2[order(getit2)],predict(model)[order(getit2)],col="red",type="l",lwd=2);
		if(path_1!="FALSE") dev.off()      
    }
  }
  ##############################################################################
  # predict -> recalibrate peaklist ############################################
  if(what=="mass"){
  		that <- data.frame("obs" = peaklist[,1], "delta" = peaklist[,1]);
		pred2 <- mgcv::predict.gam(model,newdata=that);
		newpeaks <- peaklist;
  		if(max_recal == "FALSE"){
			newpeaks[,1] <- c(peaklist[,1] + pred2);		
		}else{
			if(ppm){
				if(!any((pred2 / that[,1] * 1E6) > max_recal)){
					newpeaks[,1] <- c(peaklist[,1] + pred2);
				}else cat("\n recalibration skipped - correction off limits!")		
			}else{
				if(!any(abs(pred2) > max_recal)){
					newpeaks[,1] <- c(peaklist[,1] + pred2);
				}else cat("\n recalibration skipped - correction off limits!")
			}
		}
  }else{
		that <- data.frame("obs" = peaklist[,3], "delta" = peaklist[,3]);
		pred2 <- mgcv::predict.gam(model, newdata = that);
		newpeaks <- peaklist;
		newpeaks[,3] <- c(peaklist[,3] + pred2);
  }
  if(path_2 != "FALSE") save(model, file = path_2)
  if (intermediate_results) {
    # simulate call by reference here:
    imr <- list()
    imr[[1]] <- model;
    imr[[2]] <- peaks;
    imr[[3]] <- getit2;
    imr[[4]] <- getit3;
    names(imr) <- c("model", "matches", "x", "delta_y");
    eval.parent(substitute(intermediate_results <- imr))
  }
  return(newpeaks)
  ##############################################################################
}