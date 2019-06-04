#' @title Plot for the nowcast output function
#' @description Make plots to visualize the output of the nowcast function
#' @param out output of the nowcast function.
#' @param type 'fcst', 'factors', 'eigenvalues' or 'eigenvectors'. The 'eigenvalues' and 'eigenvectors' options are only available for the two stages methods.
#' @examples
#' \dontrun{
#' data <- Bpanel(base = USGDP$base,
#'                trans = USGDP$legend$Transformation,
#'                aggregate = FALSE)
#' frequency <- c(rep(12, ncol(data) -1), 4)
#' now2s_agg <- nowcast(formula = RGDPGR ~ ., data = data, r = 2, p = 2, q = 2, 
#'                      method = '2s_agg', frequency = frequency)
#' 
#' nowcast.plot(now2s_agg, type = "fcst")
#' nowcast.plot(now2s_agg, type = "factors")
#' nowcast.plot(now2s_agg, type = "eigenvalues")
#' nowcast.plot(now2s_agg, type = "eigenvectors")
#' }
#' @export

nowcast.plot <- function(out, type = "fcst"){
  

  
  if(type == "fcst"){
    
    data <- data.frame(date = as.Date(out$yfcst), out$yfcst)
    data[max(which(is.na(data$out))),"out"] <- data[max(which(is.na(data$out))),"in."] 
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 5), xpd = F)
    graphics::plot(data[,"y"], xaxt = "n", main = "",  bty = "l",
         col = "#FFFFFF", ylab = "y unit", xlab = "Time")
    graphics::axis(1, at = seq(1,nrow(data),frequency(out$yfcst)), labels = substr(data[seq(1,nrow(data),frequency(out$yfcst)),"date"],1,7), las=1, cex = 0.7)
    graphics::abline(v = seq(1,nrow(data),frequency(out$yfcst)), lty = 3, col = "#D9D9D9")
    graphics::grid(col = "#D9D9D9", nx = NA, ny = NULL)
    graphics::lines(data[,"y"], type = "l", lty = 1, col = 1)
    graphics::lines(data[,"in."], type = "l", lty = 1, lwd = 2, col = "dodgerblue")
    graphics::lines(data[,"out"], type = "l", lty = 4, lwd = 2, col = "orangered")
    graphics::par(xpd = T)
    add_legend("topright", legend=c("y","yhat","fcst"), bty = "n",
           lty = c(1,1,4), lwd = c(1,2,2), col = c(1,"dodgerblue","orangered"))
    graphics::title(main = list("Forecasting", font = 1, cex = 0.9))

  }else if(type == "eigenvalues"){
  
    if(!("reg" %in% names(out))){
      stop('Graph not available for EM method')
    }
    
    graphics::par(mar = c(5.1,4.1,4.1,2.1), xpd = F)
    eig <- out$factors$eigen$values/sum(out$factors$eigen$values)*100
    n <- min(20,length(eig))
    graphics::barplot(eig[1:n], col = "#ADD8E6", border = "steelblue", names.arg = 1:n, ylim = c(0, seq(0,100,5)[min(which(!(max(eig[1:n]) > seq(0,100,5))))]),
            xlab = "eigenvalues", ylab = "%")
    graphics::grid(col = "#D9D9D9")
    graphics::title(main = list("eigenvalues: percentage variance", font = 1, cex = 0.9))

  }else if(type == "eigenvectors"){
    
    if(!("reg" %in% names(out))){
      stop('Graph not available for EM method')
    }
    
    graphics::par(mar = c(5.1,4.1,4.1,5), xpd = F)
    vec <- out$factors$eigen$vectors[,1]
    pvec <- (out$factors$eigen$vectors[,1]^2)*100
    color <- ifelse(vec >= 0, "dodgerblue", "orangered")
    graphics::plot(pvec, main = "",  bty = "l", xaxt = "n", type = "h", ylab = "weight (%)", xlab = "variable", col = color)
    graphics::axis(1, at = seq(1,length(vec),1), labels = seq(1,length(vec),1), las=1, cex = 0.7)
    graphics::title(main = list("Variable Percentage Weight in Factor 1", font = 1, cex = 0.9))
    graphics::par(xpd = T)
    graphics::text(y = max(pvec), x = length(pvec)*1.1, labels = "signal weights:", col = 1, cex = 0.8)
    graphics::text(y = max(pvec)*0.94, x = length(pvec)*1.1, labels = "positive", col = "dodgerblue", cex = 0.8)
    graphics::text(y = max(pvec)*0.88, x = length(pvec)*1.1, labels = "negative", col = "orangered", cex = 0.8)

    
  }else if(type == "factors"){
    
    graphics::par(mar=c(5.1, 4.1, 4.1, 5), xpd = F)
    n <- ncol(data.frame(out$factors$dynamic_factors))
    stats::ts.plot(out$factors$dynamic_factors, col = c(1,"orangered","blue"), lty = c(1,2,3), gpars = list(bty = "l"))
    anos <- unique(as.numeric(substr(as.Date(out$factors$dynamic_factors),1,4)))
    graphics::grid()
    graphics::par(xpd = T)
    graphics::title(main = list("Estimated Factors", font = 1, cex = 1))
    
    if("month_y" %in% names(out)){
      leg <- colnames(out$factors$dynamic_factors)
    }else{
      leg <- paste("Factor", 1:n)
    }
    add_legend("topright", legend = leg, bty = "n",
           col = c(1,"orangered","blue"), lty = c(1,2,3), lwd = c(1,1,1,2,2,2), cex = 0.9)
    
  }
}

