### PLOT P MATRIX ####
plot_pmat <- function(mod, covariates = 0) {
  
  pmat <- t(pmatrix.msm(mod, t = 10, covariates = covariates, ci = "none"))
  class(pmat) <- "matrix"
  pmat <- round(pmat*100, 1)
  diag(pmat) <- 0
  
  pos.box <- cbind(c(0.5, 0.05, 0.95, 0.5),
                   c(0.95, 0.5, 0.5, 0.05))
  
  par(mar=c(.5,0.5,.5,.5))
  pm <- plotmat(pmat, pos = pos.box, curve = 0.07, name = states,
                lwd = 1.2, relsize = .9,
                box.cex = 2, cex.txt = 0, txt.col = "white",
                dtext = .3, txt.font = 2,
                box.lwd = 0.1, box.type = "rect", shadow.size = 0.005,
                box.prop = 0.35, box.size = 0.13, box.col = st_col,
                arr.length = pmat/40, arr.width = pmat/40,  arr.type ="triangle",
                arr.col ="grey40", arr.lcol = "grey40",
                arr.lwd = pmat*.8)
  pm$arr$TextX[6] <- pm$arr$TextX[6]-.03
  pm$arr$TextX[7] <- pm$arr$TextX[7]+.03
  
  text(pm$arr$TextX, pm$arr$TextY, pm$arr$Value, cex = 1.5)
}


#### function to plot risk ratio ####

plot_risk <- function(mod, 
                      varnames = c("Température", "CMI", 
                                   "Drainage", "pH", 
                                   "Naturelle 1", "Naturelle 2", 
                                   "Coupe 1", "Coupe 2")) {
  
  # Get estimates, CI and p-value
  HR <- hazard.msm(mod)
  
  n_var <- length(varnames)
  
  trans <- rownames(HR[[1]])
  
  trans_list <- list()
  
  for(i in 1:length(trans)) {
    coef <- lapply(HR, function(x) x[i,])
    coef <- do.call(rbind, coef)
    
    signif <- ifelse(coef[,2] <= 1 & coef[,3] <= 1 | coef[,2] >= 1 & coef[,3] >= 1, 1, 0)
    
    
    # Color
    col_pt <- rep(NA, n_var)
    
    col_pt[which(signif==1 & coef[,1]>1)] <- "dodgerblue4"
    col_pt[which(signif==0 & coef[,1]>1)] <- "#C3D3E2"
    col_pt[which(signif==1 & coef[,1]<1)] <- "#B41414"
    col_pt[which(signif==0 & coef[,1]<1)] <- "#ECC4C4"
    
    trans_list[[trans[i]]] <- cbind.data.frame(coef, signif, col_pt=col_pt)
  }
  
  # Layout
  mat <- matrix(c(0,1,2,0,
                  3,0,4,5,
                  6,7,0,8,
                  0,9,10,0),
                4, 4, byrow = T)
  diag(mat) <- 11:14
  mat <- rbind(mat, 15:18)
  
  layout(mat, heights = c(1,1,1,1,.7))  
  
  # Plot
  par(mar = c(.3,2.5,.3,.5), oma = c(0,2.5,1.5,0))
  for(i in trans) {
    
    tmp <- trans_list[[i]]

    labx <- "Hazard ratio"

    ylim <- range(lapply(trans_list, function(x) range(x[,1:3])))
    
    plot(tmp$HR, log = "y", ylim = ylim, col = "transparent", yaxs = "i",
         ann=F, xaxt="n", yaxt="n",  bty = "l", frame.plot = T)
    
    abline(h = 1, col = "grey65")
    
    # bars for confidence interval
    arrows(x0 = 1:n_var, 
           y0 = tmp$L,
           y1 = tmp$U, 
           angle = 90, code = 1, 
           length = 0, col = "grey45", lwd = 1.8, xpd = NA)
    
    # points
    points(tmp$HR, pch = 19, 
           col = as.character(tmp$col_pt), cex = 2)
    
    # axis
    axis(2, at=seq(0,100,1), tcl= -0.2, labels=F, col = "grey35")
    axis(2, at=c(.1,1,10,100), labels = c(.1,1,10,100), las = 1, cex.axis = 1.4)
  }
  
  # Labels
  par(mar=c(0.5,2.5,0.5,0.5))
  for(st in states) {
    plot0(text = st, fill = "grey95", font = 2, cex = 2.5)
  }
  
  
  par(mar = c(.5,2,.1,.5))
  for(i in 1:4) {
    plot0(x=1:n_var, y=rep(0,n_var))
    text(x=1:n_var, y=rep(1,n_var), labels = varnames, font = 2, cex = 1.6,
         srt = 90, xpd = NA, adj = 1) 
  }
  
  mtext(labx, 2,  line = 1, cex = 1.5, font = 2, las=0, outer=T)
  
  mtext("De", 3, adj = -0.01, line = -2, cex = 1.3, font = 2, outer=T)
  
  mtext("À", 3, adj = 0.06, line = 0, cex = 1.3, font = 2, outer=T)

}


#### steady state ####

steady_state <- function(pmat = NULL, qmat = NULL){
  
  if(!is.null(pmat)) {
    eig <- eigen(t(pmat))
    
    pi <- t(eig$vectors[,which(abs(eig$values-1) < 1e-12)])
  } 
  
  if(!is.null(qmat)) {
    eig <- eigen(t(qmat))
    
    pi <- t(eig$vectors[,which(abs(eig$values) < 1e-12)])
  }
  
  steady <- pi/sum(pi)
  colnames(steady) <- colnames(pmat)
  
  return(steady)
}

### FIND INTERSECT BETWEEN 2 CURVES ####
curve_intersect <- function(curve1, curve2) {
  
  colnames(curve1) <- colnames(curve2) <- c("x", "y")
  # Approximate the functional form of both curves
  curve1_f <- approxfun(curve1$x, curve1$y, rule = 2)
  curve2_f <- approxfun(curve2$x, curve2$y, rule = 2)
  
  # Calculate the intersection of curve 1 and curve 2 along the x-axis
  point_x <- uniroot(function(x) curve1_f(x) - curve2_f(x),
                     c(min(curve1$x), max(curve1$x)))$root
  
  # Find where point_x is in curve 2
  point_y <- curve2_f(point_x)
  
  return(list(x = point_x, y = point_y))
}

#### Plot steady state
plot_SS <- function(logging=0, lang = "fr") {
  par(mar = c(4.2,4.3,.7,9.7))
  
  col_bb <- "#158282"
  col_tt <- "#D43650"
  col_pts <- c("white", "grey", "black")
  
  if(lang == "fr") {
    lgd <- c("Peu", 
             "Modérée", 
             "Majeure")
    xlab <- "Température de la saison de croissance"
    ylab <- "Proportion d'états à l'équilibre"
    stB <- "Boréal"
    stT <- "Tempéré"
  } else {
    lgd <- c("Minor", 
             "Moderate", 
             "Major")
    xlab <- "Growing season temperature"
    ylab <- "Steady state proportion"
    stB <- "Boreal"
    stT <- "Temperate"
  }
 
  
  # Empty plot
  plot0(ylim = c(0,1), xlim = range(x), xaxs = "i", frame.plot = TRUE, yaxs = "i")
  polygon(x = c(tp_mixed, rev(tp_mixed)), y = c(0,0,1,1),
          col = alpha("grey", .2), border = NA)
  axis(1, cex.axis = 1.2)
  axis(2, cex.axis = 1.2, las = 1)
  mtext(xlab, 1, line = 3, cex = 1.8, font = 2)
  mtext(ylab, 2, line = 3, cex = 1.8, font = 2)
  
  text(10.7, .94, stB, col = col_bb, cex = 2, font = 2)
  text(13.8, .94, stT, col = col_tt, cex = 2, font = 2)
  
  
  for(i in 1:length(logging)) {
    ll <- which(df[,"logging"] == logging[i])
    lines(bb[ll] ~ x, col = col_bb, lwd = 3.5, lty = i)
    lines(tt[ll] ~ x, col = col_tt, lwd = 3.5, lty = i)
    
    # Intersect between boreal and mixed+temperate SS curves
    int <- curve_intersect(curve1 = cbind.data.frame(x, bb[ll]), 
                           curve2 = cbind.data.frame(x, tt[ll]))
    points(int$x, 1, xpd = NA, 
           pch = 21, col = "black", bg = col_pts[i], cex = 2.5, lwd = 2)
    
  }
  
  legend(14.35, .65, legend = lgd[1:length(logging)],
         pch = 21, col = "black", pt.bg = col_pts[1:length(logging)], 
         lty = 1:length(logging), lwd = 2.5, seg.len = 3.53, x.intersp = .5,
         cex = 1.4, pt.cex = 2.4, pt.lwd = 2.5, xpd = NA, bty = "n")
}


  

### BARPLOT STEADY STATE & TRANSIENT ####


barplot_index <- function(index = NULL, bars = 1:3, ylim = NULL,
                          ylab = NULL, lgd = NULL,
                          colss = c("#f1ba53","#E38451", "#b5305d"),
                          lang = "fr") {
  par(mar = c(4,4.3,.7,.5))

  if(is.null(lgd)) {
    if(lang == "fr") {
      lgd <- c("Peu ou pas de coupe", 
               "Coupe modérée", 
               "Coupe majeure")
    } else {
      lgd <- c("Minor logging", 
               "Moderate logging", 
               "Major logging")
    }
  }

  
  if(is.null(ylim)) ylim <- range(index)
  
  bp <- barplot(as.matrix(index), beside = T, plot = F)
  
  plot0(xlim = range(bp), ylim = ylim, frame.plot = TRUE, yaxs = "i")
  axis(2, las = 1, cex.axis = 1.2)
  mtext(states, 1, at = bp[1,], line = 1, xpd = NA, adj = 0, cex = 1.9, font = 2)
  mtext(ylab, 2, line = 3, font = 2, cex = 1.8)
  
  lines(index[bars,] ~ bp[bars,], type = "h", lwd = 20, lend = 2, col = colss[bars])
  
  legend("topleft", legend = lgd[bars], cex = 1.8,
         pt.cex = 0, text.col = colss[bars], text.font = 2, 
         bty = "n")
}

barplot_halflife <- function(index = NULL, ylim = NULL, bars = 1:3,
                          colss = c("#f1ba53","#E38451", "#b5305d"),
                          lang = "fr") {
  par(mar = c(4,5,.7,.5))
  
  if(lang == "fr") {
    lgd <- c("Peu ou pas de coupe", 
             "Coupe modérée", 
             "Coupe majeure")
    ylab = "Temps de convergence (années)"
  } else {
    lgd <- c("Minor logging", 
             "Moderate logging", 
             "Major logging")
    ylab = "Convergence time (years)"
  }
  
  if(is.null(ylim)) ylim <- range(index)
  
  bp <- barplot(index, plot = FALSE)
  
  plot0(xlim = c(0.3, 3.6), ylim = ylim, frame.plot = TRUE, yaxs = "i")
  
  axis(2, las = 1, cex.axis = 1.2)

  mtext(ylab, 2, line = 3, font = 2, cex = 1.8)
  
  lines(index[bars] ~ bp[bars], type = "h", lwd = 20, lend = 2, col = colss[bars])

}
