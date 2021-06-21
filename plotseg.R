

plotseg = function( x1, x2, object, label1 = "x1", label2 = "x2", Pvalue = NULL, 
                    title = "", AddLegend = T, plotOtherDirect = FALSE, 
                    gg = T, plotfit = T, plotTwoBP = T, pt_size = 1.2){
  
  ### object fitted model from function findDirec()
  
  library(ggplot2)
  
  coefs.p = object$coef.p
  coefs.c = object$coef.c
  
  if( object$parent == 1 ) {
    np = x1; nc = x2
    label_p = label1; label_c = label2; 
  } else {
    np = x2; nc = x1
    label_p = label2; label_c = label1; 
  }
  
  dat = as.data.frame(cbind(np, nc))
  # colnames(dat) = c(label_p, label_c)
  if ( plotOtherDirect == FALSE ) {
    
    if(gg){
      if(plotfit){
        
        p = ggplot(data = dat, aes(x=np, y=nc)) + 
          geom_point(color = "grey45", size = pt_size) + 
          geom_vline(aes(xintercept = object$tau.p), #, color = paste0(label_p, ": tau =", round(object$tau.p, 3))
                     color = "coral",size =1.5, linetype = "dashed") +
          # geom_hline(aes(yintercept = object$tau.c), #, color = paste0(label_c, ": tau =", round(object$tau.c, 3))
          #            color = "turquoise", size =1.5, linetype = "dashed") + #col = "skyblue2", lty = 2, lwd = 1
          geom_segment(aes(x = min(np), y = coefs.p[1] + coefs.p[2]*min(np), 
                           xend = object$tau.p, yend = coefs.p[1] + coefs.p[2]*object$tau.p), 
                       color = "coral") +
          geom_segment(aes(x = object$tau.p, y = coefs.p[3] + coefs.p[4]*object$tau.p, 
                           xend = max(np), yend = coefs.p[3] + coefs.p[4]*max(np)), 
                       color = "coral") +  
          xlab(label_p) + ylab(label_c) +
          theme_bw() +
          theme(panel.grid = element_blank(), panel.background = element_blank(), 
                axis.text=element_text(size=10),
                axis.title=element_text(size=14, colour = "grey25"))          
        
        if(plotTwoBP){
          p = p + geom_hline(aes(yintercept = object$tau.c), #, color = paste0(label_c, ": tau =", round(object$tau.c, 3))
                         color = "turquoise", size =1.5, linetype = "dashed") #col = "skyblue2", lty = 2, lwd = 1
        }
        # p = ggplot(data = dat, aes(x=np, y=nc)) + 
        #   geom_point(color = "grey45", size = pt_size) + 
        #   geom_vline(aes(xintercept = object$tau.p), #, color = paste0(label_p, ": tau =", round(object$tau.p, 3))
        #              color = "coral",size =1.5, linetype = "dashed") +
        #   geom_hline(aes(yintercept = object$tau.c), #, color = paste0(label_c, ": tau =", round(object$tau.c, 3))
        #              color = "turquoise", size =1.5, linetype = "dashed") + #col = "skyblue2", lty = 2, lwd = 1
        #   geom_segment(aes(x = min(np), y = coefs.p[1] + coefs.p[2]*min(np), 
        #                    xend = object$tau.p, yend = coefs.p[1] + coefs.p[2]*object$tau.p), 
        #                color = "coral") +
        #   geom_segment(aes(x = object$tau.p, y = coefs.p[3] + coefs.p[4]*object$tau.p, 
        #                    xend = max(np), yend = coefs.p[3] + coefs.p[4]*max(np)), 
        #                color = "coral") +  
        #   xlab(label_p) + ylab(label_c) +
        #   theme_bw() +
        #   theme(panel.grid = element_blank(), panel.background = element_blank(), 
        #         axis.text=element_text(size=10),
        #         axis.title=element_text(size=14, colour = "grey25"))
      }else{
        p = ggplot(data = dat, aes(x=np, y=nc)) + 
          geom_point(color = "grey45", size = pt_size) + 
          # geom_vline(aes(xintercept = object$tau.p), #, color = paste0(label_p, ": tau =", round(object$tau.p, 3))
          #            color = "coral",size =1.5, linetype = "dashed") +
          # geom_hline(aes(yintercept = object$tau.c), #, color = paste0(label_c, ": tau =", round(object$tau.c, 3))
          #            color = "turquoise", size =1.5, linetype = "dashed") + #col = "skyblue2", lty = 2, lwd = 1
          xlab(label_p) + ylab(label_c) +
          theme_bw() +
          theme(panel.grid = element_blank(), panel.background = element_blank(), 
                axis.text=element_text(size=10),
                axis.title=element_text(size=14, colour = "grey25"))
      }
      
      return(p)
    }else{
      plot(np, nc, cex = 0.1, xlab = label_p, ylab = label_c, main = title)
      abline(v = object$tau.p, col = "red")
      abline(h = object$tau.c, col = "green")
      segments(min(np), coefs.p[1] + coefs.p[2]*min(np),
               object$tau.p, coefs.p[1] + coefs.p[2]*object$tau.p,
               col = "blue", lty = 2, lwd = 2.5)
      segments(object$tau.p, coefs.p[3] + coefs.p[4]*object$tau.p,
               max(np), coefs.p[3] + coefs.p[4]*max(np),
               col = "blue", lty = 3, lwd = 2.5)
      
      if(AddLegend) {
        
        
        
        legend("topleft", col = c("red", "green"), lty = c(1,1), bty = 'n', cex = .5,
               legend = c(paste(label_p, ":", "tau", "=", round(object$tau.p, 3)), 
                          paste(label_c, ":", "tau", "=", round(object$tau.c, 3))))
        
        if(is.null(Pvalue)) {
          legend("topright", bty = 'n', cex = .5, 
                 legend = c(paste("parent", ":", label_p),
                            paste("testStat", ":", round(object$testStat, 3))))
        }else {
          legend("topright", bty = 'n', cex = .5, 
                 legend =  c(paste("parent", ":", label_p),
                             paste("testStat", "=", round(object$testStat, 3)),
                             paste("Pvalue", "=", Pvalue))) 
        }
      }
    }
      # annotate("text", x = quantile(np, 0.002), y = max(nc), size = 3, color = "grey25",
      #          label = paste0("parent: ", label_p, "\n", "RR: ", round(object$testStat, 3), "\n")) +
      # theme(legend.title = element_blank())
    
    } else {
      
    if(gg) {
      p = ggplot(data = as.data.frame(cbind(np,nc)), aes(x=nc, y=np)) + 
        geom_point(col = "grey45", size = pt_size) + 
        geom_vline(aes(xintercept = object$tau.c, color = label_c), 
                   color = "turquoise", size =1.5, linetype = "dashed") +
        # geom_hline(aes(yintercept = object$tau.p, color = label_p), 
        #            color = "coral", size =1.5, linetype = "dashed") +
        geom_segment(aes(x = min(nc), y = coefs.c[1] + coefs.c[2]*min(nc), 
                         xend = object$tau.c, yend = coefs.c[1] + coefs.c[2]*object$tau.c), 
                     color = "turquoise") +
        geom_segment(aes(x = object$tau.c, y = coefs.c[3] + coefs.c[4]*object$tau.c, 
                         xend = max(nc), yend = coefs.c[3] + coefs.c[4]*max(nc)), 
                     color = "turquoise") +
        xlab(label_c) + ylab(label_p) +
        theme_bw() +
        theme(panel.grid = element_blank(), panel.background = element_blank(), 
              axis.text=element_text(size=10),
              axis.title=element_text(size=14, colour = "grey25"))
      
      if(plotTwoBP){
        p = p + geom_hline(aes(yintercept = object$tau.p, color = label_p),
                       color = "coral", size =1.5, linetype = "dashed")
      }
      
      # p = ggplot(data = as.data.frame(cbind(np,nc)), aes(x=nc, y=np)) + 
      #   geom_point(col = "grey45", size = pt_size) + 
      #   geom_vline(aes(xintercept = object$tau.c, color = label_c), 
      #              color = "turquoise", size =1.5, linetype = "dashed") +
      #   geom_hline(aes(yintercept = object$tau.p, color = label_p), 
      #              color = "coral", size =1.5, linetype = "dashed") +
      #   geom_segment(aes(x = min(nc), y = coefs.c[1] + coefs.c[2]*min(nc), 
      #                    xend = object$tau.c, yend = coefs.c[1] + coefs.c[2]*object$tau.c), 
      #                color = "turquoise") +
      #   geom_segment(aes(x = object$tau.c, y = coefs.c[3] + coefs.c[4]*object$tau.c, 
      #                    xend = max(nc), yend = coefs.c[3] + coefs.c[4]*max(nc)), 
      #                color = "turquoise") +
      #   xlab(label_c) + ylab(label_p) +
      #   theme_bw() +
      #   theme(panel.grid = element_blank(), panel.background = element_blank(), 
      #         axis.text=element_text(size=10),
      #         axis.title=element_text(size=14, colour = "grey25"))
      
      return(p)
    }else{
      plot(nc, np, cex = .1, xlab = label_c, ylab = label_p, main = title)
      abline(v = object$tau.c, col = "green")
      abline(h = object$tau.p, col = "red")
      segments(min(nc), coefs.c[1] + coefs.c[2]*min(nc),
               object$tau.c, coefs.c[1] + coefs.c[2]*object$tau.c,
               col = "blue", lty = 2, lwd = 2.5)
      segments(object$tau.c, coefs.c[3] + coefs.c[4]*object$tau.c,
               max(nc), coefs.c[3] + coefs.c[4]*max(nc),
               col = "blue", lty = 3, lwd = 2.5)
      
      if(AddLegend) {
        
        
        
        legend("topleft", col = c("red", "green"), lty = c(1,1), bty = 'n', cex = .5,
               legend = c(paste(label_p, ":", "tau", "=", round(object$tau.p, 3)), 
                          paste(label_c, ":", "tau", "=", round(object$tau.c, 3))))
        
        if(is.null(Pvalue)) {
          legend("topright", bty = 'n', cex = .5, 
                 legend = c(paste("parent", ":", label_p),
                            paste("testStat", ":", round(object$testStat, 3))))
        }else {
          legend("topright", bty = 'n', cex = .5, 
                 legend =  c(paste("parent", ":", label_p),
                             paste("testStat", "=", round(object$testStat, 3)),
                             paste("Pvalue", "=", Pvalue))) 
        }
      }
    }
  }
  

}




