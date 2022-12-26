plotGradient_modify <- function (hM, Gradient, predY, measure, xlabel = NULL, 
                                 ylabel = NULL, index = 2, q = c(0.025, 0.5, 0.975),
                                 showData = TRUE, yshow = NA, jigger = 0) 
{

  Pr = NA
  if (is.null(xlabel)) {
    switch(class(hM$X)[1L], matrix = {
      xlabel = colnames(Gradient$XDataNew)[[1]]
    }, list = {
      xlabel = colnames(Gradient$XDataNew[[1]])[[1]]
    })
  }
  switch(class(hM$X)[1L], matrix = {
    xx = Gradient$XDataNew[, 1]
  }, list = {
    if (measure == "Y") {
      xx = Gradient$XDataNew[[index]][, 1]
    } else {
      xx = Gradient$XDataNew[[1]][, 1]
    }
  })
  ngrid = length(xx)
  if (measure == "S") {
    predS = abind::abind(lapply(predY, rowSums), along = 2)
    Pr = mean(predS[ngrid, ] > predS[1, ])
    qpred = apply(predS, c(1), quantile, probs = q, na.rm = TRUE)
    if (is.null(ylabel)) {
      ylabel = "Summed response"
      if (all(hM$distr[, 1] == 2)) {
        ylabel = "Species richness"
      }
      if (all(hM$distr[, 1] == 3)) {
        ylabel = "Total count"
      }
    }
  }
  if (measure == "Y") {
    tmp = abind::abind(predY, along = 3)
    Pr = mean(tmp[ngrid, index, ] > tmp[1, index, ])
    qpred = apply(tmp, c(1, 2), quantile, probs = q, na.rm = TRUE)
    qpred = qpred[, , index]
    if (is.null(ylabel)) {
      ylabel = hM$spNames[[index]]
    }
  }
  if (measure == "T") {
    if (all(hM$distr[, 1] == 1)) {
      predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)), 
                                                                      hM$nt), ncol = hM$nt))
    }
    else {
      predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a), 
                                                                 hM$nt), ncol = hM$nt))
    }
    predT = abind::abind(predT, along = 3)
    Pr = mean(predT[ngrid, index, ] > predT[1, index, ])
    qpred = apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)
    qpred = qpred[, , index]
    if (is.null(ylabel)) {
      ylabel = hM$trNames[[index]]
    }
  }
  lo = qpred[1, ]
  hi = qpred[3, ]
  me = qpred[2, ]
  lo1 = min(lo, yshow, na.rm = TRUE)
  hi1 = max(hi, yshow, na.rm = TRUE)
  if (showData) {
    switch(class(hM$X)[1L], matrix = {
      XDatacol = which(colnames(Gradient$XDataNew)[[1]] == 
                         colnames(hM$XData))
    }, list = {
      XDatacol = which(colnames(Gradient$XDataNew[[1]])[[1]] == 
                         colnames(hM$XData[[1]]))
    })
    if (measure == "S") {
      pY = rowSums(hM$Y, na.rm = TRUE)
    }
    if (measure == "Y") {
      pY = hM$Y[, index]
    }
    if (measure == "T") {
      if (all(hM$distr[, 1] == 1)) {
        tmp = (exp(hM$Y) %*% hM$Tr)/matrix(rep(rowSums(exp(hM$Y)), 
                                               hM$nt), ncol = hM$nt)
      }
      else {
        tmp = (hM$Y %*% hM$Tr)/matrix(rep(rowSums(hM$Y), 
                                          hM$nt), ncol = hM$nt)
      }
      pY = tmp[, index]
    }
    if (!is.numeric(pY)) 
      pY <- as.numeric(pY)
    switch(class(hM$X)[1L], matrix = {
      pX = hM$XData[, XDatacol]
    }, list = {
      pX = hM$XData[[1]][, XDatacol]
    })
    hi1 <- max(hi1, max(pY, na.rm = TRUE))
    lo1 <- min(lo1, min(pY, na.rm = TRUE))
  }
  if (is.factor(xx)) {
    toPlot = data.frame(xx, me, lo, hi, stringsAsFactors = TRUE)
    pl = ggplot(toPlot, aes_string(x = xx, y = me)) + geom_bar(position = position_dodge(), 
                                                               stat = "identity") + xlab(xlabel) + ylab(ylabel) + 
      geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, 
                    position = position_dodge(0.9))
    if (showData) {
      if (jigger > 0) {
        pX = as.numeric(pX)
        pX = pX + runif(n = length(pY), min = -jigger, 
                        max = jigger)
      }
      dataToPlot = data.frame(pX = pX, pY = pY, stringsAsFactors = TRUE)
      pl = pl + geom_point(data = dataToPlot, aes_string(x = pX, 
                                                         y = pY), size = 1)
    }
    return(pl)
  } else {
    if (inherits(hM$X, "list") && !measure == "Y") {
      plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", 
           xaxt = "n", xlab = xlabel, ylab = ylabel)
      axis(1, c(min(xx), (min(xx) + max(xx))/2, max(xx)), 
           c("min", "mean", "max"))
    }
    else {
      df <- cbind(xx, v1 = qpred[2,])
      dataToPlot <- data.frame(pX = pX, pY = pY)
      ct <- cor.test(df[, "xx"], df[, "v1"])
      
      pred.response <- ggplot(as.data.frame(df), aes(x = xx, y = v1)) +
        ylim(lo1, hi1) + geom_point(data = dataToPlot, aes(x = pX, pY), color = "gray60") +
        labs(caption = bquote(r==.(round(ct$estimate, 3))),
      x = gsub("_", " ", xlabel), y = ylabel) +  
        geom_smooth(method = "lm", size = .7, color = "black", fill = "#AA3377") +
        theme(legend.position = "none") + theme_classic() 
      }
return(pred.response)
  }
}
