"plotltr" <-
function(x, which="dist", pch = 16, cex = 0.7, addlines = TRUE, ...)
  {
    if (!inherits(x,"ltraj"))
        stop("x should be of class ltraj")
    opar <- par(mfrow=n2mfrow(length(x)))
    on.exit(par(opar))
    toto <- lapply(x, function(i) {
        if (!is.null(attr(i, "infolocs")))
            i <- cbind(i, attr(i, "infolocs"))
        ex<- parse(text=which)
        coin <- eval(ex, envir=i)
        plot(i$date, coin, main=attr(i,"burst"), xlab="Time",
             ylab=which, pch=pch, cex=cex,...)
        if (addlines)
            lines(i$date, coin)
    })
    invisible()
  }

