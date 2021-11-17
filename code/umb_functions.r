
sp_short <- function(x) {
  sp_short_mid <- function(x) {
    if (length(grep("_",x))>0 | length(grep(" ",x))>0) {
      x <- sapply(strsplit(x, "[_ ]"),function(x)paste0(substr(x[1],1,1),". ",paste(x[-1], collapse="_")))
    }
    return(x)
  }
  return(as.character(sapply(x,function(x)ifelse(length(grep("noname",x))>0,x,sp_short_mid(x)))))
}

breakup <- function(string, sep, elements=NULL, 
                    fromend=NULL, join=NULL, seq=FALSE) {
  if (seq) {
    ll <- lapply(strsplit(string, sep[1]), function(x){x[elements[1]]})
  } else {
    ll <- strsplit(string, sep[1])    
  }
  if (length(sep)>1) {
    for (i in 2:length(sep)) {
      if (seq) {
        ll <- lapply(ll, function(x){unlist(strsplit(x, sep[i]))[elements[i]]})
      } else {
        ll <- lapply(ll, function(x){unlist(strsplit(x, sep[i]))})
      }
    }
  }
  if (seq) {
    return(sapply(ll, function(x){paste(x,collapse=join)}))
  } else {
    return(sapply(ll, function(x){paste(x[c(elements, length(x)-fromend+1)],collapse=join)}))
  }
}

boxer <- function(data, x, width=1, quants=c(0.05,0.25,0.5,0.75,0.95), border="black", 
                  col=NA, flicks=0, lwd=1, horiz=FALSE, ...) {
  #qtl <- quantile(data,quants)
  if (!is.list(data)) {
    data <- list(data)
  }
  if (length(x)!=length(data)) {
    stop("'data' must be same length as 'x'")
  }
  if (length(quants)!=5) {
    stop("'quants' must be vector of length 5")
  }
  qtl <- as.matrix(sapply(data,quantile, quants))
  if (!horiz) {
    rect(x-width/2,qtl[2,],x+width/2,qtl[4,], col=col, border=border, lwd=lwd, ...)
    segments(rep(x,2),c(qtl[1,],qtl[4,]),rep(x,2),c(qtl[2,],qtl[5,]), col=border, lwd=lwd, ...)
    segments(x-width/2, qtl[3,], x+width/2, qtl[3,], col=border, lwd=lwd, ...)
    if (flicks>0) {
      segments(rep(x-flicks*width/2,2),c(qtl[1,],qtl[5,]),rep(x+flicks*width/2,2),c(qtl[1,],qtl[5,]), col=border, lwd=lwd, ...)
    }
  } else {
    rect(qtl[2,],x-width/2,qtl[4,],x+width/2, col=col, border=border, lwd=lwd, ...)
    segments(c(qtl[1,],qtl[4,]),rep(x,2),c(qtl[2,],qtl[5,]),rep(x,2), col=border, lwd=lwd, ...)
    segments(qtl[3,],x-width/2, qtl[3,], x+width/2, col=border, lwd=lwd, ...)
    if (flicks>0) {
      segments(c(qtl[1,],qtl[5,]),rep(x-flicks*width/2,2),c(qtl[1,],qtl[5,]),rep(x+flicks*width/2,2), col=border, lwd=lwd, ...)
    }
  }
}

