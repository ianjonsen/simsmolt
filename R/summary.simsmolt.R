##' @importFrom dplyr %>% summarise group_by n
##' @method summary simsmolt
##' @export
summary.simsmolt <- function(x, data = NULL, ...) {
  if (length(list(...)) > 0) {
    warning("additional arguments ignored")
  }
  
  if(length(class(x)) == 1) {
    stop("summary not yet implemented for single simulation runs...\n")
  }
  
  else if(class(x)[2] == "rowwise_df" || class(x)[2] == "grouped_df") {
    compl <- sapply(x$rep, function(.) !inherits(., "try-error"))
    cat(sprintf("dropping %i failed runs\n\n", sum(!compl)))
    x <- x[compl, ]
  
    ## count smolts passing each receiver line/array (whether or not detected)
    r <- data$recs
    yrec <- unique(r$y)
    maxy <- sapply(x$rep, function(.) nrow(.$sim)) %>% max()
    na.fill <- function(., maxy) {
      y <- .$sim$y
      ly <- length(y)
      if(ly < maxy) {
        y <- c(y, rep(NA, maxy - ly))
      }
      y
    }
    sm.y <- lapply(x$rep, na.fill, maxy) %>% do.call(cbind, .)
    
    fn <- function(x) sapply(1:length(yrec), function(i) sum(sum((yrec[i] - x) < 0, na.rm = TRUE) > 0))
    num.cross <- apply(sm.y, 2, fn) %>% apply(., 1, sum)

    ## summarise multiple replicates
    if(names(x)[2] != "rep") stop("expecting simulation output objects to be named 'rep'")
    all.tr <- lapply(x$rep, function(.) .$trans) %>% do.call(rbind, .)
    all.dt <- lapply(x$rep, function(.) .$detect) %>% do.call(rbind, .)
    dt.by.line <- all.dt %>% group_by(line) %>% summarise(n = n())
    tr.by.line <- all.tr %>% group_by(line) %>% summarise(n = n())
    dt.by.recline <- all.dt %>% group_by(line,recv_id) %>% summarise(n = n()) %>% group_by(line) %>% summarise(n = n())
  
    ## count individual smolts detected at each line
    num.ind.smolt.dt.by.line <- all.dt %>% group_by(line, trns_id) %>% summarise(n = n())
    num.smolt.dt.by.line <-  num.ind.smolt.dt.by.line %>% summarise(n = n())
    
    lines <- paste0("l",1:length(yrec))
    if(nrow(dt.by.line != length(yrec))) {
      ex <- which(!lines %in% dt.by.line$line)
      dt.by.line <- rbind(dt.by.line, data.frame(line = lines[ex], n = ex*0)) %>%
        arrange(line)
    }
    if(nrow(tr.by.line != length(yrec))) {
      ex <- which(!lines %in% tr.by.line$line)
      tr.by.line <- rbind(tr.by.line, data.frame(line = lines[ex], n = ex*0)) %>%
        arrange(line)
    }
    if(nrow(dt.by.recline != length(yrec))) {
      ex <- which(!lines %in% dt.by.recline$line)
      dt.by.recline <- rbind(dt.by.recline, data.frame(line = lines[ex], n = ex*0)) %>%
        arrange(line)
    }
    if(nrow(num.smolt.dt.by.line) != length(yrec)) {
      ex <- which(!lines %in% num.smolt.dt.by.line$line)
      num.smolt.dt.by.line <- rbind(num.smolt.dt.by.line, data.frame(line = lines[ex], n = ex*0)) %>%
        arrange(line)
    }
  
  
  n <- c(num.cross, nrow(x))
  ndt <- c(num.smolt.dt.by.line$n, length(unique(num.ind.smolt.dt.by.line$trns_id)))

  structure(list(
    n = n,
    ndt = ndt,
    dt.num = c(dt.by.line$n, nrow(all.dt)),
    dt.by.recline = c(dt.by.recline$n, sum(dt.by.recline$n)),
    p.dt = c(dt.by.line$n, nrow(all.dt)) / c(tr.by.line$n, nrow(all.tr)),
    p.smolt = ndt/n
  ),
  class = "summary.simsmolt")
  }
  
}

##' @method print summary.simsmolt
##' @export
print.summary.simsmolt <- function(x)
{

  xx <- t(cbind(c(NA, x$n), 
                c(NA, x$ndt),
                   c(NA, x$dt.num), 
                   c(NA, x$dt.by.recline),
                   c(NA, round(x$p.dt, 4)),
                   c(NA, round(x$p.smolt, 4))
                   ))
  
  xx <- cbind(c("smolts", "smolts detected", "detections", "receivers with detections", "p(transmissions detected)", "p(smolts detected)"),
              xx)

  dimnames(xx) <- list(rep("", 6), 
                       c("", " ", paste0("line.", 1:(length(x$n) - 1)), "total")
                       )
  
  print(xx, justify="right", right=TRUE, quote=FALSE, na.print="")

 # invisible(x)
}## print.summary.simsmolt
