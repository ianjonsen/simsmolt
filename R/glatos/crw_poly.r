crwr <- function(rstr, theta=c(0,10), stepLen=1, initPos=c(NA,NA), 
                           initHeading=45, nsteps=30){          
  source("fn/crw.r")
  source("fn/vector_heading.r")
  source("fn/rotate_points.r")
  
  ## define SoBI as an attracting point for salmon approaching & 
  ##  a replusion point after arrival
  Z <- cbind(740, 785)
  ## Weibull parameters for random displacements rather than fixed stepLen
  a <- 3
  b <- 5
  
  #randomly select one point in the study area
  if(any(is.na(initPos))){
    inRstr <- FALSE #logical flag; preallocate
    ex <- raster::extent(rstr)
    while(!inRstr){
      init <- c(runif(1, ex[1],ex[2]),
                runif(1, ex[3],ex[4]))
      inRstr <- raster::extract(rstr, cbind(init[1],init[2]), method="simple") %>%
        as.logical()
    } #end while
  } #end if

  if(all(!is.na(initPos))) {
    init <- initPos
    if(raster::extract(rstr,cbind(init[1],init[2])) == 0) {
      stop("initPos is outside viable region")
    }
  } 
  
  #randomly select heading
  if(is.na(initHeading)) initHeading <- runif(1,0,360)
  
  winSz <- 10 #size of window
  
  mod <- nsteps%/%winSz #number of times to loop through full window
  rmd <- nsteps%%winSz #size of last window
  
  path.fwd <- data.frame(x=rep(NA,nsteps+1),y=NA) #preallocate
  path.fwd[1,] <- init
  rows.i <- 1 #rows for ith window iteration; preallocate
  
  #initialize progress bar
  pb <- txtProgressBar(min=0, max=mod+1, initial=0, style=3)	  
  
  nwin <- mod+(rmd>0) #number of windows
  for(i in 1:nwin){
    #simulate track forward
    #update starting point and heading
    init <- as.vector(path.fwd[max(rows.i),]) #start at previous end
    if(i>1) {
      initHeading <- vector_heading(path.fwd$x[max(rows.i)-1:0],
                                    path.fwd$y[max(rows.i)-1:0])
    }
    if(i==nwin & rmd>0) winSz <- rmd #size of last window
    rows.i <- max(rows.i)+(1:winSz) #update rows for this iteration
    path.fwd[rows.i,] <- crw(theta=theta,stepLen=stepLen,initPos=init,
                             initHeading,nsteps=winSz)

    #check if in polygon
    inRstr <- raster::extract(rstr, cbind(path.fwd$x[rows.i],path.fwd$y[rows.i]), method = "simple") %>% 
      as.logical()
    
    while(!all(inRstr)){
      outside <- match(0,inRstr) #identify first step outside
      #truncate turn angle distribution based on encounter with boundary
      
      constrainThetaToRaster <- function(x,y,len,rstr,theta){
        heading_deg <- vector_heading(x,y)
        heading_deg2 <- ((heading_deg-180):(heading_deg+180)) %% 360
        heading_rad <- heading_deg2 * pi/180 #convert to radians
        xlen <- sin(heading_rad)*len #x-component vector
        ylen <- cos(heading_rad)*len #y-component vector
        bar <- data.frame(x=x[2]+xlen,y=y[2]+ylen,theta=-180:180)
        #check if sampled points in polygon
        bar$inRstr <- raster::extract(rstr, cbind(bar$x,bar$y), method="simple")
        if(sum(bar$inRstr) == 0) stop("Path stuck at a boundary.")
        probs <- bar$inRstr * dnorm(-180:180,theta[1],theta[2])
        theta2 <- sample(bar$theta, 1, prob=probs)
        return(theta2) #new turn angle			
      } 
     
      theta2 <- constrainThetaToRaster(x=path.fwd$x[rows.i[outside]+(-2:-1)],
                                        y=path.fwd$y[rows.i[outside]+(-2:-1)],stepLen,rstr,theta)	  
      
      #rotate track from first step outside polygon to end
      path.fwd[rows.i,][outside:length(rows.i),] <- rotate_points(
        x=path.fwd$x[rows.i][outside:length(rows.i)],
        y=path.fwd$y[rows.i][outside:length(rows.i)],
        theta=theta2, focus=path.fwd[rows.i[outside]-1,])
      
      #identify points in and out of polygon	
      inRstr <- raster::extract(rstr, cbind(path.fwd$x[rows.i],path.fwd$y[rows.i]), method = "simple") %>%
        as.logical()
      
    } #end while
    
    #update progress bar
    setTxtProgressBar(pb, i)
    if(i==(nwin)) close(pb)		
  } #end i
  
  return(path.fwd)
} #end crwInRect def
