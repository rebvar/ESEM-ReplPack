options(check.names = FALSE)
library(ScottKnott)
library(ScottKnottESD)


Partition <- function(g,
                      means,
                      mMSE,
                      dfr,
                      sig.level,
                      av,
                      k,
                      group,
                      ngroup,
                      markg,
                      g1=g,
                      sqsum=rep(0, g1))
{
  for(k1 in k:(g-1)) {
    t1 <- sum(means[k:k1])
    
    k2 <- g-k1
    
    t2 <- sum(means[(k1+1):g])
    
    # SK between groups sum of squares
    sqsum[k1] <- t1^2/(k1-k+1) + t2^2/k2 - (t1+t2)^2/(g-k+1)
  }
  
  # the first element of this vector is the value of k1 which maximizes sqsum (SKBSQS)
  ord1 <- order(sqsum, decreasing=TRUE)[1]
  
  ############################################################
  # Compute the magnitude of the difference of all treatments in a group
  # Use effect size test to identify whether to split metrics into two groups of {k:i} and {j:g}.
  # The split is accepted if effect.size({k}, {g}) != negligible
  diff <- function(k, g, av, means) {
    if(k==g){ # if k and g are the same treatment
      return(TRUE)
    }
    
    a <- av$model[av$model[,2] == names(means[k]),1]
    b <- av$model[av$model[,2] == names(means[g]),1]  
    
    magnitude <- as.character(cohen.d(a,b)$magnitude)
    return(magnitude == "negligible")
  }
  ############################################################
  
  # if true it returns one node to the right if false it goes forward one node to the left
  if(diff(k, g, av, means) | 
     (ord1 == k)) {
    # In the case of a single average left (maximum)
    if(!diff(k, g, av, means)) {
      # it marks the group to the left consisting of a single mean
      ngroup <- ngroup + 1
      
      # it forms a group of just one mean (the maximum of the group)
      group[k] <- ngroup
      
      # lower limit on returning to the right
      k <- ord1 + 1
    }
    if(diff(k, g, av, means)) {
      # it marks the groups
      ngroup <- ngroup + 1
      
      # it forms a group of means
      group[k:g] <- ngroup
      
      # if this group is the last one
      if (prod(group) > 0)
        # If the upper limit of the latter group formed is equal to the total
        # number of treatments than  the grouping algorithm is ended
        return (group)
      
      # it marks the lower limit of the group of means to be used in the
      # calculation of the maximum sqsum on returning one node to the right
      k <- g + 1
      
      # it marks the upper limit of the group of means to be used in the
      # calculation of the maximum sqsum on returning one node to the right
      g <- markg[g]
    }
    while(k == g) {
      # there was just one mean left to the right, so it becomes a group
      ngroup   <- ngroup + 1
      
      group[g] <- ngroup
      
      if(prod(group) > 0)
        # If the upper limit of the latter group formed is equal to the total
        # number of treatments than  the grouping algorithm is ended
        return(group)
      
      # the group of just one mean group had already been formed, a further
      # jump to the right and another check whether there was just one mean
      # left to the right
      k <- g + 1
      
      g <- markg[g]
    }
  } else {
    # it marks the upper limit of the group split into two to be used on
    # returning to the right later
    markg[ord1] <- g
    
    g <- ord1
  }
  
  Partition(g,
            means,
            mMSE,
            dfr,
            sig.level,
            av,
            k,
            group,
            ngroup,
            markg)
}


scottknott <- function(x,
                       which=NULL,
                       id.trim=3,
                       sig.level=.05,
                       dispersion=c('mm', 's', 'se'), ...)   
{
  if(is.null(which))
    which <- names(x$model)[2]
  
  mt <- model.tables(x, "means")  # summary tables for model fits
  
  if(is.null(mt$n))
    stop("No factors in the fitted model!")
  
  r   <- mt$n[names(mt$tables)][[which]] # groups and its number of replicates
  
  MSE <- deviance(x)/df.residual(x)
  
  m   <- as.vector(mt$tables[[which]])   # means
  
  nms <- names(mt$tables[[which]])
  
  ord <- order(m, decreasing=TRUE)
  
  m.inf <- m.inf.1a(x,
                    which,
                    dispersion)
  
  rownames(m.inf) <- nms  
  
  m.inf <- m.inf[order(m.inf[,1],
                       decreasing=TRUE),]
  
  mMSE <- MSE / r
  
  dfr  <- x$df.residual  # residual degrees of freedom
  
  g    <- nrow(m.inf)
  
  groups <- Partition(g,
                      m.inf[, 1],
                      mMSE,
                      dfr,
                      sig.level=sig.level,
                      av=x,
                      1,
                      rep(0, g),
                      0,
                      rep(0, g))
  
  res <- list(av=x,
              groups=groups,
              nms=nms,
              ord=ord,
              m.inf=m.inf,
              sig.level=sig.level)
  
  class(res) <- c('sk_esd',
                  'list')
  
  invisible(res)
}


sk_esd2 <- function(x, alpha=0.05, ...){
  x <- data.frame(x,check.names = FALSE)
  av <- aov(value ~ variable, data=reshape2::melt(x, id.vars=0 )) 
  sk <- scottknott(av, which='variable',  dispersion='s', sig.level=alpha) 
  names(sk$groups) <- rownames(sk$m.inf)
  class(sk) <- c(class(sk),"sk_esd")
  return(sk)
}


mpath = 'OutF/'
savepath = "GC3/"
setwd(mpath)

files = list.files(path = ".")


framesLst = list()


i = 1
for(file in files)
{
  
  framesLst[[i]] <- as.data.frame(read.csv(paste(mpath,file,sep=''), check.names=FALSE), check.names=FALSE)
  i=i+1
}

mar.default <- c(9,4,0.8,0.5)

countFiles = length(files)


fileCnt = countFiles/4
mCnt = 4

mes = c('Accuracy','F-measure','Precicion','Recall')
mtext = ""
for (i in 1:mCnt)
{
  df = data.frame() 
  
  
  for (j in 1:fileCnt)
  {
    findex = (i-1)*fileCnt+j
    print (findex)
    df = rbind(df, framesLst[[findex]])
  }
  names(df) = names(framesLst[[1]])
  #setEPS()
  #postscript(paste(savepath,"0-",mes[i],".eps", sep=""), width=400, height = 150)
  png(paste(savepath,"0-",mes[i],".png", sep=""), width=2000, height = 1100, res=200) 
  par(mar = mar.default)
  
  sk <- sk_esd2(df, title = "")
  
  plot(sk,las=2, main="",ylim=c(0, 1), xlab = "", ylab=mes[i], sub="", title = "")
  
  #rect(5, -0.4, 15, 0.2, border = TRUE)
  
  dev.off()
  graphics.off()
  
  
  pdf(paste(savepath,"0-",mes[i],".pdf", sep=""), width=2000, height = 1100)
  par(mar = mar.default)
  
  sk <- sk_esd(df, title="")
  
  plot(sk,las=2, main="",ylim=c(0, 1), xlab = "", title="")
  
  
  title(ylab=mes[i], line=0, cex.lab=0.1)
  dev.off()
  graphics.off()
  
  
}

