# Functions ---------------------------------------------------------------

# 1
viable <- function(x){
  Viable <- vector()
  Missings <- vector()
  for (i in 1:ncol(x)) {
    Missings[i] <- sum(is.na(x[i]))
    Viable[i] <- nrow(x) - Missings[i]
  }
  data.frame(Var=colnames(x), Viable = Viable, Missings = Missings)
}


# 2
cv <- function(ave, std){
  std/ave * 100; 
  #  print(cv)
}


# 3
is.outlier = function (x) {
  # See: Davies, P.L. and Gather, U. (1993).
  # "The identification of multiple outliers" (with discussion)
  # J. Amer. Statist. Assoc., 88, 782-801.
  
  x <- na.omit(x)
  lims <- median(x) + c(-1, 1) * 2.5 * mad(x)
  x < lims[1] | x > lims[2]
}


# 4
CIwidth <- function(mean, sd, N){
  se <- sd / sqrt(N);
  upper <- mean + 1.96 * se # Get upper limit for 95% CI, with z=1.96. Or a 99% one with 2.575
  lower <- mean - 1.96 * se # Get lower limit for 95% CI.
  #width <- upper - lower
  upper - lower
  #print(width)
}


# 5
## STEP 1
getStability <- function(perc=perc, dat=dat, iter=iter, method=method){
  require(mclust)
  require(clValid)
  
  classif2 <- matrix(ncol=3, nrow=iter)
  classif <- list(1,2,3,4,5)
  for (i in 1:5){
    classif[[i]] <- vector("list", iter)
  }

  if (method=="kmeans"){
    for (j in 1:iter){
      subsamp <- dat[sample(row.names(dat), round((perc*nrow(dat))/100, 0), replace=F), ]
      opt.sc.k <- optimalScores(clValid(subsamp, 2:8, clMethods = "kmeans", validation="internal", metric="euclidean", maxitems=900))
      classif2[j, ] <- c(opt.sc.k[1, 3], opt.sc.k[2, 3], opt.sc.k[3, 3])
      print(j)
    }
     
  } else if (method=="hierarchical"){
    for (l in 1:iter){
      subsamp <- dat[sample(row.names(dat), round((perc*nrow(dat))/100, 0), replace=F), ]
      opt.sc.h <- optimalScores(clValid(subsamp, 2:8, clMethods = "hierarchical", validation="internal", method="average", metric="correlation", maxitems=900))
      classif2[l, ] <- c(opt.sc.h[1, 3], opt.sc.h[2, 3], opt.sc.h[3, 3])
      print(l)
    }
    
    } else if (method=="modelbased") {
    for (i in 1:iter){
      subsamp <- dat[sample(row.names(dat), round((perc*nrow(dat))/100, 0), replace=F), ]
      cls <- Mclust(subsamp)
      classif[[1]][[i]] <- cls$classification
      classif[[2]][[i]] <- cls$G
      classif[[3]][[i]] <- cls$modelName
      classif[[4]][[i]] <- cls$loglik
      classif[[5]][[i]] <- cls$parameters$mean
      print(i)
    }
    
    
  } else {
      stop("Unrecognized method.")
    }

  if (method %in% c("kmeans", "hierarchical")){
    classif2
  } else {
    classif
  }
  
}


## STEP 2
getOverlapMB <- function(dat, perc, iter){
  require(mclust)
  require(vcd)
  
  classif1 <- list(1,2,3,4,5)
  for (i in 1:5){
    classif1[[i]] <- vector("list", iter)
  }
  classif2 <- list(1,2,3,4,5)
  for (i in 1:5){
    classif2[[i]] <- vector("list", iter)
  }
  
  crosstabs <- vector("list", iter)
  
  for (i in 1:iter){
    subsamp <- dat[sample(row.names(dat), round((perc*nrow(dat))/100, 0), replace=F), ]
    cls1 <- Mclust(subsamp, G=3)
    classif1[[1]][[i]] <- cls1$classification
    classif1[[2]][[i]] <- cls1$G
    classif1[[3]][[i]] <- cls1$modelName
    classif1[[4]][[i]] <- cls1$loglik
    classif1[[5]][[i]] <- cls1$parameters$mean
    
    cls2 <- Mclust(subsamp, G=5)
    classif2[[1]][[i]] <- cls2$classification
    classif2[[2]][[i]] <- cls2$G
    classif2[[3]][[i]] <- cls2$modelName
    classif2[[4]][[i]] <- cls2$loglik
    classif2[[5]][[i]] <- cls2$parameters$mean
    
    crosstabs[[i]]<- assocstats(table(classif1[[1]][[i]], classif2[[1]][[i]]))
      #table(classif1[[1]][[i]], classif2[[1]][[i]])
      #
    
    
    print(i)
  }
  
  list(classif1,
       classif2,
       crosstabs)
  
}


# 6
extractInfo <- function(MBoutput){
  Groups <- vector()
  for (i in 1:length(MBoutput[[2]])){
    Groups[i]  <- MBoutput[[2]][[i]]
  }
  table(Groups)  
}


# 7
extractInfo2 <- function(MBoutput){
  require(psych)
  Phis <- vector()
  for (i in 1:length(MBoutput[[3]])){
    Phis[i]  <- MBoutput[[3]][[i]]$cramer
  }
  describe(Phis)  
}


# 8
plotCompDens <- function(dat, clus.no){
  x <- princomp(dat, scores=TRUE)$scores
  maxdens <- max(density(x[,1])$y, density(x[,2])$y, density(x[,3])$y)
  plot(density(x[,1]), col="blue",  main="", ylim=c(0,maxdens), xlim=c(-4, 4), lwd=2)
  lines(density(x[,2]), col="magenta",ylim=c(0,maxdens), xlim=c(-4, 4), lwd=2)
  lines(density(x[,3]), col="green",ylim=c(0,maxdens), xlim=c(-4, 4), lwd=2)
  legend("topleft", c("Component 1", "Component 2", "Component 3"), fill=c("blue", "magenta", "green"), cex=0.7)
  title(paste("Cluster", clus.no))
}


# 9
## From clusteval package:
random_clustering <- function(x, K, prob = NULL) {
  if (!is.null(prob)) {
    if (!is.numeric(prob)) {
      stop("The vector 'prob' must be 'numeric'.")
    }
    if (K != length(prob)) {
      stop("The length of 'prob' must equal 'K'.")
    }
    if (sum(prob) != 1) {
      stop("The sum of the probabilities must sum to 1.")
    }
    if (any(prob <= 0) || any(prob >= 1)) {
      stop("The cluster probabilties must be between 0 and 1.")
    }
  }
  sample(x = seq_len(K), size = nrow(x), replace = TRUE, prob = prob)
}


# 10
compare_stimulus_selections <- function(stimulus_lists){
  format_image_codes <- function(codes){
    as.character(formatC(codes, digits=1, format="f"))
  }
  
  for(i in 1:length(stimulus_lists)){
    reformatted_image_codes <- format_image_codes(stimulus_lists[[i]])
    nb_common_stims <- length(which(reformatted_image_codes %in% forclus2$code))
    nb_outliers <- length(which(reformatted_image_codes %in% outlier_index))
    nb_wide_CI <- length(which(reformatted_image_codes %in% wide_CI_index))
    nb_missing_dom <- length(which(reformatted_image_codes %in% missing_dominance_index))
    which_clusters <- table(forclus2[which(forclus2$code %in% reformatted_image_codes), "classif"])
    print(paste("----- For stimulus group - ", names(stimulus_lists)[i], ": -----", sep=""))
    print(paste("Nb. common stimuli: ", 
                nb_common_stims,
                " of the total number the study used, i.e., ",
                length(stimulus_lists[[i]]),
                sep=""))
    print(paste("Of the stimuli outside our clustering solution: ",
                nb_missing_dom, " missing Dom scores + ",
                nb_outliers, " outliers + ",
                nb_wide_CI, " with wide CIs, leaving ", 
                length(stimulus_lists[[i]]) - sum(nb_common_stims, nb_outliers, nb_wide_CI, nb_missing_dom),
                " stimuli unaccounted for.", sep=""))
    print("Common stimuli distributed across my clusters as:")
    print(which_clusters)
  }
}


# 11
grpdist <- function(X)
{
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}
