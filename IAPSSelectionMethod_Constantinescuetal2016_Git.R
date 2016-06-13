
# 0. Necessary packages ---------------------------------------------------

AlreadyInstalled <- row.names(installed.packages())

PackagesRequiredForThisCode <- c("amap", "car", "cclust", "clue", "clusteval", "clValid", "descr", "dismo", "fpc", "GenABEL", "ggplot2", "GMD", "gridExtra", "ICC", "igraph", "ks", "lattice", "latticeExtra", "lsr", "mclust", "MVN", "mvnormtest", "NbClust", "plyr", "psych", "RColorBrewer", "rgl", "shape", "stargazer", "TeachingDemos", "vcd", "vegan", "wordcloud", "xtable")


NewPackages <- PackagesRequiredForThisCode[which(! PackagesRequiredForThisCode %in% AlreadyInstalled)]

install.packages(NewPackages) # *Install* all *completely new* packages.

lapply(PackagesRequiredForThisCode,function(x){require(x,character.only=TRUE)}) 


rm(list=ls())
options(scipen=1000, bitmapType="cairo") # to support transparency in eps files

# Get helper functions from additional script:
source("/home/caterina/Documents/PhD/Publications/ClustAn_Manuscript/IAPSHelpFUNs.R")







# 1. Data prep -----------------------------------------------------------


all <- as.matrix(read.table("/home/caterina/Documents/PhD/Publications/ClustAn_Manuscript/IAPSAll.txt", header=T, sep="\t"))
NA -> all[which(all == "0.00")]

all <- data.frame(all)
all <- data.frame(all[,1:2], apply(all[,3:11], 2, function(x) as.numeric(as.character(x)))) 

dim(all)
# viable(all) # More missing data in the case of dom1mn.

tail(all)
table(all$set)






# 2. Remove duplicates -------------------------------------------------------

table(duplicated(all))
table(duplicated(all$IAPS)) # Iaps codes re-occur for 12 cases.

all$ID <- as.numeric(as.factor(all$IAPS))
tail(all) # Odd: ID counts fewer rows than dim(). Turns out all$IAPS sometimes includes doubles!

doubles <- which(cbind(table(all$ID))==2)
duplicates <- all[all$ID %in% doubles, ]

stargazer(duplicates[, c(1:3,5,7,11)], summary=F)

uniques1 <- matrix(ncol=2, nrow=0)
for (i in seq(0, 22, by=2)){
  uniques1 <- rbind(uniques1, unique(duplicates[ c(1, 2) + i, 1:2]) )
}

uniques2 <- matrix(ncol=10, nrow=0)
for (i in seq(0, 22, by=2)){
  uniques2 <- rbind(uniques2, colMeans(duplicates[ c(1, 2) + i, 3:12]) )
}

uniques <- data.frame(uniques1, uniques2)

# And now to substitute in real data:
all.unique <- all[! all$ID %in% doubles, ]
all.unique <- rbind(all.unique, uniques)
all.unique <- all.unique[order(all.unique$ID), ]

summary(as.matrix(all.unique[, 3:10]))


# 3. Remove missing vals from dom1mn data version----------------------------------------------

# Which to use ? Dom1mn or Dom2mn?
all2 <- all.unique

# Are dom1mn and dom2mn correlated?
for.cor <- na.exclude(all2); dim(for.cor) # N=60
cor(for.cor$dom1mn, for.cor$dom2mn) # highly correlated.

missing_dominance_index <- all2[which(is.na(all2$dom1mn)), "IAPS"]

dim(na.exclude(all2[, -c(9:10)]) ) # dims for data with NAs excluded, when Dom2mn is not considered
all3 <- na.exclude(all2[, -c(9:10)]) 



# 4. Outliers ------------------------------------------------------------

table(is.outlier(all3$valmn))
table(is.outlier(all3$aromn))
table(is.outlier(all3$dom1mn)) # 32 here.

outlier_index <- which(is.outlier(all3$dom1mn))
outlier_index <- all3[outlier_index, "IAPS"]

all3 <- all3[- which(is.outlier(all3$dom1mn)), ]



# 5. 95% CI -------------------------------------------------------------

imprecise.v <- which(CIwidth(all3$valmn, all3$valsd, 100)>1) #7
imprecise.a <- which(CIwidth(all3$aromn, all3$arosd, 100)>1) #46
imprecise.d <- which(CIwidth(all3$dom1mn, all3$dom1sd, 100)>1) #25

table(table(c(imprecise.v, imprecise.a, imprecise.d)))

span.over.1 <- sort(unique(c(imprecise.v, imprecise.a, imprecise.d)))
wide_CI_index <- all3[span.over.1, "IAPS"]

all3<-all3[-span.over.1,]


dim(all3)




forclus<-all3[, c(3,5,7)]

norval<- rntransform(forclus$valmn)
noraro<- rntransform(forclus$aromn) 
nordom<- rntransform(forclus$dom1mn)
nor<-cbind(norval, noraro, nordom)



# 6. Coefficient of variation --------------------------------------------

length(which(cv(all3$valmn, all3$valsd) > 30)) # not good, 406.
length(which(cv(all3$aromn, all3$arosd) > 30)) # 817
length(which(cv(all3$dom1mn, all3$dom1sd) > 30)) # 730.

representative.cases <- which(cv(all3$dom1mn, all3$dom1sd) < 30 & cv(all3$aromn, all3$arosd) < 30 & cv(all3$valmn, all3$valsd) < 30)
length(representative.cases) # Just 1!
all3[representative.cases, ]




# Get measure of multivariate distances from centroid:
mahD <- mahalanobis(forclus, center=colMeans(forclus), cov(forclus) )
summary(mahD)
quantile(mahD, probs=seq(0, 1, 0.1))





#7. Plotting the data ------------------------------------------------------------


#1
P1 <- ggplot(all3, aes(x=valmn, y=aromn)) + geom_point() + stat_smooth(method="lm", color="magenta") + stat_smooth(method="loess") + xlab("Valence") + ylab("Arousal") + ggtitle("(Non-linear) relationship between\nValence and Arousal") + theme(text = element_text(size=20)) + scale_x_continuous(expression(atop("Valence", atop(italic("r = -0.21***"), ""))), breaks=1:9, labels=1:9) + scale_y_continuous("Arousal", breaks=1:9, labels=1:9)
#round(cor(all3$valmn, all3$aromn), 2)


P2 <- ggplot(all3, aes(x=valmn, y=dom1mn)) + geom_point() + stat_smooth(method="lm", color="magenta") + stat_smooth(method="loess") + xlab("Valence") + ylab("Dominance") + ggtitle("Relationship between\nValence and Dominance")  + theme(text = element_text(size=20)) + scale_y_continuous("Dominance", breaks=1:9, labels=1:9) + scale_x_continuous(expression(atop("Valence", atop(italic("r = 0.83***"), ""))), breaks=1:9, labels=1:9)
#round(cor(all3$valmn, all3$dom1mn), 2)


P3 <- ggplot(all3, aes(x=aromn, y=dom1mn)) + geom_point() + stat_smooth(method="lm", color="magenta") + stat_smooth(method="loess") + xlab("Arousal") + ylab("Dominance") + ggtitle("Relationship between\nArousal and Dominance") + theme(text = element_text(size=20)) + scale_x_continuous(expression(atop("Arousal", atop(italic("r = -0.55***"), ""))), breaks=1:9, labels=1:9) + scale_y_continuous("Dominance", breaks=1:9, labels=1:9)
#round(cor(all3$aromn, all3$dom1mn), 2)

postscript("~/Desktop/ScatterPADCorrs.eps", width=6, height=18, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
grid.arrange(P1, P2, P3, ncol=1)
dev.off()


#2
## a)
print(cloud(valmn ~ aromn * dom1mn, data=all3, col=all3$ID, perspective=T, distance=0.4, type="p", shade=T, lex=2, xlab="Arousal", ylab="Dominance", zlab="Valence", main="3D Scatterplot", R.mat = diag(4), screen = list(z=-40, x=-42), aspect = c(1,1)),  position = c(0, 0, 0.5, 1), more = TRUE)

## b)
print(splom( ~ data.frame(valmn, aromn, dom1mn), groups=ID, data=all3, main="2D visualization", varnames = c("Valence", "Arousal", "Dominance"),
             panel = function(x, y, ...) {
               panel.xyplot(x, y, ...)
               ok <- (x > 0) & (y > 0)
               fm <- lm(  log(exp(y[ok]) - 1) ~  log(exp(x[ok]) - 1) )
               panel.abline(fm, ...)
             }),  position = c(0.5, 0, 1, 1), more=F)


#3
PANA <- kmeans(forclus, 2)$cluster
PANAsym <- ifelse(PANA==1, 2, 3)
PANAcols <- ifelse(PANA==1, "#fdb863", "#e66101")
Cl2 <- cloud(valmn ~ aromn * dom1mn, data=all3, 
             col=PANAcols, pch=PANA, lex=1.5, 
             main=list("3D structure of IAPS data", cex=2), 
             perspective=T, distance=0.4, 
             type="p",  cex=0.7, shade=T, lex=2, 
             xlab=list("Arousal", rot=43, cex=2), 
             ylab=list("Dominance", rot=-50, cex=2), 
             zlab=list("Valence", rot=96, cex=2),
             R.mat = diag(4), screen = list(x=-42, y=35, z=25), aspect = c(1,1) )

PANA <- kmeans(forclus, 3)$cluster
PANAcols <- ifelse(PANA==1, "#fdb863", ifelse(PANA==2, "#e66101", "#5e3c99"))
Cl3 <- cloud(valmn ~ aromn * dom1mn, data=all3, 
             col=PANAcols, pch=PANA, lex=1.5, 
             main=list("3D structure of IAPS data", cex=2), 
             perspective=T, distance=0.4, 
             type="p",  cex=0.7, shade=T, lex=2, 
             xlab=list("Arousal", rot=43, cex=2), 
             ylab=list("Dominance", rot=-50, cex=2), 
             zlab=list("Valence", rot=96, cex=2),
             R.mat = diag(4), screen = list(x=-42, y=35, z=25), aspect = c(1,1) )


postscript("~/Desktop/PANA.eps", width=6, height=12, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
c(Cl2, Cl3, layout=c(1, 2))
dev.off()




# 8. Exploring the data ---------------------------------------------------

# Some descriptives
summary(all3[,c('valmn', 'aromn', 'dom1mn')])
describe(all3$valmn); describe(all3$aromn); describe(all3$dom1mn)


# Assessing univariate normality
histkdnc(all3$valmn)
shapiro.test(all3$valmn)

histkdnc(all3$aromn)
shapiro.test(all3$aromn)

histkdnc(all3$dom1mn)
shapiro.test(all3$dom1mn)

# Correlations between dimensions
round(cor(all3[, c(3,5,7)]), 3)  # same corrs as shown in plots above.


# How about multivariate distrib?
rgl.viewpoint(plot(kde(forclus), 
                   display="persp", cont=c(seq(0, 100, by=15)),
                   drawpoints=F, labcex=1.5,
                   xlim=c(1,9),ylim=c(1,9),zlim=c(1,9),
                   xlab="Valence", ylab="Arousal", zlab="Dominance")
              , fov = 50, zoom = 5)
grid3d(side=c("x", "y", "z"), at=c(1:9))
#rgl.snapshot( "/home/caterina/Documents/PhD/Publications/3Ddens", fmt="png", top=TRUE)


# and as opposed to random data generated to be normal:
a <- rnorm(900, mean=5, sd =1.5)
b <- rnorm(900, mean=5, sd =1.5)
c <- rnorm(900, mean=5, sd =1.5)
rgl.viewpoint(plot(kde(data.frame(a,b,c)), 
                   display="persp", cont=c(seq(0, 100, by=15)),
                   drawpoints=F, labcex=1.5,
                   xlim=c(1,9),ylim=c(1,9),zlim=c(1,9),
                   xlab="X", ylab="Y", zlab="Z")
              , fov = 50, zoom = 5)
grid3d(side=c("x", "y", "z"), at=1:9)
# rgl.snapshot( "/home/caterina/Documents/PhD/Publications/3Ddensnormal", fmt="png", top=TRUE)

# Also:
# as opposed to:
P1 <- cloud(forclus[,1] ~ forclus[,2]*forclus[,3], main="Comparison of original, normalised, and random normal data",
            xlab="Arousal", ylab="Dominance", zlab="Valence") # main="Original data"
P2 <- cloud(nor[,1] ~ nor[,2]*nor[,3], main="Normalised data",
            xlab="Arousal", ylab="Dominance", zlab="Valence")
P3 <- cloud(a ~ b*c, main="Random normal",
            xlab="Arousal", ylab="Dominance", zlab="Valence")
c(P1,P2,P3, layout=c(3, 1))




# Multivariate normality for the full data:
mardiaTest(forclus, qqplot=T)
mardiaTest(nor)
hzTest(forclus, qqplot=T)
hzTest(nor)
roystonTest(forclus, qqplot=T)
mshapiro.test(t(forclus))

result1 <- hzTest(forclus[, 1:2]) 
result2 <- hzTest(forclus[, 2:3]) 
result3 <- hzTest(forclus[, c(1,3)]) 
mvnPlot(result1, type = "persp", default = TRUE ) 
mvnPlot(result2, type = "persp", default = TRUE )
mvnPlot(result3, type = "persp", default = TRUE ) 





# 9. *** \\// K-MEANS CLUSTERING \\// *** --------------------------------------------------


# How many clusters...? 

NbClust(forclus, method="kmeans") # Max nb of clusters by default = 15
NbClust(forclus, method="kmeans", max.nc=30)

# 1. cl_validity
K <- list()
K.res <- vector()
for (i in 2:21){
  print(i)
  K[[i-1]] <- Kmeans(as.matrix(forclus), i, iter.max=2000, nstart = 1000, method="euclidean")
}

for(i in 1:20){
  print(paste("k=", i+1, sep="")); print(cl_validity(K[[i]], dist(forclus, "euclidean")))
  K.res <- c(K.res, cl_validity(K[[i]], dist(forclus, "euclidean") ) ) 
} # Only seems to increase with k, which is not unexpected. 

# Check bend:
Output <- data.frame(NbClus=2:21, Res=as.numeric(K.res))
postscript("~/Desktop/DissimilarityElbow.eps", width=8, height=8, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
ggplot(Output, aes(x=NbClus, y=Res)) + geom_point(size=4) + geom_path() + xlab(expression(italic("k"))) + ylab("Dissimilarity accounted for") + ggtitle("Amount of dissimilarity explained by clusters\n") + theme(text = element_text(size=20)) + ylim(c(0.3, 0.8))
dev.off()

"--------------------------------------------------------------------------"

# 2. clustIndex

nReps <- 100000
nMaxCl <- 8
criteria <- vector("list", nReps)
nCrit <- 4
critNames <- c( "calinski", "hartigan", "ball", "ssi")

for(rep in 1:nReps){
  criteria[[rep]] <- matrix(ncol=nCrit, nrow=nMaxCl-1, dimnames=list(2:nMaxCl, critNames))
}


for(rep in 1:nReps){
  print(rep)
  
  K <- list()
  
  for (ncl in 2:nMaxCl){
    K[[ncl-1]] <- cclust(as.matrix(forclus), ncl, dist="euclidean", method= "kmeans")
    criteria[[rep]][ncl-1, ] <- cclust::clustIndex(K[[ncl-1]], forclus, index=c("calinski", "hartigan", "ball", "ssi")) 
  }
}

finalUnsc <- data.frame(Reduce("+", criteria) / length(criteria)) # average matrix

Cal <- ggplot(finalUnsc, aes(x=2:8, y=calinski)) + geom_point(size=3.5) + geom_path() + ggtitle("Calinski") + ylab("Index value") + xlab("k") + theme(text = element_text(size=20)) + scale_x_continuous(expression(italic("k")), breaks=2:8, labels=2:8) + ylim(c(800, 1200))

Hart <- ggplot(finalUnsc, aes(x=2:8, y=hartigan)) + geom_point(size=3.5) + geom_path() + ggtitle("Hartigan") + ylab("Index value") + xlab("k") + theme(text = element_text(size=20)) + scale_x_continuous(expression(italic("k")), breaks=2:8, labels=2:8) + ylim(c(0, 2))

Ball <- ggplot(finalUnsc, aes(x=2:8, y=ball)) + geom_point(size=3.5) + geom_path() + ggtitle("Ball")+ ylab("Index value") + xlab("k") + theme(text = element_text(size=20)) + scale_x_continuous(expression(italic("k")), breaks=2:8, labels=2:8) + ylim(c(0, 1000))

SSI <- ggplot(finalUnsc, aes(x=2:8, y=ssi)) + geom_point(size=3.5) + geom_path() + ggtitle("SSI")+ ylab("Index value") + xlab("k") + theme(text = element_text(size=20)) + scale_x_continuous(expression(italic("k")), breaks=2:8, labels=2:8) + ylim(c(0.4, 0.6))


postscript("~/Desktop/ClustCritSelection.eps", width=9, height=9, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
grid.arrange(Cal, Hart, Ball, SSI, ncol=2)
dev.off()




# Or :
Connect <- vector()
Dunn <- vector()
Sil <- vector()
for(rep in 1:100){
  print(rep)
  optSc <- optimalScores(clValid(forclus, 2:8, clMethods="kmeans", validation="internal", maxitems=1000, metric="euclidean"))
  Connect[rep] <- as.character(optSc[which(row.names(optSc) == "Connectivity"), "Clusters"])
  Dunn[rep] <- as.character(optSc[which(row.names(optSc) == "Dunn"), "Clusters"])
  Sil[rep] <- as.character(optSc[which(row.names(optSc) == "Silhouette"), "Clusters"])
}

table(Connect)
table(Dunn)
table(Sil)
# 8 or 2 clusters


# For future use:
Kmeans.5 <-  Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="euclidean")

# Checking other methods for computing distances
Kmeans.51 <- Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="correlation")
Kmeans.52 <- Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="pearson")
Kmeans.53 <- Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="manhattan")
Kmeans.54 <- Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="pearson")
Kmeans.55 <- Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="maximum")
# However, euclidean is considered the most appropriate option for k-means

cluster_similarity(Kmeans.5$cluster, Kmeans.51$cluster, similarity="rand")
cluster_similarity(Kmeans.5$cluster, Kmeans.52$cluster, similarity="rand")
cluster_similarity(Kmeans.5$cluster, Kmeans.53$cluster, similarity="rand")
cluster_similarity(Kmeans.5$cluster, Kmeans.54$cluster, similarity="rand")
cluster_similarity(Kmeans.5$cluster, Kmeans.55$cluster, similarity="rand")
#They seem to not vary very much according to the distance method used


cl_validity(Kmeans.55, dist(forclus, "euclidean") )
cl_validity(Kmeans.55, dist(forclus, "manhattan") ) # correlation distances not available. Kmeans with euclidean distances is doing well however.
cl_validity(Kmeans.55, dist(forclus, "maximum") )



# 10. *** \\// HIERARCHICAL CLUSTERING \\// *** ------------------------------------------------------------

hdist <- Dist(forclus, method="correlation", upper=T)
NbClust(forclus, method="average", diss=hdist, distance=NULL) # Max nb of clusters by default = 15
NbClust(forclus, method="average", max.nc=30, diss=hdist, distance=NULL)



# 10A. Cophenetic correlations / Gower distances --------------------------
# From: Borcard, D., Gillet, F., and Legendre, P. (2011). Numerical ecology with R. New York: Springer.


# Create clusters with various settings:
hc.w.euc <- hcluster(forclus, method="euclidean", link="ward"); plot(hc.w.euc)
hc.w.cor <- hcluster(forclus, method="correlation", link="ward"); plot(hc.w.cor)

hc.av.euc <- hcluster(forclus, method="euclidean", link="average"); plot(hc.av.euc)
hc.av.cor <- hcluster(forclus, method="correlation", link="average"); plot(hc.av.cor) #*** Best.

hc.com.euc <- hcluster(forclus, method="euclidean", link="complete"); plot(hc.com.euc)
hc.com.cor <- hcluster(forclus, method="correlation", link="complete"); plot(hc.com.cor)

hc.sin.euc <- hcluster(forclus, method="euclidean", link="single"); plot(hc.sin.euc)
hc.sin.cor <- hcluster(forclus, method="correlation", link="single"); plot(hc.sin.cor)


"--------------------------------------------------------------------------"

# With Ward linkage:
cor(Dist(forclus, method="euclidean"), cophenetic(hc.w.euc)) # r=0.72
cor(Dist(forclus, method="correlation"), cophenetic(hc.w.cor)) # r=0.82 

# Average linkage:
cor(Dist(forclus, method="euclidean"), cophenetic(hc.av.euc)) # r=0.71
cor(Dist(forclus, method="correlation"), cophenetic(hc.av.cor)) # r=0.91 !

# Complete linkage:
cor(Dist(forclus, method="euclidean"), cophenetic(hc.com.euc)) # r=0.72
cor(Dist(forclus, method="correlation"), cophenetic(hc.com.cor)) # r=0.77

# Single linkage:
cor(Dist(forclus, method="euclidean"), cophenetic(hc.sin.euc)) # r=0.22
cor(Dist(forclus, method="correlation"), cophenetic(hc.sin.cor)) # r=0.87


# Average linkage with distances based on correlations seem to be the best option here (stays closest to the original data). 
# Also, it seems distance measures based on corrs (rather than Euclidean distances) present higher corrs overall.


# 10B. Gower distances -----------------------
sum((Dist(forclus, method="correlation") - cophenetic(hc.w.cor))^2)
sum((Dist(forclus, method="correlation") - cophenetic(hc.av.cor))^2) # Smallest value! Average linkage + correlation = best match
sum((Dist(forclus, method="correlation") - cophenetic(hc.com.cor))^2)
sum((Dist(forclus, method="correlation") - cophenetic(hc.sin.cor))^2)



# 10C. The "Elbow" Method for Clustering Evaluation -----------------------

Step1 <- dist(forclus,method="euclidean",diag=FALSE,upper=FALSE)
Step2<- hclust(Step1, method = "average", members=NULL)
Step3<- css.hclust(Step1, hclust.obj=Step2)
Step4 <- elbow.batch(Step3, 0.05, 0.80, precision=3)
print(Step4) # Suggests k = 7
# And:
plot(Step3$k, Step3$ev, main="Variance explained by clusters", xlab="Nb of clusters", ylab="Explained variance", sub="(Valence, arousal & dominance used for clustering)")


#Normalised version
Step1 <- dist(nor, method="euclidean",diag=FALSE,upper=FALSE)
Step2<- hclust(Step1, method = "average", members=NULL)
Step3<- css.hclust(Step1, hclust.obj=Step2)
Step4 <- elbow.batch(Step3, 0.05, 0.80, precision=3)
print(Step4) # Suggests k = 15 ! The 'normalization' process seems to rather have interfered with finding a clustering structure here.
# And:
plot(Step3$k, Step3$ev, main="Variance explained by clusters", xlab="Nb of clusters", ylab="Explained variance", sub="(Valence, arousal & dominance used for clustering)")



# 10D. Cluster indices:: clValid ------------------------

hierk1 <- clValid(forclus, 2:8, clMethods=c("hierarchical"), validation=c("internal"), maxitems=10000, method="average", metric="correlation")
optimalScores(hierk1) # 2 or 3 
plot(hierk1)



# 10E. Average silhouette widths for most appropriate model----------------
# (i.e. Ave. linkage + correlation distaces)

asw <- numeric(nrow(forclus))
for (k in 2:(nrow(forclus)-1)) {
  sil <- silhouette(cutree(hc.av.cor, k=k), Dist(forclus, method="correlation"))
  asw[k] <- summary(sil)$avg.width
}


k.best <- which.max(asw)
# windows(title="Silhouettes - Average - k = 2 to n-1")
plot(1:nrow(forclus), asw, type="h", xlim=c(0,20), 
     main="Silhouette-optimal \n number of clusters \n Average", 
     xlab="k (number of groups)", ylab="Average silhouette width") 
axis(1, k.best, paste("optimum", k.best,sep="\n"), col="red", font=2, col.axis="red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k =", k.best, "\n", 
    "with an average silhouette width of", max(asw), "\n")
# Same as clValid output. However, not reasonable given large amount of data.


# 10F. Mantel Optimality ---------------------------------------------------

kt <- data.frame(k=1:nrow(forclus), r=0)
for (i in 2:30) {
  gr <- cutree(hc.av.cor, i)
  distgr <- grpdist(gr)
  mt <- cor(Dist(forclus, method="correlation"), distgr, method="pearson")
  kt[i,2] <- mt
}


k.best <- which.max(kt$r)

plot(kt$k, kt$r, type="h", 
     main="Mantel-optimal number of clusters - Average", 
     xlim=c(0,30),
     xlab="k (number of groups)", 
     ylab="Pearson's correlation")
axis(1, k.best, paste("optimum", k.best, sep="\n"), col="red", font=2, col.axis="red")
points(k.best, max(kt$r), pch=16, col="red", cex=1.5)
cat("", "Mantel-optimal number of clusters k =", k.best, "\n", "with a matrix linear correlation of", max(kt$r), "\n")
# Suggests k=2.


# 11. *** \\// MODEL-BASED CLUSTERING \\// *** -----------------------------------------------
PADdims <- forclus
names(PADdims) <- c("Valence", "Arousal", "Dominance")
  
fit <- Mclust(forclus)
fit2 <- Mclust(PADdims) # This is with proper column names
t(fit$parameters$mean) # The means for each component
# write.table(round(t(fit$parameters$mean), 3), sep="&", file="IAPSClustercenters.txt")

dens <- densityMclust(PADdims)
dens$density

summary(densityMclust(PADdims))

postscript("~/Desktop/UncertPlot.eps", width=6, height=6, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
plot(fit2, what="uncertainty") # plot results ; 
dev.off()


plot(fit2, what="BIC") # plot results ; 
BIC <- fit$BIC[1:9, ]
min(BIC)
max(BIC)

plot(fit2, what="classification") # plot results ; 
print(fit) 
table(fit$classification) # Cluster N.
clPairs(forclus, classification = fit$classification, symbols = 1:5)
# The 5 clusters from Mclust()





forclus2 <- data.frame(forclus, unc=fit$uncertainty, classif=fit$classification)
round(aggregate(unc ~ classif, data=forclus2, FUN=mean), 3)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C1 <- subset(forclus, fit$classification==1)
C2 <- subset(forclus, fit$classification==2)
C3 <- subset(forclus, fit$classification==3)
C4 <- subset(forclus, fit$classification==4)
C5 <- subset(forclus, fit$classification==5)


# Component normality: rgl package

grp1 <- c(rep("C1", dim(C1)[1]), rep("C5", dim(C5)[1]))
grp2 <- c(rep("C2", dim(C2)[1]), rep("C3", dim(C3)[1]), rep("C4", dim(C4)[1]))

# After manual rotation of rgl plot, when suitable angle has been found:
pp <- par3d(no.readonly=TRUE)
pp$cex <- 1.5

rgl.viewpoint(plot(kde(rbind(C1, C5)),
                   approx.cont=F,
                   y.group=grp1,
                   # colors=heat.colors(6),
                   display="persp", cont=c(seq(0, 100, by=15)),
                   drawpoints=F, cex=10,
                   xlim=c(1,9),ylim=c(1,9),zlim=c(1,9),
                   zoom=3,
                   xlab="Valence", ylab="Arousal", zlab="Dominance") )
grid3d(side=c("x", "y", "z"), at=1:9)
par3d(pp) # retrieve appropriate angle consistent with previous selection

# rgl.snapshot( "/home/caterina/Desktop/3DdensClus1and5.png", fmt="png", top=TRUE)
rgl.postscript("/home/caterina/Desktop/3DdensClus1and5.eps","eps",drawText=FALSE)


rgl.viewpoint(plot(kde(rbind(C2,C3, C4)),
                   approx.cont=F,
                   y.group=grp2,
                   # colors=heat.colors(6),
                   display="persp", cont=c(seq(0, 100, by=15)),
                   drawpoints=F, cex=10,
                   xlim=c(1,9),ylim=c(1,9),zlim=c(1,9),
                   xlab="Valence", ylab="Arousal", zlab="Dominance"))
grid3d(side=c("x", "y", "z"), at=1:9)
par3d(pp)

# rgl.snapshot( "/home/caterina/Desktop/3DdensClus2,3,4.png", fmt="png", top=TRUE) 
rgl.postscript("/home/caterina/Desktop/3DdensClus2,3,4.eps","eps", drawText=FALSE)






# Getting results from normality tests, for each component:
normalityTests <- list(1,2,3,4,5)
for (i in 1:5){
  normalityTests[[i]] <- list(1:5)  
}

Clusters <- list(C1, C2, C3, C4, C5)

for (k in 1:5){
  normalityTests[[1]][[k]] <- mardiaTest(Clusters[[k]], qqplot=F)
  normalityTests[[2]][[k]] <- mshapiro.test(t(Clusters[[k]]))   
  normalityTests[[3]][[k]] <- hzTest(Clusters[[k]], qqplot=F)  
  normalityTests[[4]][[k]] <- roystonTest(Clusters[[k]], qqplot=F)  
  normalityTests[[5]][[k]] <- pca(Clusters[[k]])   
}



# Now to print out a table:
mardia <- matrix(ncol=5, nrow=4)
for (k in 1:5){
  # Mardia test
  mardia[1,k] <- normalityTests[[1]][[k]]@chi.skew
  mardia[2,k] <- normalityTests[[1]][[k]]@p.value.skew
  mardia[3,k] <- normalityTests[[1]][[k]]@z.kurtosis
  mardia[4,k] <- normalityTests[[1]][[k]]@p.value.kurt
}; mardia <- round(mardia, 3)


shapiro <- matrix(ncol=5, nrow=2)
for (k in 1:5){
  # S-W
  shapiro[1,k] <- normalityTests[[2]][[k]]$statistic
  shapiro[2,k] <- normalityTests[[2]][[k]]$p.value
}; shapiro <- round(shapiro,3)


hz <- matrix(ncol=5, nrow=2)
for (k in 1:5){
  # Henze-Zirkler
  hz[1,k] <- normalityTests[[3]][[k]]@HZ
  hz[2,k] <- normalityTests[[3]][[k]]@p.value
}; hz <- round(hz,3)


royston <- matrix(ncol=5, nrow=2)
for (k in 1:5){
  # Royston
  royston[1,k] <- normalityTests[[4]][[k]]@H
  royston[2,k] <- normalityTests[[4]][[k]]@p.value
}; royston <- round(royston,3)


PCA <- matrix(ncol=5, nrow=3)
for (k in 1:5){
  # PCA
  PCA[1,k] <- normalityTests[[5]][[k]]$eig[1]
  PCA[2,k] <- normalityTests[[5]][[k]]$eig[2]
  PCA[3,k] <- normalityTests[[5]][[k]]$eig[3]
}; PCA <- round(PCA,3)

# write.table(rbind(mardia, shapiro, hz, royston, PCA), sep="&" , row.names=F, file="MVNTests.txt")



postscript("~/Desktop/DensitiesForEachComponent.eps", width=5, height=5, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
plotCompDens(C1,1)
plotCompDens(C2,2)
plotCompDens(C3,3)
plotCompDens(C4,4)
plotCompDens(C5,5)
dev.off()









# Best reps per cluster:

# Display only 5 (or x) cases / cluster :
dat.uncert <- data.frame(fit$data, fit$uncertainty, fit$classification)
dat.uncert <- dat.uncert[order(dat.uncert$fit.classification, dat.uncert$fit.uncertainty), ]
# x <- 20
x <- 5
prototypes <- rbind(subset(dat.uncert, fit.classification == 1)[1:x, ],
                    subset(dat.uncert, fit.classification == 2)[1:x, ],
                    subset(dat.uncert, fit.classification == 3)[1:x, ],
                    subset(dat.uncert, fit.classification == 4)[1:x, ],
                    subset(dat.uncert, fit.classification == 5)[1:x, ] )

colIndex2 <- mapvalues(prototypes$fit.classification, 1:5, c("#d7191c", "#fdae61", "darkgreen", "#abdda4", "#2b83ba") )

prototypes$img.descr <- paste(as.character(all3[ row.names(prototypes), "desc"]), 
                              as.character(all3[ row.names(prototypes), "IAPS"]), sep="\n")


# Most "likely" cases per cluster:
postscript("~/Desktop/MclustSmallClusters.eps", width=7, height=7, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
cloud(valmn ~ aromn * dom1mn, data=prototypes, 
      col=colIndex2, 
      pch=prototypes$fit.classification, perspective=T, distance=0.4, 
      type="p",  cex=1.5, shade=T, lex=2, 
      xlab=list("Arousal", rot=-25, cex=2), 
      ylab=list("Dominance", rot=38, cex=2), 
      zlab=list("Valence", rot=95, cex=2),
      R.mat = diag(4), screen = list(z=-40, x=-80), aspect = c(1,1) )
dev.off()



# postscript("~/Desktop/WordplotX.eps", width=12, height=11, colormodel="cmyk", onefile=F, paper="special", horizontal=F)
pdf("~/Desktop/Wordplot.pdf", width=12, height=11, colormodel="cmyk", onefile=F, paper="special")
par(mar= c(5, 5, 4, 2) )
plot(prototypes[,1], prototypes[,3], type="n", ylim=c(1,9), xlim=c(1,9),
     xlab="Valence", ylab="Arousal", main="Scatterplot of best 5 representatives from each cluster",
     cex.lab=2.5, cex.main=2.6, cex.axis=2.2)
abline(v=seq(0,10), h=seq(0,10), col="black", lty="dotted")
textplot(prototypes[,1], prototypes[,2], 
         as.character(prototypes$img.descr),
         col=colIndex2,
         lwd=2,
         cex=2.2, font=2,
         ylim=c(1,9), xlim=c(1,9),
         new=F)
legend("bottomright", legend=c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), fill=unique(colIndex2), bty="n", cex=2, text.font=1)
dev.off()





pdf("~/Desktop/Boxplots.pdf", width=4.5, height=10, colormodel="cmyk", onefile=F, paper="special")
# postscript("~/Desktop/Boxplots.eps", width=4.5, height=10, colormodel="CMYK", onefile=F, paper="special", horizontal=F)
par(mfrow=c(3,1), mar=c(3.8, 5, 4, 2)) # A numerical vector of the form c(bottom, left, top, right)
boxplot(forclus$valmn ~ fit$classification, las=1, 
        ylab="Valence", cex.lab=1.8, cex.axis=1.5, 
        col=unique(colIndex2), varwidth=TRUE); grid(col="darkgray")
boxplot(forclus$valmn ~ fit$classification, las=1, 
        xaxt="n", yaxt="n", col=unique(colIndex2), varwidth=TRUE, add=TRUE)

boxplot(forclus$aromn ~ fit$classification, las=1, 
        cex.lab=1.8, cex.axis=1.5, ylab="Arousal", 
        col=unique(colIndex2), varwidth=TRUE); grid(col="darkgray")
boxplot(forclus$aromn ~ fit$classification, las=1, 
        xaxt="n", yaxt="n", col=unique(colIndex2), varwidth=TRUE, add=TRUE)

boxplot(forclus$dom1mn ~ fit$classification, las=1,
        cex.lab=1.8, cex.axis=1.5, ylab="Dominance", 
        col=unique(colIndex2), xlab="Cluster", varwidth=TRUE); grid(col="darkgray")
boxplot(forclus$dom1mn ~ fit$classification, las=1, 
        xaxt="n", yaxt="n", col=unique(colIndex2), varwidth=TRUE, add=TRUE)
par(mfrow=c(1,1))
title("Dimensions broken down by cluster", outer=F)
dev.off()




# 12. ***\\//  VALIDATION \\// *** ------------------------------------------------------


# 12A. ANOVA / Comparison between methods-----------------------------------------


Kmeans.5 <-  Kmeans(as.matrix(forclus), 5, iter.max=2000, nstart = 1000, method="euclidean")
hc.av.cor <- hcluster(forclus, method="correlation", link="average")
fit <- Mclust(forclus)

modelling.dat <- cbind(all3, 
                       alldims = all3$valmn*all3$aromn*all3$dom1mn,
                       clasMB=as.factor(fit$classification),
                       clasKM=as.factor(Kmeans.5$cluster),
                       clasHIER=as.factor(cutree(hc.av.cor, 5)) )


# Set the desired outcome:
outcome <- "alldims" # To be changed as necessary.

# Single predictor models:
M1a <- lm(get(outcome) ~ factor(clasMB), modelling.dat); summary(M1a)$r.squared
M1b <- lm(get(outcome) ~ factor(clasKM), modelling.dat); summary(M1b)$r.squared
M1c <- lm(get(outcome) ~ factor(clasHIER), modelling.dat); summary(M1c)$r.squared
# 

                    #model-b;  #k-means;   #hierarchical
meanRsq <- matrix(c(0.7175725, 0.8084096, 0.5947314,  # alldims
                    0.8852297, 0.8850551, 0.7539337,    # valmn
                    0.4304445, 0.7156515, 0.60237,      # aromn 
                    0.743523,  0.8163968, 0.73796), ncol=3, byrow=T) # dom1mn

mean(meanRsq)

colMeans(meanRsq) # hmm.

sd(meanRsq)


# K-means/Hierarchical vs. model-based clustering:
M2a <- lm(get(outcome) ~ factor(clasMB) + factor(clasKM), modelling.dat)
M2b <- lm(get(outcome) ~ factor(clasMB) + factor(clasHIER), modelling.dat)
M2c <- lm(get(outcome) ~ factor(clasKM) + factor(clasHIER), modelling.dat)

summary(M2a)$r.squared 
summary(M2b)$r.squared 
summary(M2c)$r.squared 

summary(M2a)$r.squared - summary(M1a)$r.squared 
anova(M1a)$'Pr(>F)'[1]  ; anova(M2a)$'Pr(>F)'[1]

summary(M2b)$r.squared - summary(M1c)$r.squared 


## kmeans vs model-based
anova(M1a, M2a) # Adding KMeans improves fit 
anova(M1b, M2a) # Adding MB to Kmeans also improves fit, so they are complementary methods

## hierarchical vs model-based
anova(M1a, M2b) # Adding hierarchical improves fit 
anova(M1c, M2b) # Adding MB to hierarchical also improves fit, so they are complementary methods


# lsr / effect sizes
etaSquared(aov(get(outcome) ~ factor(clasMB), modelling.dat), type=3, anova=T)
etaSquared(aov(get(outcome) ~ factor(clasKM), modelling.dat), type=3, anova=T)
etaSquared(aov(get(outcome) ~ factor(clasHIER), modelling.dat), type=3, anova=T)


# Do the different types of clustering provide correlated results though?
tab1 <- table(fit$classification, Kmeans.5$cluster, deparse.level=2); tab1
chisq.test(tab1,simulate.p.value = T, B = 2000)
assocstats(tab1) 


tab2 <- table(fit$classification, cutree(hc.av.cor, 5), deparse.level=2); tab2
chisq.test(tab2,simulate.p.value = T, B = 2000)
assocstats(tab2) 






# 12B. Randomly removing 10% and reassessing cluster number --------------------



# This section is currently under reconstruction and will be added shortly.
# Thank you for your patience.




# 12C.  Making 2 equal-sized samples: one to compute k, the other to validate it. -------


half1 <- forclus[sample(row.names(forclus), as.integer(nrow(forclus)/2), replace=F), ]
remainder <- row.names(forclus)[! row.names(forclus) %in% row.names(half1)]
half2 <- forclus[sample(remainder, as.integer(nrow(forclus)/2), replace=F), ]

h1 <- Mclust(half1, G=5)
h2 <- Mclust(half2, G=5)

# Predicting cluster structure in one half based on the other half of the data:
pred1<-cl_predict(h1, newdata=half2, type="class_ids")
assocstats(table(pred1$classification, h2$classification))$cramer
  
pred2<-cl_predict(h2, newdata=half1, type="class_ids")
assocstats(table(pred2$classification, h1$classification))$cramer

mean(assocstats(table(pred1$classification, h2$classification))$cramer,
assocstats(table(pred2$classification, h1$classification))$cramer)

# Predicting cluster structure in sample 2, based on sample 1 (where k=5)

pred.cl1 <- predict(h1, newdata=half2)
plot(half2, col=pred.cl1$classification)
assocstats(table(pred.cl1$classification, h1$classification))$cramer

pred.cl2 <- predict(h2, newdata=half1)
plot(half1, col=pred.cl2$classification)
assocstats(table(pred.cl2$classification, h2$classification))$cramer



# 12D. Overlap ----------------------------------------------------------------

clusteval::cluster_similarity(fit$classification, Kmeans.5$cluster, similarity = "rand")
clusteval::cluster_similarity(fit$classification, cutree(hc.av.cor, 5), similarity = "rand")
cluster_similarity(fit$classification, random_clustering(forclus, 5, prob = NULL) )


comembership_table(fit$classification, Kmeans.5$cluster)
comembership_table(fit$classification, cutree(hc.av.cor, 5))



mclust::adjustedRandIndex(fit$classification, Kmeans.5$cluster)
mclust::adjustedRandIndex(fit$classification, cutree(hc.av.cor, 5))
mclust::adjustedRandIndex(fit$classification, random_clustering(forclus, 5, prob = NULL))

igraph::compare(fit$classification, Kmeans.5$cluster, method = c("vi")) #  1.186996
igraph::compare(fit$classification, cutree(hc.av.cor, 5), method = c("vi")) # 1.462122

# max = log(nclust)*2 = 3.219
1.186996 / (log(5)*2) # OR
1.186996 / log(nrow(forclus)) # 0.1760062

1.462122 / (log(5)*2) # OR
1.462122  / log(nrow(forclus)) # 0.2168015


igraph::compare(Kmeans.5$cluster, cutree(hc.av.cor, 5), method = c("vi")) # 1.293943
1.293943 / log(nrow(forclus)) # 0.1918641




# 13. ***\\// Comparing partition with random research paper \\//***---------------------


# Study 1: 
# http://www.biomedcentral.com/1471-2202/8/83/
# Effects of affective picture viewing on postural control
# John F Stins and Peter J Beek	 

forclus2 <- data.frame(forclus, unc=fit$uncertainty, classif=fit$classification, code=all3$IAPS)
forclus2$code <- as.character(all3$IAPS)
image_codes_with_decimals <- forclus2$code[which(substr(forclus2$code, 6, 6) != 0)] 
image_codes_with_decimals

t(fit$parameters$mean) # The means for each component

# NEUTRALS
neutral_faces <- c(2190, 2200, 2210, 2214, 2215, 2221, 2270, 2271, 2280, 2383, 2440, 2512, 2516, 2570)
neutral_household <- c(7000, 7002, 7004, 7009, 7010, 7025, 7030, 7035, 7050, 7052, 7060, 7090, 7150, 7175, 7211)
# POSITIVES
erotic <- c(4002, 4180, 4210, 4232, 4250, 4255, 4460, 4510, 4520, 4531, 4572, 4607, 4608, 4652, 4659, 4670, 4800)
family <- c(2057, 2070, 2080, 2165, 2260, 2311, 2340, 2341, 2360, 2387, 2391, 2660)
# NEGATIVES
mutilation <- c(3000, 3010, 3053, 3060, 3064, 3080, 3100, 3110, 3130, 3150, 3170)
fear <- c(1050, 1120, 1200, 1201, 1300, 1301, 1930, 1932, 3022, 3550, 6190, 6230, 6250, 6260, 6300, 6313, 6350, 6560)

Stins_and_Beek_2007 <- list(neutral_faces=neutral_faces,
                            neutral_household=neutral_household,
                            erotic=erotic,
                            family=family,
                            mutilation=mutilation,
                            fear=fear)

compare_stimulus_selections(Stins_and_Beek_2007)



# Study 2:
# Emotional reactivity in nonsuicidal self-injury: Divergence between self-report and startle measures
# Catherine R. Glenn, Terry D. Blumenthal, E. David Klonsk, Greg Hajcak
# http://www.sciencedirect.com/science/article/pii/S0167876011000791

pleasant <- c(1463, 1710, 1811, 2070, 2080, 2092, 2165, 2311, 2340, 4180, 4460, 4651, 4659, 4660, 4669, 4810, 7325, 8461)
neutral <- c(2320, 2570, 2580, 2870, 5390, 5410, 5532, 5534, 5731, 7009, 7010, 7025, 7041, 7140, 7175, 7224, 7235, 7550)
unpleasant <- c(1050, 1300, 3261, 3500, 3530, 6230, 6250, 6313, 6510, 6560, 6571, 9250, 9253, 9400, 9405, 9410, 9420, 9433)

Glenn_et_al_2011 <- list(pleasant=pleasant, neutral=neutral, unpleasant=unpleasant)
compare_stimulus_selections(Glenn_et_al_2011)



# Study 3:
# Divergent Trajectories in the Aging Mind: Changes in Working Memory for Affective Versus Visual Information With Age
# Joseph A. Mikels, Gregory R. Larkin, Patricia A. Reuter-Lorenz, and Laura L. Carstensen

negative <- c(1030, 9594,  9582, 1274, 9584, 1051,  1220, 9180,  1080, 9560,  2810, 1300,  1070, 3160,  1090, 9620,  1270, 1111,  1201, 9571,  9101, 6360,  1230, 9181,  1390, 1019,  1110, 5972,  1302, 1052,  7380, 9300,  9830, 7360,  1945, 7361,  3220, 9050,  2700, 3300,  2271, 2141,  3280, 1113,  9290, 8230,  2490, 9561,  9331, 3230,  9373, 3250,  9110, 6840,  2120, 9230,  6200, 6213,  6010, 6260,  9440, 3500,  2682, 5971,  8231, 6370,  3022, 5940,  6410, 3210,  2692, 6300,  6930, 9160,  9280, 9480,  6210, 6211, 9404, 9622)
positive <- c(60204, 8185, 60018, 5629, 8501, 4599, 1650, 8180, 60112, 5594, 93052, 2050, 5830, 2070, 1942, 5910, 2345, 8490, 60134, 60267, 8497, 5621, 5991, 5950, 5982, 4598, 93060, 2260, 1340, 60273, 8034, 8117, 7195, 1810, 60262, 1850, 1500, 2550, 60084, 5890, 8220, 8116, 60152, 60288, 2560, 2311, 8600, 1811, 5201, 2540, 60025, 8030, 2040, 8130, 1660, 4601, 60099, 8503, 2240, 60224, 8250, 2209, 8170, 2170, 8041, 5831,  8460, 8300, 60316, 1721, 1610, 1440, 1590, 8280, 1540, 1710, 1604, 2057, 60220, 1460)
# Mixed up 2 and 5.

Mikels_et_al_2005 <- list(pleasant=pleasant, unpleasant=unpleasant)
compare_stimulus_selections(Mikels_et_al_2005)


#
# Study 4
# Attentional rubbernecking: Cognitive control and personality in emotion-induced blindness
# STEVEN B. MOST, MARVIN M. CHUN, and DAVID M. WIDDERS & D. ZALD


negative <- c(1050, 3100, 6313, 2800, 3102, 6350, 3000, 3110, 6560, 3010, 3120, 7361, 3015, 3130, 9040, 3030, 3140, 9253, 3053, 3168, 9405, 3060, 3170, 9410, 3061, 3261, 9433, 3062, 3266, 9570, 3063, 3301, 9571, 3064, 3350, 9800, 3071, 3550, 9810)
neutral <- c(1450, 2480, 2870, 1640, 2485, 2890, 1670, 2487, 4100, 1942, 2495, 4233, 2020, 2500, 4533, 2190, 2515, 4536, 2200, 2516, 4571, 2210, 2520, 4605, 2214, 2560, 4631, 2215, 2570, 5410, 2221, 2575, 7503, 2230, 2580, 7550, 2250, 2590, 8040, 2270, 2600, 8160, 2372, 2620, 8311, 2383, 2702, 9070, 2385, 2749, 9210, 2410, 2840, 9331, 2850)

Most_et_al_2005 <- list(neutral=neutral, negative=negative)
compare_stimulus_selections(Most_et_al_2005)


# Study 5
# Neural Correlates of Using Distancing to Regulate Emotional Responses to Social Situations 
#Harold W. Koenigsberg, M.D. 1,2 , Jin Fan, Ph.D. 1 , Kevin N. Ochsner, Ph.D. 3 , Xun Liu, Ph.D. 1 , Kevin Guise 1 , Scott Pizzarello 4 , Christine Dorantes 1 , Lucia Tecuta , Stephanie Guerreri 1 , Marianne Goodman, M.D. 1,2 , Antonia New, M.D. 1 , Janine Flory, Ph.D. , and Larry J Siever, M.D.

Negative <- c(2053,2490,9810,6821,6360,2691, 6838,3550, 3500, 3181,6311,2710,2900,9433,6313,2800,3160,2095,9910,3530,6312,9800,3350,8485,3180,6315,6200,6242,6510,2455,3301,6230,6530, 6350,9252,3170,3230,9050,6370,2683,2205,3022,6212,9040,6020,6540,3300)

Neutral <- c(5875,4605,2215,2393,2870,4000,8160,2780,2749,2575,2810,2518,2570,2441,9210,2210,2385,2487,2516,8232,2394,8311,2410,9070,2495,8060,5455,2372,2890,2745.1,2580,2635,2702,2880,2514,2850,2493,2440,8010,2485,2499,2830,2235,2980,2383,5410,2515,7550,2020)

Koenigsberg_et_al_2010 <- list(Negative=Negative, Neutral=Neutral)

compare_stimulus_selections(Koenigsberg_et_al_2010)
