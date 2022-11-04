# load geomorph package library
# also load all functions in the wingmorphofunctions.R script
# (that needs to be done using that script)
library(geomorph)

###############################################
##### Create 3d array of 2 original triangles
###############################################
# The 3d array is the standard input to geomorph
# Create two vectors
T2vec <- c(1.55, 3.10, 1.07, 0.82, -0.72, -1.65)
T1vec <- c(0,1,-1,1,-1,-1)
# pass these vectors as input to the array.
# to make 3d array input
result <- array(c(T1vec, T2vec), dim = c(3,2,2))
# subset to the two original triangles for future use
T1<-result[,,1]
T2<-result[,,2]
plot(x=T1[,1], y = T1[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T2[,1], y=T2[,2], col=c('black', 'blue', 'gray'))

#######################################
##### Generalized Procrustes Analysis
#######################################
gpatri<-gpagen(result)
plot(gpatri)
# subset results to two aligned triangles for future use
T1post<-gpatri$coords[,,1]
T2post<-gpatri$coords[,,2]
plot(x=T1post[,1], y = T1post[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T2post[,1], y=T2post[,2], col=c('black', 'blue', 'gray'))

# Focus on T2: plot pre and post
plot(x=T2[,1], y = T2[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T2post[,1], y=T2post[,2], col=c('black', 'blue', 'gray'))

##############################################################################
##### Find centroids of both original triangles and the target post triangle
##############################################################################
# will be helpful later
# see "findCentroid2d" function in wingmorphofunctions.R
T1cent <- findCentroid2d(T1)
T2cent <- findCentroid2d(T2)
T1postcent <- findCentroid2d(T1post)
T2postcent<-findCentroid2d(T2post)

plot(x=T1post[,1], y = T1post[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T1postcent[1], y=T1postcent[2], col=c('pink'))

############################################
##### Find translation between old and new
############################################
# see "findTransl" and "Translate" functions in wingmorphofunctions.R
transl <- findTransl(T2, T2post)

# check by applying translation
T2check<-Translate(T2, transl)
plot(x=T2[,1], y = T2[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T2post[,1], y=T2post[,2], col=c('black', 'blue', 'gray'))
points(x=T2check[,1], y=T2check[,2], col=c('pink', 'purple', 'green'))
# the centroid of 'T2check' should now be the same as 'T2post' =~0
# because we shifted T2 to lie at T2post
T2checkcent <- findCentroid2d(T2check)
# here the accumulation of floating point error is evident

# note the order of the point sets in the functions:
# the first set of points in the transl function 
# the set that gets moved by translate
# works if you reverse and want to move the other
# triangle, too:
posttransl <- findTransl(T2post, T2)
T2postcheck <- Translate(T2post, posttransl)
plot(x=T2[,1], y = T2[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
points(x=T2post[,1], y=T2post[,2], col=c('black', 'blue', 'gray'))
points(x=T2postcheck[,1], y=T2postcheck[,2], col=c('black', 'blue', 'gray'))
# Centroids should be the same
# 'T2postcheckcent' =~ 'T2cent'
# because we shifted T2post to lie at T2
T2postcheckcent<-findCentroid2d(T2postcheck)

#######################################
##### Find Centroid Sizes for scaling
#######################################
# see "findCentroidSize" function in wingmorphofunctions.R
T1cs <- findCentroidSize(T1, T1cent)
T2cs <- findCentroidSize(T2, T2cent)
T2postcs<-findCentroidSize(T2post,T2postcent)
T1postcs<-findCentroidSize(T1post,T1postcent)
# T1postcs and T2postcs should both be =~1 as
# the GPA has rescaled the point clouds to unit centroid size.
# Floating point error (maybe some from GPA algorithm) is here.
# so not exactly one. 

#######################################################################
##### Find the rotation angle between original and post point clouds
#######################################################################
# see "findAng" function in wingmorphofunctions.R

T2angle <- FindAng(T2, T2post, transl, T2cs)

T2angle

# check by applying the translation, scaling, rotation
TransformPoints<- function(points, trans, size, ang){
  # Translate
  new <- Translate(points, trans)
  # Scale
  newsize <- new/size
  # Rotate
  mark <- length(points[,1])
  x = 1
  newang <- data.frame()
  while (x <= mark){
    oldx <- newsize[x, 1]
    oldy <- newsize[x, 2]
    newx <- oldx*cos(-ang) - oldy*sin(-ang)
    newy <- oldx*sin(-ang) + oldy*cos(-ang)
    coords <- c(newx, newy)
    newang <- rbind(newang, coords)
    x = x+1
  }
  return(newang)
}

T2trans <- TransformPoints(T2, transl, T2cs, T2angle)

#plot post-gpa landmarks
plot(x=T2post[,1], y = T2post[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
#plot pre-gpa landmarks
points(x=T2[,1], y = T2[,2], col=c('black','blue','gray'))
#plot post-gpa-transformed extra points (sensors)
points(x=T2trans[,1], y=T2trans[,2], col=c('pink','purple', 'green'))

# Try for T1 and T1post
translT1<-findTransl(T1, T1post)
T1angle<-FindAng(T1, T1post, translT1, T1cs)
T1trans<-TransformPoints(T1, translT1, T1cs, T1angle)

#plot post-gpa landmarks
plot(x=T1post[,1], y = T1post[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
#plot pre-gpa landmarks
points(x=T1[,1], y = T1[,2], col=c('black','blue','gray'))
#plot post-gpa-transformed extra points (sensors)
points(x=T1trans[,1], y=T1trans[,2], col=c('pink','purple', 'green'))

#####################################
###### Test on non-landmark points
#####################################
# Points are associated with T2 but are not used in the GPA.
# Now, we want to take the transforms applied to T2 in the GPA 
# and apply them to these non-landmark points so that they are 
# moved to the same relative locations on T2post. 

pointsvec <- c(2, 2, 0.5, -0.5)
points <- array(pointsvec, dim = c(2,2,1))[,,1]
#print(points)

new <- TransformPoints(points, transl, T2cs, T2angle)

#plot post-gpa landmarks
plot(x=T2post[,1], y = T2post[,2], xlim=c(-2,4), ylim=c(-2,2), col=c('black', 'blue', 'gray'))
#plot pre-gpa landmarks
points(x=T2[,1], y = T2[,2], col=c('black','blue','gray'))
#plot pre-gpa extra points (sensors)
points(x=points[,1], y=points[,2], col=c('green','purple'))
#plot post-gpa-transformed extra points (sensors)
points(x=new[,1], y=new[,2], col=c('green','purple'))

# Check to see if relative distances are preseved, as they should be if "shape"
# is preserved
grTblT2post <- sqrt((new[1,1]-T2post[1,1])^2 + (new[1,2]-T2post[1,2])^2) 
grTblT2 <- sqrt((points[1,1]-T2[1,1])^2 + (points[1,2]-T2[1,2])^2) 
grTpurpT2post <- sqrt((new[1,1]-new[2,1])^2 + (new[1,2]-new[2,2])^2) 
grTpurpT2 <- sqrt((points[1,1]-points[2,1])^2 + (points[1,2]-points[2,2])^2)
# these should be the same
T2postratio <- grTblT2post/grTpurpT2post
T2ratio <- grTblT2/grTpurpT2
