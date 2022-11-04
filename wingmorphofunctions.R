# finds centroid for a 2d set of points
# a helper function for findTransl
findCentroid2d <- function(coords) {
  x = sum(coords[,1])/length(coords[,1])
  y = sum(coords[,2])/length(coords[,2])
  return(c(x,y))
}

# finds the translation vectors between an original 
# and an aligned set of points
findTransl <- function(coords1, coords2) {
  cent1 <- findCentroid2d(coords1)
  cent2 <- findCentroid2d(coords2)
  #map coords1 onto coords2
  newx <- cent1[1] - cent2[1]
  newy <- cent1[2] - cent2[2]
  coords <-cbind(newx, newy)
  return(coords)
}

# Translates coordinates using the result from findTransl
Translate <- function(coords, transl){
  # first translate and scale
  mark <- length(coords[,1])
  x = 1
  new <- data.frame()
  while (x <= mark) {
    transx <- coords[x, 1] - transl[1]
    transy <- coords[x, 2] - transl[2]
    x = x + 1
    coord <- c(transx, transy)
    new<-rbind(new, coord)
  }
  return(new)
}

# finds the centroid size of a set of points
findCentroidSize <- function(list, centroid) {
  n = length(list[,1])
  centsizelist = c()
  while (n>0) {
    x = list[n,]
    dist = sqrt((x[1] - centroid[1])^2 + (x[2] - centroid[2])^2)
    distsq = dist^2
    centsizelist = append(centsizelist, distsq[[1]])
    n = n-1
  }
  CS = sqrt(sum(centsizelist))
  return(CS)
}

# finds the rotation angle between an original and the 
# post-aligned set of points
FindAng <- function(coords1, coords2, transl, cs){
  new <- Translate(coords1, transl)
  #scale
  newsized <- new/cs
  # doing this for one coordinate pair can yield inexact results, 
  # here we take the mean of all coordinates of all points
  mark <- length(coords1[,1])
  x = 1
  ang<-c()
  while (x <= mark) {
    Xold <- coords2[x,1]
    Yold <- coords2[x,2]
    Xnew <- newsized[x,1]
    Ynew <- newsized[x,2]
    cosang1 <- ((Yold*Ynew+Xold*Xnew)/(Xold^2+Yold^2))
    angle<-acos(cosang1)
    ang <- append(ang, angle)
    x = x +1
  }
  return(mean(ang))
}
