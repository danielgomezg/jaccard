#define Jaccard Similarity function 
jaccardManual <- function (a, b) {
  intersección <- length( intersect(a, b))
  union <- length(a) + length (b) 
    return (intersección / union)
}