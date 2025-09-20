#APU functions ----

myrecoding <- function(vect, recoding) {
  vect <- factor(vect)
  levs <- setNames(levels(vect), levels(vect))
  stopifnot(all(recoding %in% levs))
  levs[recoding] <- names(recoding)
  levels(vect) <- levs
  return(as.character(vect))
}
