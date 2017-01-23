#' Clear working session
#' 
#' Function to transform P-values into symbols of significance (***)
#' @param s A vector of P-values
#' @return TRUE for clear
#' 
#' @author C.J.
#' @export

clearSession <- function(){
  lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
}
