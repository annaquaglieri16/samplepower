# Return true only if the object exists and has a specific class check existence of an object is


exists_with_class <- function(x, check_class,silent=FALSE){

  # check existence
  exist <- exists(x = x) # t/F

  if(exist){

    try_class <- class(get0(x)) == check_class

    if(exist & try_class){
        return(TRUE)
    }else{

      if(!silent){
        message(paste0("Object ",x, " exists but is not of class ",check_class))
      }
      return(FALSE)
    }
  }else{
    if(!silent){
      message(paste0("Object ",x, " does not exist in your R environment"))
    }
    return(FALSE)
  }

}
