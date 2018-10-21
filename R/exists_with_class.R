# Return true only if the object exists and has a specific class check existence of an object is


exists_with_class <- function(obj, check_class,silent=FALSE,where){

  # check existence
  exist <- exists(obj,where = where,inherits=TRUE) # t/F

  if(exist){

    try_class <- class(get0(obj,envir=where)) == check_class

    if(exist & try_class){
        return(TRUE)
    }else{

      if(!silent){
        message(paste0("Object ",obj, " exists but is not of class ",check_class))
      }
      return(FALSE)
    }
  }else{
    if(!silent){
      message(paste0("Object ",obj, " does not exist in your R environment"))
    }
    return(FALSE)
  }

}
