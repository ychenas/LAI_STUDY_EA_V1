fun_read_nc <- function(arg1, var_st=1) {
# load  ncdf library
  library(ncdf4)
# arg1: filepath fo the nc file from James multilayer output
  print(paste("arg1: for reading file path ;", arg1))
# open the read in file and copy the variables to the dataframe for analysis
  input_nc <- nc_open(arg1)
  result <- list()
  
  for (i in 1:input_nc$ndims ) {
       # store each  variable with respect of the dim <- name
       result[[input_nc$dim[[i]]$name]] <- input_nc$dim[[i]]$vals
    }

  for (i in var_st:length(input_nc$var) ) {
       # store each variable with respect of the var <- name
       result[[input_nc$var[[i]]$name]] <- ncvar_get(input_nc,input_nc$var[[i]]$name)
  }
   nc_close(input_nc)
 # show result structure
 #  print(str(result))
 # export the datatable
   return(result)
  }
