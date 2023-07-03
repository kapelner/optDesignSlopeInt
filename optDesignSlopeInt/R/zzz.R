.onAttach = function(libname, pkgname){
  packageStartupMessage(paste(
	  "Welcome to optDesignSlopeInt v", utils::packageVersion("optDesignSlopeInt"), ".\n", 
	  sep = ""
  ))  
}