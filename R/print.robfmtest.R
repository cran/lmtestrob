print.robfmtest <-
function(x, digits = 4, ...)
{
  cat(x$method)
  cat(" \n")
  cat("data: ")
  cat(x$data.name)
  cat(" \n")
  cat("FMrob = ")
  cat(drop(x$statistic), digits = digits)
  cat(", Degrees of freedom = ")
  cat(drop(x$dof))
  cat(", p-value = ")
  cat(x$p.value)
}
