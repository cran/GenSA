# File: GenSA.R
# 
# Author: Sylvain Gubian
# Aim: Function for General Simulated Annealing

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

#########################################################################################


GenSA <- function(x, param, lb, ub, fn, jc=NULL, control=list())
{
	# Do some checks
	if (!is.function(fn) || is.null(fn)) {
		stop("'fn' has to be a R function.")
	}
	
	# Create an environment for sharing data between R and C
	genSA.env <- new.env(hash=TRUE, parent=emptyenv())
	assign("params", param, envir=genSA.env)
	
	fn1 <- function(x) {
		ret <- try(fn(x,get("params", genSA.env)))
		if (is.character(ret))
			return(10e9)	
		else
			return(ret)
	}
	
	if (!is.null(jc)) {
		if (!is.function(jc)) {
			stop("'jc' has to be a R function.")
		}
		jc1 <-  function(x) {
			return(jc(x, get("params", genSA.env)))
		}
	}
	else
	{
		jc1 <- function(x) {
			return(TRUE)
		}
	}
	
	con <- list(
			rseed = 1000,
			seed.init = -1333333,
			seed.random = -100377,
			max.step = 5000,
			interval = 1,
			know.real.energy = FALSE,
			error.real.energy = 0.01,
			real.energy = 0.0,
			has.judge.function = FALSE,
			temp.init = 3000,
			visiting.param = 2.62,
			acceptance.param = -5.0,
			component.change = 2,
			markov.length = length(lb),
			verbose = FALSE
	)
	
	
	# Perform some checks before callinc C code
	nmsC <- names(con)
	con[(namc <- names(control))] <- control
	if(length(noNms <- namc[!namc %in% nmsC]))
		warning("unknown names in control: ", paste(noNms,collapse=", "))
	
	if (is.null(jc)) {
		con$has.judge.function = FALSE ;
	}
	else con$has.judge.function = TRUE ;
	
	if (length(x)==0 && (length(lb)==0 || length(ub)==0)) {
		stop("There is no x or no lower/upper bounds defined")
	}
	
	if (length(lb) != length(ub)) {
		stop("Lower and upper bounds vector do not have the same length")
	}
	
	if (length(x) > 0 ) {
			
		if (length(lb)==0 || length(lb) != length(x)) {
			if (con$verbose) {
				message("Lower bounds vector size does not match with x size, using -Inf")
			}
			lb <- rep(-Inf, length(x))
		}
		if (length(ub)==0 || length(ub) != length(x)) {
			if (con$verbose) {
				message("Upperr bounds vector size does not match with x size, using -Inf")
			}
			ub <- rep(Inf, length(x))
		}
	}
	
	else {
		if (con$verbose) {
			message("Initializing x with random data inside bounds")
		}
		x <- vector()
		#initialize x with random values in the bounds
		set.seed(con$rseed)
		x <- runif(length(lb), min=min(lb), max=max(ub)) 
	}
	
	if (is.null(param))
		param <- 0
	
	ret <- list()
	# Create instance of the GenSACaller
	instance <- .Call(createInstance)
	if (is.null(instance)) {
		stop("Can not create GenSACaller instance!")
	}
	
	# Call execute on the instance
	res <- .Call(execute, x, param, lb, ub, fn1, jc1, con, genSA.env, instance)
	if (is.null(res)) {
		stop("Can not call execute function on instance")
	}
	
	# Get the results in a list
	res <- .Call(getREnergy, instance)
	if (is.null(res)) {
		message("Can not get energy value")
	}
	else
	{
		ret$energy <- res
	}
	
	res <- .Call(getRXMiniVector, instance)
	if (is.null(res)) {
		message("Can not get calculated x values")
	}
	else
	{
		ret$x <- res
	}
	
	res <- .Call(getRTraceMat, instance)
	if (is.null(res)) {
		message("Can not get Matrice of traces")
	}
	else
	{
		ret$trace.mat <- t(as.matrix(res))
		ret$trace.mat <- ret$trace.mat[1:which(ret$trace.mat[,1]==0)[2]-1,,drop=FALSE]
		#ret$trace.mat <- ret$trace.mat[1:(con$markov.length * 6 * con$max.step),]
		colnames(ret$trace.mat) <- c("record.index", "nb.steps", "temperature", "nb.calls", "current.energy", "current.minimum")
	}
	
	res <- .Call(getRNbFuncCall, instance)
	if (is.null(res)) {
		message("Can not get number of function calls")
	}
	else
	{
		ret$nb.calls <- res
	}
	.Call(releaseInstance, instance)
	return(ret)
}


