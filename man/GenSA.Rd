%- File: GenSA.Rd
%- 
%- Author: sylvain.gubian@pmintl.com
%-
%- This program is distributed in the hope that it will be useful,
%- but WITHOUT ANY WARRANTY; without even the implied warranty of
%- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%- GNU General Public License for more details.
%-
%-#########################################################################################

\name{GenSA}
\Rdversion{1.0.0}
\alias{GenSA}
\title{
Generalized Simulated Annealing Function
}
\description{
This function searches for global minimum of a very complex non-linear objective function with a very large number of optima. 
Another attractive capability of  this function is that it can automatically handle all types of 
constraints, ranging from the simplest box-constraints to even the most complex non-linear constraints.
}

\note{
In this package, the value of the objective function is called energy.
}
\usage{
GenSA(par, lower, upper, fn, control=list(), ...)
}

\details{
	The default values of the control components are set for a complex optimization problem. For usual optimization problem with medium
	complexity, \code{GenSA} can find a reasonable solution quickly so the user is strongly recommended to let \code{GenSA} stop
	earlier by setting \code{threshold.stop} if \code{threshold.stop} is the targeted function value, or by setting \code{max.time}
	if the user just want to run \code{GenSA} \code{max.time} seconds, or by setting \code{max.call} if the user just want to run 
	\code{GenSA} within \code{max.call} function calls. Please refer to the examples below. For very complex function, the user is recommended to increase \code{maxit} and 
	\code{temp}. 
	   
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{par}{
		Vector. Initial values for the components to be optimized. Default is \code{NULL}, in which case, default values will be generated automatically.
	}
	\item{lower}{
		Vector with length of \code{par}. Lower bounds for components.
	}
	\item{upper}{
		Vector with length of \code{par}. Upper bounds for components.
	}
	\item{fn}{
		A function to be minimized, with first argument the vector of parameters over which minimization is to take place. It should return
		a scalar result.
	}
	\item{...}{
		allows the user to pass additional arguments to the function \code{fn}.
	}
	\item{control}{
		The argument is a list that can be used to control the behavior of the algorithm:
		\describe{
			\item{maxit:}{
				Integer. Maximum number of iterations of the algorithm. 
			}
			\item{threshold.stop:}{
				Numeric. The program will stop when the targeted objective function value \code{threshold.stop} is reached. Default value is \code{NULL}
			}
			\item{nb.stop.improvement:}{
				Integer. Number of iteration to let run when the algorithm does not give any improvement in the minimum value found of the objective function.
			}
			\item{smooth}{
				Logical.\code{TRUE} when the objective function is smooth, or differentiable almost everywhere in the
				region of \code{par}, \code{FALSE} otherwise. Default value is \code{TRUE}.
			}
			\item{max.call}{
				Integer. Maximum number of call of the objective function. Default is set to 1e7.
			}
			\item{max.time}{
				Numeric. Maximum running time in seconds.
			}
			\item{temp:}{
				Numeric. Initial value for temperature.
			}
			\item{visiting.param:}{
				Numeric. Parameter for visiting distribution.
			}
			\item{acceptance.param:}{
				Numeric. Parameter for acceptance distribution.
			}
			\item{verbose:}{
				Logical. \code{TRUE} means that messages from the algorithm are shown. Default is \code{FALSE}.
			}
		}
	}
}

	
\value{
	The returned value is a list with the following fields:
	\describe{
		\item{value:}{
			Numeric. The value of \code{fn} corresponding to \code{par}.
		}
		\item{par:}{
			Vector. The best set of parameters found.
		}
		\item{trace.mat:}{
			A matrix which contains the history of the algorithm. (By columns: Step number, temperature,
			current objective function value, current minimum objective function value).
		}
		\item{counts:}{
			Integer. Total number of calls of the objective function.
		}
	}
}

\keyword{ package, Generalized Simulated Annealing, optimization }
\seealso{
optim
}

\references{
Xiang Y, Gubian S, Suomela B, Hoeng (2011). "Generalized Simulated Annealing for Efficient Global Optimization: the GenSA Package for R". The R Journal, Forthcoming. \url{http://journal.r-project.org/}.

Tsallis C (1988). "Possible generalization of Boltzmann-Gibbs statistics." Journal of Statistical Physics, 52, 479--487.

Tsallis C, Stariolo DA (1996). "Generalized Simulated Annealing." Physica A, 233, 395--406.
 	
Xiang Y, Sun DY, Fan W, Gong XG (1997). "Generalized Simulated Annealing Algorithm and Its Application to the Thomson Model." Physics Letters A, 233, 216--220.

Xiang Y, Gong XG (2000a). "Efficiency of Generalized Simulated Annealing." PHYSICAL REVIEW E, 62, 4473.

Xiang Y, Sun DY, Gong XG (2000). "Generalized Simulated Annealing Studies on Structures and Properties of Nin (n=2-55) Clusters." Journal of Physical Chemistry A, 104, 2746--2751.
}


\examples{
library(GenSA)
# Try Rastrgin function (The objective function value for global minimum is 0 with all components of par are 0.)
rastrigin <- function(x) {
  return(sum(sapply(1:(length(x)), function(i){ x[i]^2 - 10*cos(2*pi*x[i]) + 10 })))  
}
# Perform the search on a 30 dimensions rastrigin function. Rastrigin function with dimension 30 is known as the most difficult optimization
# problem according to
# GenSA will stop after finding the targeted function value 0 with absolute tolerance 1e-13
set.seed(1234) # The user can use any seed.
dimension <- 30
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = rastrigin, control=list(threshold.stop=global.min+tol,verbose=TRUE))
out[c("value","par","counts")]

# GenSA will stop after running during 10 seconds
set.seed(1234) # The user can use any seed.
dimension <- 30
global.min <- 0
tol <- 1e-13
lower <- rep(-5.12, dimension)
upper <- rep(5.12, dimension)
out <- GenSA(lower = lower, upper = upper, fn = rastrigin, control=list(max.time=10))
out[c("value","par","counts")]
}



\author{
Yang Xiang, Sylvain Gubian, Brian Suomela, Julia Hoeng, PMP SA.
}


