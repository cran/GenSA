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
\Rdversion{0.99.2}
%- Also NEED an '\alias' for EACH other topic documented here.
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
GenSA(x, param, lb, ub, fn, jc=NULL, control=list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{x}{
		Vector. Initial values for the components to be optimized. Default is \code{NULL}, in which case, default values will be generated automatically.
	}
	\item{param}{
		Vector. Parameters that can be used in the function that is being optimized. If parameters are not needed, it has to be set to \code{NULL}.
		Please refer to the Schubert function example for details.
	}
	\item{lb}{
		Vector with length of \code{x}. Lower bounds for each component.
	}
	\item{ub}{
		Vector. with length of \code{x}. Upper bounds for each component.
	}
	\item{fn}{
		Objective function which the user would like to optimize. The objective function has two arguments, \code{x} and \code{param}. \code{fn} will return a scalar result. Please refer to the given example.
	}
	\item{jc}{
		Constraint function. A function that describes the linear/non-linear constraints with first argument \code{x} and second argument \code{param}.
        \code{jc} has to return a logical value. Default value is \code{NULL} which means there is not any further constraint except \code{lb} and \code{ub}.
        Please refer to the given example.
	}
	\item{control}{
		The argument is a list that can be used to control the behavior of the algorithm for advanced usage, the following controls are available:
		\describe{
			\item{rseed:}{
				Integer used by the \code{set.seed} function for random number generation in R in case initial values of \code{x} are not
				provided.
			}
			\item{max.step:}{
				Integer. Maximum number of run of the algorithm. Default value is 5000.
			}
			\item{interval:}{
				Integer. Number of results output per interval steps when \code{verbose} control option is set to \code{TRUE}.
			}
			\item{know.real.energy:}{
				Logical. \code{TRUE} means that we know the real energy in advance which can stops the algorithm earlier after the it has found the real
				global minimum. This control option can be used for testing purpose. Default value is \code{FALSE} which means we don't know the real energy value.
			}
			\item{error.real.energy:}{
				In case \code{know.real.energy} is \code{TRUE}, percentage value corresponding to the accepted error. The algorithm will stop when:\cr
			 	abs(\code{current.energy-real.energy})*100/abs(\code{real.energy}) <= \code{error.real.energy} is satisfied.
				Default value is 0.01.
			}
			\item{real.energy:}{
				Numeric. In case \code{know.real.energy} is \code{TRUE}, it corresponds to the real energy value (real global minimum of the objective function).
			}
			\item{has.judge.function}{
				Logical. \code{TRUE} means that the constraint function passed as \code{jc} argument will be called by the algorithm. Default is \code{FALSE}.
			}
			\item{temp.init:}{
				Numeric. Initial value for temperature.	Default value is 3000.
			}
			\item{visiting.param:}{
				Numeric. Parameter for visiting distribution. Default value is 2.62.
			}
			\item{acceptance.param:}{
				Numeric. Parameter for acceptance distribution. Default value is -5.0.
			}
			\item{component.change:}{
				Numeric. 1 means that components of \code{x} will be changed one by one, 2 means that they will be changed simultaneously. Default value is 2.
			}
			\item{seed.init:}{
				Negative integer used for internal behavior of the algorithm.  Default value is -1333333. In general, this should not be modified.
			}
			\item{seed.random:}{
				Negative integer used for internal behavior of the algorithm.  Default value is -100377. In general, this should not be modified.
			}
			\item{markov.length:}{
				Integer. Length of the Markov chain. Default value is the length of \code{x}.
			}
			\item{verbose:}{
				Logical. \code{TRUE} means that messages from the algorithm are shown. Default is \code{FALSE}.
				Note that displaying algorithm messages on the console decreases a lot the performances. Default option is recommended.
			}
		}
	}
}

	
\value{
	The returned value is a list with the following fields:
	\describe{
		\item{energy:}{
			Numeric. Global minimum identified of the objective function.
		}
		\item{x:}{
			Coordinates of the global minimum found.
		}
		\item{trace.mat:}{
			A matrix which contains the history of the algorithm. (By columns: Record index, number of steps, temperature, number of function calls,
			current energy value, current minimum energy value). Each row corresponds to a run record.
		}
		\item{nb.calls:}{
			Integer. Total number of calls of the objective function.
		}
	}
}

\keyword{ package, Generalized Simulated Annealing, optimization }
\seealso{
optim
}

\references{
Xiang Y, Gubian S, Suomela B, Hoeng (2011). "Generalized Simulated Annealing for Efficient Global Optimization: the GenSA Package for R". Journal of Statistical Software, Forthcoming. \url{http://www.jstatsoft.org/}.

Tsallis C (1988). "Possible generalization of Boltzmann-Gibbs statistics." Journal of Statistical Physics, 52, 479--487.

Tsallis C, Stariolo DA (1996). "Generalized Simulated Annealing." Physica A, 233, 395--406.
 	
Xiang Y, Sun DY, Fan W, Gong XG (1997). "Generalized Simulated Annealing Algorithm and Its Application to the Thomson Model." Physics Letters A, 233, 216--220.

Xiang Y, Gong XG (2000a). "Efficiency of Generalized Simulated Annealing." PHYSICAL REVIEW E, 62, 4473.

Xiang Y, Sun DY, Gong XG (2000). "Generalized Simulated Annealing Studies on Structures and Properties of Nin (n=2-55) Clusters." Journal of Physical Chemistry A, 104, 2746--2751.
}


\examples{
library("GenSA")
# Example of Rosenbrock banana function with param = 2
# The global minimum is 0 at c(1,1).
Banana <- function(x,param) {   
    x1 <- x[1]
    x2 <- x[2]
    exponential <- param 
    100 * (x2 - x1 * x1)^exponential + (1 - x1)^exponential
}

# Identify the global minimum of Rosenbrock banana function by using GenSA.
# The param of Rosenbrock banana function is 2.
lb <- c(-2,-1); ub <- c(2,3)
result <- GenSA(x=NULL, param=2, lb=lb, ub=ub, fn=Banana)
cat("The global minimum of Rosenbrock banana function is: ",result$energy,
"\n",fill=TRUE)
cat("The coordinates of global minimum of Rosenbrock banana function are: ",
 result$x,"\n",fill=TRUE)
 
# Sometimes there is a linear or nonlinear constraint for the objective function. 
# Here we add a nonlinear constraint which is an annulus surrounding c(1,1), 
# the global minimum of Banana function. Optimization algorithm can only search 
# the area which is not covered by this annulus constraint. 

# Example of constraint function
judge_constraint <- function(x,param) {
  # judge function  
  # Below give a annulus constraint, which surround c(1,1), the global 
  # minimum of Banana function.
  in_constraint <- ( (x[1] -1)^2 + (x[2] -1)^2 <= 0.5^2 ) || 
   ( (x[1] -1)^2 + (x[2] -1)^2 >= 0.6^2 )
  return(in_constraint)
}

# Identify the global minimum of Rosenbrock Banana function with a non-linear
# constraint by using GenSA.
result <- GenSA(x=NULL, param=2, lb=lb, ub=ub, fn=Banana,jc=judge_constraint)
cat("The global minimum of Rosenbrock Banana function is: ",result$energy,
"\n",fill=TRUE)
cat("The coordinates of global minimum of Rosenbrock Banana function are:",
 result$x,"\n",sep=" ",fill=TRUE)


# Example of function Shubert (SH) with param=NULL since this function
# has no parameter. Please note: though there is no parameter, the function
# should still have two input variables as shown below. This function 
# has 760 local minima, 18 of which are global minima with function 
# value -186.7309. -10 <= x$_{i}$ <=10
Shubert <- function(x,param) {      
  x1 <- x[1]
  x2 <- x[2]
  sum1 <- 0 
  for (i in 1:5) {
    tmp <- i*cos((i+1)*x1+i) 
    sum1 <- sum1+tmp
  }
  sum2 <- 0 
  for (i in 1:5) {
    tmp <- i*cos((i+1)*x2+i) 
    sum2 <- sum2+tmp
  }
  sum1*sum2
}

# Identify the global minimum of Shubert function by using GenSA.
lb <- rep(-10,2)
ub <- rep(10,2)
result <- GenSA(x=NULL, param=NULL, lb, ub, fn=Shubert)
cat("The global minimum of Shubert function is: ",result$energy,"\n")
cat("The coordinates of global minimum of Shubert function are: ",
 result$x,"\n",fill=TRUE)
}

\author{
Yang Xiang, Sylvain Gubian, Brian Suomela, Julia Hoeng, PMP SA.
}


