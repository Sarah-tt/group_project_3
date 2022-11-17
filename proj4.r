#Qin Wang :s2304682
#Yanting Wang :s2327787
#Yanyan Qiu :s2261604

#address of github repo: https://github.com/Sarah-tt/group_project_3.git
#Qin Wang finish the step 1,2
#Yanting Wang finish the step 3,4
#Yanyan Qiu finish the step 5,6

#All of the students contribute the same proportion for this project.

##This function use Newton's method for minimization of functions.
#In summary, starting with n = 0(n is the number of iterate), and a guesstimate theta[0],iterate these steps:
#1. We need to check whether the objective or derivatives are finite at initial theta[0], if not the function will stop.
#2. Test whether theta[n] is a minimum value, if it is, then stop iterate and find the minimum value.
#3. Test whether Hessian matrix is positive definite using the eigen value larger than 0, if not transfer it to a positive definite matrix.
#4. Solve hessian matrix*direction = -gradient for the search direction
#5. If f(theta + d) is not < f(theta), repeatedly halve d until it is. If it's more than a feasible times to halve, stop function.
#6. Set theta[n+1] = theta[n] + d, n = n + 1 and return to step 2
#If the number of iterations is exceeded the function will stop. 



# The newt function is to implement Newton’s method for minimization of functions.
# The input of the function is as follows：
# (1) theta: a vector of initial values for the guesstimate
# (2) func: the objective function to minimize. Its first argument is the vector of guesstimate. 
#     Remaining arguments will be passed from newt using ‘...’.
# (3) grad: the gradient function. It has the same arguments as func but returns the gradient vector of the
#     objective w.r.t. the elements of parameter vector.
# (4) hess: the Hessian matrix function.  It has the same arguments as func but returns the Hessian matrix of the 
#     objective w.r.t. the elements of parameter vector.
# (5) tol: the convergence tolerance.
# (6) fscale: a rough estimate of the magnitude of func near the optimum - used in convergence testing.
# (7) maxit: the maximum number of Newton iterations to try before giving up.
# (8) max.half: the maximum number of times a step should be halved before concluding that the step has failed to improve the objective.
# (9) eps: the finite difference intervals to use when a Hessian function is not provided.

# And the this function return a list containing:
# f the value of the objective function at the minimum.
# theta the value of the parameters at the minimum.
# iter the number of iterations taken to reach the minimum.
# g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
# Hi the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).

newt = function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  #vector of initial values for the guesstimate
  x = theta
  #set up number of iterations to 1
  n = 1
  #check whether objective is not finite at the initial theta
  #if objective is infinite at the initial theta, the function will issue error
  if (func(x,...) == Inf | func(x,...) == -Inf)
    stop("The function has no domain at the current value.")
  #check whether derivatives is not finite at the initial theta
  #for each element in gradient vector, if it is infinite
  #then the function will issue error
  for(k in 1:length(grad(x,...))){
    if (grad(x,...)[k] == Inf | grad(x,...)[k] == -Inf)
      stop("The gradient function has no derivative at the current value.")
  }
  
  #we can iterate up to maxit numbers in order to reach the minimum
  for (n in 1:maxit){
    #set up gradient vector at vector x
    g = grad(x,...)
    #check whether the user enters the Hessian matrix, 
    #if not supplied, we will get the Hessian matrix by finite differencing of the gradient vector
    if(missing(hess)){ #if the Hessian matrix is not supplied
      #count the length of gradient vector
      l = length(grad(x,...))
      #set up a empty matrix to store the element in Hessian matrix
      #it has l row and l column
      hs = matrix(nrow = l,ncol = l)
      #get the gradient vector at the elements of parameter vector
      derivative = grad(x,...)
      #for i row in Hessian matrix, for j column in Hessian matrix
      for (i in 1:l){
        for (j in 1:l){
          #obtain an approximation to the Hessian by finite differencing of the gradient vector
          #the difference interval is eps, which is a small number close to 0
          adjust = x
          adjust[j] = adjust[j]+eps
          adjust_derivative = grad(adjust,...)
          hs[i,j] = (adjust_derivative[i]-derivative[i])/eps
        }
      }
    }
    #if the Hessian matrix is provided, we can use it directly
    else{ 
      hs = hess(x,...)
    }