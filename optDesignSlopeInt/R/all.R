

#' Plots a standard error estimate of thetahat (slope over intercept) over a range of possible theta0 values
#' in order to investigate robustness of the the initial theta0 guess.
#' 
#' @param n 			The number of experimental runs.
#' @param xmin 			The minimum value of the independent variable.
#' @param xmax 			The maximum value of the independent variable.
#' @param theta			The putative true value. This is used to see how much efficiency given up by designing it for \code{theta0}.
#' @param theta0 		The guess used to construct the experimental design. Specify only if you wish to see this 
#' 						value plotted. Default is \code{NULL}.
#' @param theta0_min 	Simulating over different guesses of theta0, this is the minimum guess.
#' @param theta0_max 	Simulating over different guesses of theta0, this is the maximum guess.
#' @param sigsq 		A guess to be used for the homoskedastic variance of the measurement errors. If known accurately,
#' 						then the standard errors (i.e. the y-axis on the plot) will be accurate. Otherwise, the standard
#' 						errors are useful only when compared to each other in a relative sense. Defaults to \code{1}.
#' @param RES 			The number of points on the x-axis to simulate. Higher numbers will give smoother results. Default is \code{20}.
#' @param Nsim 			The number of models to be simulated for estimating the standard error at each value on the x-axis. Default is \code{1000}.
#' @param theta_logged	Should the values of theta be logged? Default is \code{TRUE}.
#' 
#' @return 				A list with original parameters as well as data from the simulation
#' 
#' @author 				Adam Kapelner
#' @export
std_err_vs_theta0_plot_for_homo_design = function(n, xmin, xmax, theta, theta0_min, theta0_max, theta0 = NULL, beta0 = 1, sigsq = 1, RES = 20, Nsim = 1000, theta_logged = TRUE){
	theta0s = seq(theta0_min, theta0_max, length.out = RES)

	homo_designs = matrix(NA, nrow = 0, ncol = n)
	for (i in 1 : RES){		
		homo_designs = rbind(homo_designs, oed_for_slope_over_intercept_homo(n, xmin, xmax, theta0s[i]))
	}
	
	std_errors = array(NA, RES)
	for (d in 1 : RES){
		ests_opt = array(NA, Nsim)
		for (nsim in 1 : Nsim){
			xs_opt = homo_designs[d, ]			
			ys_opt = beta0 + theta  * beta0 * xs_opt + rnorm(n, 0, sqrt(sigsq))			
			mod_opt = lm(ys_opt ~ xs_opt)			
			ests_opt[nsim] = coef(mod_opt)[2] / coef(mod_opt)[1]
		}
		std_errors[i] = sd(ests_opt)
	}	
	
	plot(theta0s, std_errors, type = "l", col = "grey", ylab = "standard error of theta-hat", xlab = "theta0 hypothesis")
	opt_design = oed_for_slope_over_intercept_homo(n, xmin, xmax, theta0)
	true_rhostar = sum(opt_design == xmin) / n
	abline(v = true_rhostar, col = "blue")

	lines(theta0s, predict(loess(std_errors ~ theta0s)), col = "black", lwd = 2)
	#return the simulated results only if user cares
	invisible(list(
		n = n,
		xmin = xmin,
		xmax = xmax,
		theta = theta, 
		res = cbind(theta0s, std_errors)
	))		
}

#' Create an optimal design for measuring the slope divided by the intercept
#' 
#' @param n 			The number of experimental runs.
#' @param xmin 			The minimum value of the independent variable.
#' @param xmax 			The maximum value of the independent variable.
#' @param theta0 		The guess of the true value of the slope / intercept.
#' @param f_hetero		Specification of heteroskedasticity: the h(x) which relates the value of the 
#' 						independent variable to the variance in the response around the line at that place
#' 						or the proportional variance at that point. If \code{NULL}, homoskedasticity is
#' 						assumed (this is the default behavior).
#' 
#' @return 				An n-vector of x-values which specifies the optimal design
#' 
#' @author 				Adam Kapelner
#' @export
oed_for_slope_over_intercept = function(n, xmin, xmax, theta0, f_hetero = NULL, ...){
	#throw errors here
	if (is.null(f_hetero)){
		oed_for_slope_over_intercept_homo(n, xmin, xmax, theta0)
	} else {
		oed_for_slope_over_intercept_hetero(n, xmin, xmax, theta0, f_hetero, ...)
	}
}

#private
oed_for_slope_over_intercept_homo = function(n, xmin, xmax, theta, f_hetero = NULL, ...){
	rho_to_design(n, xmin, xmax, compute_rho_star_homo(xmin, xmax, theta))
}

#private
rho_to_design = function(n, xmin, xmax, rho){
	num_xmin = round(n * rho)
	if (num_xmin == n){ #just in case...
		num_xmin = num_xmin - 1
	} else if (num_xmin == 0){
		num_xmin = 1
	}
	c(rep(xmin, num_xmin), rep(xmax, n - num_xmin))	
}

#private
plot_rho_star_by_theta = function(xmin, xmax, thetas = seq(0, 10, length.out = 100), ...){
	rhostars = array(NA, length(thetas))		
	for (i in 1 : length(thetas)){
		rhostars[i] = compute_rho_star_homo(xmin, xmax, theta = thetas[i])
	}		
	plot(thetas, rhostars, type = "l", ...)
}

#private
compute_rho_star_homo = function(xmin = 1, xmax = 10, theta = 1){
	(1 + theta * xmax) / (2 + theta * (xmax + xmin))
}

##########heteroskedastic stuff


#private
oed_for_slope_over_intercept_hetero = function(n, xmin, xmax, theta, f_hetero, MaxIter = 6000, MaxFunEvals = 6000, TolFun = 1e-6, NUM_RAND_STARTS = 50){
	#define the objection function in local scope so I can use theta and f_hetero without fear of conflict elsewhere
	Q_prop = function(xs){
		xs = as.numeric(xs)
		#create design matrix
		X = cbind(rep(1, length(xs)), xs)
		#do some intermediate calculations
		XtXinv = solve(t(X) %*% X)
		XtSigmaX = t(X) %*% diag(f_hetero(xs)) %*% X
		var_B = XtXinv %*% XtSigmaX %*% XtXinv
		del_g_beta = t(as.matrix(c(-theta, 1)))
		#asymptotic variance
		as.numeric(del_g_beta %*% var_B %*% t(del_g_beta))
	}
#	f_hetero = function(x){2 - (x - xmin) / (xmax - xmin)}
#	f_hetero = function(x){1}
	options = optimset(MaxIter = MaxIter, MaxFunEvals = MaxFunEvals, TolFun = TolFun)
	#use the homoskedastic as a starting point
	
	x_starts = rbind(
			seq(xmin, xmax, length.out = n)
	)
	for (rho in seq(0, 1, length.out = 2 * n)){
		x_starts = rbind(x_starts, rho_to_design(n, xmin, xmax, rho))
	}
	for (i in 1 : NUM_RAND_STARTS){
		x_starts = rbind(x_starts, runif(n, xmin, xmax))
	}
	x_starts = unique(x_starts)
	
	sols = list()
	#run the Nelder-Mead search from different starting points
	for (i in 1 : nrow(x_starts)){
#		cat(i, "\n")
		sols[[i]] = fminbnd(Q_prop, as.numeric(x_starts[i, ]), rep(xmin - TolFun, n), rep(xmax + TolFun, n), options, verbose = F)
	}
	
	sols_vals = lapply(sols, function(sol){neldermead.get(this = sol, key = "fopt")})
	sol = sols[[which.min(sols_vals)]]
	
	sols_vecs = round(matrix(unlist(lapply(sols, function(sol){sort(neldermead.get(this = sol, key = "xopt")[,1])})), nrow = length(sols), byrow = TRUE), 2)
	sols_vecs
	#return the result as a sorted array for convenience
	sort(neldermead.get(this = sol, key = "xopt")[, 1])	
}


#' A visualiation for comparing slope-divided-by-intercept estimates
#' for a number of designs
#' 
#' @param xmin 						The minimum value of the independent variable.
#' @param xmax 						The maximum value of the independent variable.
#' @param designs 					A d x n matrix where each of the d rows is a design (the x values
#' 									used to run the experiment).
#' @param gen_resp 					A model for the response which takes the design as its parameter.
#' @param Nsim 						The number of estimates per design. Default is \code{1000}.
#' @param l_quantile_display 		The lowest quantile of the simulation estimates displayed. Default is \code{0.025}.
#' @param u_quantile_display 		The highest quantile of the simulation estimates displayed. Default is \code{0.975}.
#' @param error_est 				The error metric for the estimates. The sample standard deviation (i.e. \code{sd}) 
#' 									is unstable at low sample sizes. The default is the 90%ile minus the 10%ile.
#' @param draw_theta_at 			If the user wishes to draw a horizontal line marking theta (to checked biasedness)
#' 									it is specified here. The default is \code{NULL} with no line being drawn.
#' 
#' @return 							A list with the simulated estimates and error estimates for each design.
#' 
#' @author 							Adam Kapelner
#' @export
design_bakeoff = function(xmin, xmax, designs, 
		gen_resp = function(xs){1 + 2 * xs + rnorm(length(xs), 0, 1)}, 
		Nsim = 1000, l_quantile_display = 0.025, u_quantile_display = 0.975,
		error_est = function(est){quantile(est, 0.9) - quantile(est, 0.1)}, #90%ile - 10%ile
		draw_theta_at = NULL){
	num_designs = nrow(designs)
	n = ncol(designs)
	
	ests = matrix(NA, nrow = Nsim, ncol = num_designs)
	
	for (nsim in 1 : Nsim){
		for (d in 1 : num_designs){
			xs = designs[d, ]
			ys = gen_resp(xs)
			mod = lm(ys ~ xs)
			ests[nsim, d] = coef(mod)[2] / coef(mod)[1]
		}
	}
	error_ests = apply(ests, 2, error_est)
	
	l = list()
	for (d in 1 : num_designs){
		l[[d]] = ests[, d]
	}
	boxplot(l, ylim = quantile(c(ests), c(l_quantile_display, u_quantile_display)), 
			main = "Bakeoff: Error Estimates for Many Designs",  
			ylab = "theta-hat", xlab = paste("Error Results:", paste(round(error_ests, 2), collapse = ", ")), names = 1 : num_designs)
	if (!is.null(draw_theta_at)){
		abline(h = draw_theta_at, col = "blue")
	}
	invisible(list(ests = l, error_ests = error_ests))	
}

#' Report the results of the experiment as well as confidence intervals.
#' 
#' @param ys 		The measurements of the response
#' @param xs 		The design
#' @param alpha		\code{1 - alpha} is the confidence of the computed intervals. Default is \code{0.05}.
#' 
#' @return			A list of the estimate as well as confidence intervals 
#' 
#' @author			Adam Kapelner
#' @export
experimental_results = function(ys, xs, alpha = 0.05){
	mod = lm(ys ~ xs)
	est = coef(mod_opt)[2] / coef(mod_opt)[1]

	list(
		est = est,
		
	)
}