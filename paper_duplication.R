library(optDesignSlopeInt)

xmin = 5 / 15
xmax = 19 / 1
n = 10
theta0 = 0.5

#sec 3.1
des = oed_for_slope_over_intercept(n, xmin, xmax, theta0)

#Fig1
designs = rbind(
		seq(xmin, xmax, length.out = n),
		runif(n, xmin, xmax),
		c(rep(xmin, n / 2), rep(xmax, n / 2)),
		oed_for_slope_over_intercept(n, xmin, xmax, theta0)
)
gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, 1)
}
res = design_bakeoff(xmin, xmax, designs, draw_theta_at = theta0,
		gen_resp = gen_resp, Nsim = 10000)
#to get the rmse (similar to standard error)
lapply(res$ests, function(ests){sqrt(sum((ests - theta0)^2) / 10000)})



#Fig 2
fig2A_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 1, theta0_min = 0.1, theta0_max = 3, 
				theta0 = theta0, Nsim = 10000)

fig2B_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 0.05, theta0_min = 0.001, theta0_max = 1, 
				theta0 = theta0, Nsim = 10000)
		

		
### CI STUFF

ys = gen_resp(des)
mod = lm(ys ~ des)
b0 = coef(mod)[1]
b1 = coef(mod)[2]
thetahat = b1 / b0 
es = mod$residuals



experimental_results(des, ys)

#evaluate coverage probs
Nsim = 1000

cis_cover = matrix(NA, nrow = 8, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(des)
	cis = experimental_results(des, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
}
