library(optDesignSlopeInt)

xmin = 5 / 15
xmax = 19 / 1
n = 10
theta0 = 0.5

#sec 3.1
opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0)

#Fig1
Nsim_figs = 10000
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
		gen_resp = gen_resp, Nsim = Nsim_figs)
#to get the rmse (similar to standard error)
lapply(res$ests, function(ests){sqrt(sum((ests - theta0)^2) / Nsim_figs)})



#Fig 2
fig2A_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 1, theta0_min = 0.1, theta0_max = 3, 
				theta0 = theta0, Nsim = Nsim_figs)

fig2B_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 0.05, theta0_min = 0.001, theta0_max = 1, 
				theta0 = theta0, Nsim = Nsim_figs)
		
##opt design hetero
f_hetero = function(x){2 * (xmax - x) / (xmax - xmin)}
opt_hetero_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0, f_hetero = f_hetero)

n = 100
opt_homo_design_mock = oed_for_slope_over_intercept(n, xmin, xmax, theta0)
table(opt_homo_design_mock)
opt_hetero_design_mock = oed_for_slope_over_intercept(n, xmin, xmax, theta0, f_hetero = f_hetero)
table(opt_hetero_design_mock)

#Figure 3
designs = rbind(
		seq(xmin, xmax, length.out = n),
		c(rep(xmin, n / 2), rep(xmax, n / 2)),
		opt_homo_design_mock,
		opt_hetero_design_mock
)
gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, 1 * f_hetero(xs))
}
res = design_bakeoff(xmin, xmax, designs, draw_theta_at = theta0,
		gen_resp = gen_resp, 
		Nsim = Nsim_figs)
#to get the rmse (similar to standard error)
lapply(res$ests, function(ests){sqrt(sum((ests - theta0)^2) / Nsim_figs)})

### CI STUFF OPTIMAL HOMO
experimental_results(opt_hetero_design, gen_resp(opt_hetero_design))

#evaluate coverage probs
Nsim = 1000

cis_cover = matrix(NA, nrow = 7, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(opt_hetero_design)
	cis = experimental_results(opt_hetero_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim

### CI STUFF NAIVE HOMO
naive_design = seq(from = xmin, to = xmax, length.out = n)
experimental_results(naive_design, gen_resp(naive_design))

cis_cover = matrix(NA, nrow = 7, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(naive_design)
	cis = experimental_results(naive_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim

### CI STUFF OPTIMAL HETERO
experimental_results(opt_homo_design, gen_resp(opt_homo_design))

#evaluate coverage probs
Nsim = 1000

cis_cover = matrix(NA, nrow = 7, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(opt_homo_design)
	cis = experimental_results(opt_homo_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim

### CI STUFF NAIVE HETERO

f_hetero = function(x){.2 * (xmin - x)^2 + 1}

gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, sqrt(1 * f_hetero(xs)))
}
xs = seq(xmin, xmax, length.out = 2000)
ys = gen_resp(xs)
plot(xs, ys, xlim = c(0, xmax), ylim = c(0, max(ys)), xlab = "Vol Liq / Vol Gas", ylab = " 1 / GC")

naive_design = seq(from = xmin, to = xmax, length.out = n)
experimental_results(naive_design, gen_resp(naive_design))

cis_cover = matrix(NA, nrow = 7, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(naive_design)
	cis = experimental_results(naive_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim

naive_design = seq(from = xmin, to = xmax, length.out = n)
experimental_results(naive_design, gen_resp(naive_design))

cis_cover = matrix(NA, nrow = 7, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(naive_design)
	cis = experimental_results(naive_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim



