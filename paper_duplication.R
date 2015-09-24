library(optDesignSlopeInt)



##### Section 3.2

xmin = 5 / 15
xmax = 19 / 1
n = 10
theta0 = 0.5
opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0)

##Fig 1
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



##Fig 2
fig2A_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 1, theta0_min = 0.1, theta0_max = 3, 
				theta0 = theta0, Nsim = Nsim_figs)

fig2B_res = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, 
				theta = 0.05, theta0_min = 0.001, theta0_max = 1, 
				theta0 = theta0, Nsim = Nsim_figs)

##### Section 3.2

#opt design hetero
f_hetero = function(x){2 * (xmax - x) / (xmax - xmin)}
opt_hetero_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0, f_hetero = f_hetero)

#note that opt_hetero_design = opt_homo_design, thus do mock case opt design hetero
n = 100
opt_homo_design_mock = oed_for_slope_over_intercept(n, xmin, xmax, theta0)
table(opt_homo_design_mock)
opt_hetero_design_mock = oed_for_slope_over_intercept(n, xmin, xmax, theta0, f_hetero = f_hetero)
table(opt_hetero_design_mock)

##Fig 3
sigsq = 1

designs = rbind(
		seq(xmin, xmax, length.out = n),
		c(rep(xmin, n / 2), rep(xmax, n / 2)),
		opt_homo_design_mock,
		opt_hetero_design_mock
)
gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, sqrt(sigsq) * f_hetero(xs))
}
res = design_bakeoff(xmin, xmax, designs, 
		draw_theta_at = theta0,
		gen_resp = gen_resp, 
		Nsim = Nsim_figs)
#to get the rmse (similar to standard error)
lapply(res$ests, function(ests){sqrt(sum((ests - theta0)^2) / Nsim_figs)})

##Table 1

gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, sqrt(sigsq))
}

#example interval
experimental_results(opt_homo_design, gen_resp(opt_homo_design))

#evaluate coverage probs
Nsim = 1000
cis_cover = matrix(NA, nrow = 6, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(opt_homo_design)
	cis = experimental_results(opt_homo_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim









### NOT IN PAPER

##Table 2

gen_resp = function(xs, beta0 = 1){
	beta0 + beta0 * theta0 * xs + rnorm(length(xs), 0, sqrt(sigsq) * f_hetero(xs))
}
#example interval
experimental_results(opt_hetero_design, gen_resp(opt_hetero_design))

#evaluate coverage probs
Nsim = 1000
cis_cover = matrix(NA, nrow = 6, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(opt_hetero_design)
	cis = experimental_results(opt_hetero_design, ys)$cis
	cis_cover[, i] = apply(cis, 1, function(ci){ci[1] <= theta0 && ci[2] >= theta0})
	cat(".")
}
rowSums(cis_cover == TRUE, na.rm = TRUE) / Nsim

##Table 2'

#evaluate coverage probs
Nsim = 1000
cis_cover = matrix(NA, nrow = 6, ncol = Nsim)
for (i in 1 : Nsim){
	ys = gen_resp(opt_hetero_design_mock)
	cis = experimental_results(opt_hetero_design_mock, ys)$cis
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



library(optDesignSlopeInt)

#the "true value" of the k_H of toluene
theta = 0.05
#load the raw data in and plot
X = as.data.frame(read.csv("~/Downloads/toluene.csv"))
plot(X$x, X$y)
#make a judgment call on what the x-range is for GC instrument
xmin = 0.3379
xmax = 14.437
#only consider the data in this range
Xr = X[X$x > xmin & X$x < xmax, ]
plot(Xr$x, Xr$y)
#run a regression to get a guess of what sigma is only
mod = lm(Xr$y ~ Xr$x)
summary(mod)
beta0 = coef(mod)[1]
sigma = summary(mod)$sigma

#in the future, we will be running experiments with 10 vials
n = 10

#investigate how our guess of theta0 will affect our results. Use a range of 0.01 to 10 since
#we are almost absolutely sure the k_H must fall in this window.
err_det = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, theta = theta, 
		Nsim = 10000, sigma = sigma,
		theta0_min = 0.01, theta0_max = 10, plot_rhos = TRUE)

#as an example show what happens if theta is 0.1 instead and graph the truth in green
err_det = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, theta = 0.1, 
		Nsim = 10000, sigma = sigma, theta0 = theta,
		theta0_min = 0.01, theta0_max = 10, plot_rhos = TRUE)


opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta)


designs = rbind(
		seq(xmin, xmax, length.out = n),
		c(rep(xmin, n / 2), rep(xmax, n / 2)),
		opt_homo_design
)


gen_resp = function(xs){
	beta0 + beta0 * theta * xs + rnorm(length(xs), 0, sigma)
}
par(mfrow = c(1, 1))
res = design_bakeoff(xmin, xmax, designs, 
		draw_theta_at = theta,
		gen_resp = gen_resp, 
		Nsim = 10000)

sd(res$ests[[3]]) / sd(res$ests[[1]])




###################

#the "true value" of the k_H of naphthalene
theta = 0.05
#load the raw data in and plot
X = as.data.frame(read.csv("C:\\Users\\Kapelner\\Desktop\\Dropbox\\oed_kh_project\\paper\\naphthalene.csv"))
plot(X$x, X$y)
#make a judgment call on what the x-range is for GC instrument
xmin = 0.3379
xmax = 14.437
#only consider the data in this range
Xr = X[X$x > xmin & X$x < xmax, ]
plot(Xr$x, Xr$y)
#run a regression to get a guess of what sigma is only
mod = lm(Xr$y ~ Xr$x)
summary(mod)
beta0 = coef(mod)[1]
sigma = summary(mod)$sigma



#in the future, we will be running experiments with 10 vials
n = 10

library(optDesignSlopeInt)

#investigate how our guess of theta0 will affect our results. Use a range of 0.01 to 10 since
#we are almost absolutely sure the k_H must fall in this window.
err_det = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, theta = theta, 
		Nsim = 10000, sigma = sigma,
		theta0_min = 0.01, theta0_max = 10, plot_rhos = TRUE)

#as an example show what happens if theta is 0.5 instead and graph the truth in green
err_det = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, theta = 0.5, 
		Nsim = 10000, sigma = sigma, theta0 = theta,
		theta0_min = 0.01, theta0_max = 10, plot_rhos = TRUE)


opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta)


designs = rbind(
		seq(xmin, xmax, length.out = n),
		c(rep(xmin, n / 2), rep(xmax, n / 2)),
		opt_homo_design
)


gen_resp = function(xs){
	beta0 + beta0 * theta * xs + rnorm(length(xs), 0, sigma)
}
par(mfrow = c(1, 1))
res = design_bakeoff(xmin, xmax, designs, 
		draw_theta_at = theta,
		gen_resp = gen_resp, 
		num_digits_round = 4,
		Nsim = 10000)

sd(res$ests[[3]]) / sd(res$ests[[1]])

