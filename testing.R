pacman::p_load(optDesignSlopeInt)

xmin = 5 / 15
xmax = 19 / 1
n = 10
theta0 = 0.053

opt_hetero_log_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0, f_hetero = function(x){log(x)})
#throws error

opt_homo_design = oed_for_slope_over_intercept(n, xmin, xmax, theta0)
table(opt_homo_design)


plot_info = err_vs_theta0_plot_for_homo_design(n, xmin, xmax, theta0, theta0_min = 0.001, theta0_max = 1)
plot_info

designs = rbind(
  c(rep(xmin, n / 2), rep(xmax, n / 2)),       #design A
  seq(from = xmin, to = xmax, length.out = n)  #design B
)
design_bakeoff_info = design_bakeoff(xmin, xmax, designs) #design A should win
design_bakeoff_info

xs = runif(n, xmin, xmax)
ys = 2 + 3 * xs + rnorm(n)
experimental_results_info = experimental_results(xs, ys)
experimental_results_info
