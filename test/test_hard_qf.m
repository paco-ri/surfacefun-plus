clear

[morton, C, L_max] = example_quadforest_2();
roots = 1:length(morton);
qf = quadforest(morton, L_max, C, roots);
qf.plot_quadforest()