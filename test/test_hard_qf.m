clear

[morton, C, L_max] = example_quadforest_2();
qf = quadforest(morton, L_max, C);
qf = qf.balance_quadforest();
qf.plot_quadforest()