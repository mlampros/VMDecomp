
#..................
# Parameter setting  (1-dimensional)
#..................

N = 750           # different compared to the examples section

set.seed(1)
rand_unif = runif(n = N, min = 0, max = 1.0)

f_sig1 = 6 * rand_unif
f_sig2 = cos(x = 8 * pi * rand_unif)
f_sig3 = 0.5 * cos(x = 40 * pi * rand_unif)

f_sig = f_sig1 + f_sig2 + f_sig3


#..................
# Parameter setting  (2-dimensional)
#..................

rows_cols = 20    # different compared to the examples section

set.seed(2)
data = matrix(runif(rows_cols^2), rows_cols, rows_cols)

