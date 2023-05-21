from scipy.stats import poisson
for i in range(1000000):
    p_value = poisson.cdf(100, 0.04)
print(p_value)