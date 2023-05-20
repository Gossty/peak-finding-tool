from scipy.stats import poisson
for i in range(10000):
    p_value = poisson.pmf(25, 0.14)
print(p_value)