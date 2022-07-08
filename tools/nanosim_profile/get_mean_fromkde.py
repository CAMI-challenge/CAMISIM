import joblib
import sys
import numpy as np
from scipy.integrate import quad
from sklearn.neighbors import KernelDensity

def integrate_mean(f):
    read_length_file = f
    #default is {prefix}_aligned_reads.pkl

    kde = joblib.load(read_length_file) # length is stored as joblib pkl

    # the kd has a probability density function from which we can get mean and variance via integration
    # it is the log density function though, need to np.exp
    pdf = lambda x : np.exp(kde.score_samples([[x]]))[0]
    mean = quad(lambda x: x * pdf(x), a=-np.inf, b=np.inf)[0]

    return mean

