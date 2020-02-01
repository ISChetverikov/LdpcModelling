import numpy as np


class Awgn:
    def __init__(self, sigma):
        self.sigma = sigma
        self.MU = 0

        return

    def simulate(self, codeword):
        e = np.random.normal(self.MU, self.sigma, len(codeword))

        return codeword + e

