import numpy as np


class Awgn:
    def __init__(self):
        self._sigma = 1
        self.MU = 0

        return

    def getsigma(self):  # метод для получения значения
        return self._sigma

    def setsigma(self, value):  # метод для присваивания нового значения
        self._sigma = value

    def delsigma(self):  # метод для удаления атрибута
        del self._sigma

    sigma = property(getsigma, setsigma, delsigma, "sigma property")

    def simulate(self, codeword):
        e = np.random.normal(self.MU, self._sigma, len(codeword))

        return codeword + e

