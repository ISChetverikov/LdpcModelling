
class AwgnLlrAdapter:

    def __init__(self):
        self._sigma = 1

    def getsigma(self):  # метод для получения значения
        return self._sigma

    def setsigma(self, value):  # метод для присваивания нового значения
        self._sigma = value

    def delsigma(self):  # метод для удаления атрибута
        del self._sigma

    sigma = property(getsigma, setsigma, delsigma, "sigma property")

    def adapt(self, vector):
        return -2 * vector / (self.sigma ** 2)
