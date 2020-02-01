
class AwgnLlrAdapter:

    def __init__(self, sigma):
        self.sigma = sigma
        return

    def adapt(self, vector):
        return -2 * vector / (self.sigma ** 2)

