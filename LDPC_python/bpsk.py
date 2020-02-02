import numpy as np


class Bpsk:
    def modulate(self, codeword):

        result = np.array(codeword)
        result[result == 0] = -1
        return result
