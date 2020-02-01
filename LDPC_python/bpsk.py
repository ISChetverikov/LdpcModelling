import numpy as np


class Bpsk:
    def modulate(self, codeword):

        result = codeword[codeword == 0] = -2
        return result
