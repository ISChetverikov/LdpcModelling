import numpy as np
from exceptions import DimensionsError

class OnmsDecoder:
    def __init__(self, matrix, multiplier=0.72, offset=0, max_iteration=20):
        self.c = multiplier
        self.a = offset

        m = matrix.shape[0]
        n = matrix.shape[1]

        bits = [[] for i in range(n)]
        checks = [[] for i in range(m)]

        for i in range(m):
            for j in range(n):
                if matrix[i][j] != 1:
                    continue

                checks[i].append(j)
                bits[j].append(i)

        self.m = m
        self.n = n

        self.checks = np.array(checks)
        self.bits = np.array(bits)

        self.max_iteration = max_iteration

        return

    def __min_sum_func(self, values):
        value = self.c * min(values) - self.a
        return value if value > 0 else -value

    def __horizontal_step(self, alpha, beta, gamma):
        for j in range(self.m):

            for i in self.checks[j]:

                sign = 1.0
                values = []

                for k in self.checks[j]:
                    if k == i:
                        continue

                    sign *= alpha[j][k]
                    values.append(beta[j][k])

                gamma[j][i] = sign * self.__min_sum_func(values)

        return

    def decode(self, llr):

        n = self.n
        m = self.m

        if self.n != len(llr):
            raise DimensionsError("The codeword is not from a code with given check matrix")

        alpha0 = np.copysign(np.ones(n), llr)
        beta0 = np.abs(llr)

        alpha = [{} for _ in range(m)]
        beta = [{} for _ in range(m)]
        gamma = [{} for _ in range(m)]

        for j in range(m):
            for i in self.checks[j]:
                alpha[j][i] = np.copysign(1, llr[i])
                beta[j][i] = np.abs(llr[i])

        iterations = 0
        is_success = False

        while True:
            iterations += 1

            self.__horizontal_step(alpha, beta, gamma)

            bits_values = alpha0 * beta0
            for i in range(n):
                for j in self.bits[i]:
                    bits_values[i] += gamma[j][i]

            alpha0 = np.copysign(np.ones(n), bits_values)
            beta0 = np.abs(bits_values)

            result = np.zeros(n)
            result[alpha0 == -1.0] = 1

            is_success = True
            for check in self.checks:
                if result[check].sum() % 2 != 0:
                    is_success = False
                    break
            
            if is_success:
                break

            if iterations >= self.max_iteration:
                break

        for i in range(n):
            value = bits_values[i]

            for j in self.bits[i]:
                new_value = value - gamma[j][i]
                alpha[j][i] = np.copysign(1, new_value)
                beta[j][i] = np.abs(new_value)

        return is_success, result
