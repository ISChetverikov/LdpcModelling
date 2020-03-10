import matplotlib.pyplot as plt


f = open('sim_res.txt', 'r')
snrs = []
fers = []
for line in f:
    values = line.split(' ')
    snrs.append(float(values[0]))
    fers.append(float(values[1]))
plt.plot(snrs, fers)
plt.title('FER dependency on SNR')
plt.xlabel('SNR')
plt.ylabel('FER')
plt.savefig('fer_curve.png')