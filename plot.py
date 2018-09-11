import matplotlib.pyplot as plt

plt.figure(1)
plt.title("Product of two 32 by 32 matrices with 2^16 bits integers")
plt.ylabel('Time (s)')
plt.xlabel('Level 1 moduli bitsize upper bound')
input_bitsize = [11,12,13,14,15]
time = [15.2, 9.2, 8.5, 10.1, 16.3]
plt.plot(input_bitsize, time)
plt.show()
