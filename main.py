import numpy as np

n = 10  # number of periods
r_00 = 0.05  # initial short rate (5%)
u = 1.1  # upward factor
d = 0.9  # downward factor
q_u = 0.5  # upward probability
q_d = 0.5  # downward probability
F = 100  # face value of the bond
R = 20  # recovery value
a = 0.01  # hazard rate parameter
b = 1.01  # hazard rate growth factor

r_lattice = np.zeros((n + 1, n + 1))
r_lattice[0, 0] = r_00

for i in range(1, n + 1):
    for j in range(i + 1):
        r_lattice[i, j] = r_00 * (u ** j) * (d ** (i - j))

h_lattice = np.zeros((n + 1, n + 1))
for i in range(n + 1):
    for j in range(i + 1):
        h_lattice[i, j] = a * (b ** (j - i/2))

Z_lattice = np.zeros((n + 1, n + 1))

Z_lattice[n, :] = F

for i in range(n - 1, -1, -1):
    for j in range(i + 1):
        r_ij = r_lattice[i, j]
        h_ij = h_lattice[i, j]
        Z_lattice[i, j] = (
            1 / (1 + r_ij) * (
                q_u * (1 - h_ij) * Z_lattice[i + 1, j + 1] +
                q_d * (1 - h_ij) * Z_lattice[i + 1, j]
            ) +
            1 / (1 + r_ij) * (
                q_u * h_ij * R +
                q_d * h_ij * R
            )
        )

Z_00 = Z_lattice[0, 0]
print(Z_00)