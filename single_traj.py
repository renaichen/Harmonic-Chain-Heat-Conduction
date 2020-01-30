import numpy as np
import matplotlib.pyplot as plt
import random


def accel_update(x, v, m, K, dt):
    """This is to update the acceleration for every atom at specific time,
    some parameters (ie. kb templ) are not passed but just as global
    defined in the main function"""

    x_forward = np.roll(x, 1)  # push the position list a step forward, now it can be treated as x_{n-1}
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    acc = -(K[:-1] / m) * (x - x_forward) - (K[1:] / m) * (x - x_backward)
    acc[0] = -(K[1] / m[0]) * (x[0] - x[1]) - (K[0] / m[0]) * x[0] - gammal * v[0] \
        + np.sqrt(2 * kb * templ * gammal / (m[0] * dt)) * random.gauss(0, 1)
    acc[-1] = -(K[-2] / m[-1]) * (x[-1] - x[-2]) - (K[-1] / m[-1]) * x[-1] - gammar * v[-1]  \
        + np.sqrt(2 * kb * tempr * gammar / (m[-1] * dt)) * random.gauss(0, 1)

    return acc


def vv_update(xprev, vprev, accprev, m, K, dt):
    x = xprev + vprev * dt + 0.5 * accprev * dt ** 2
    a = accel_update(x, vprev, m, K, dt)
    v = vprev + 0.5 * (accprev + a) * dt
    print(v)

    return x, v, a


def heat_update(vprev, vnow, dt):
    ql = (-gammal * 0.5 * (vprev[0] + vnow[0]) * dt +
          np.sqrt(2 * kb * templ * gammal * (mass[0] * dt)) * xil * 0.5 * (vprev[0] + vnow[0]))
    qr = (-gammar * 0.5 * (vprev[-1] + vnow[-1]) * dt +
          np.sqrt(2 * kb * tempr * gammar * (mass[-1] * dt)) * xir * 0.5 * (vprev[-1] + vnow[-1]))

    return ql, qr


def kinetic_update(m, v):
    k = 0.5 * m * v ** 2

    return np.sum(k)


def potential_update(m, x, omega):
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    u = 0.5 * m * omega[1:] ** 2 * (x_backward - x) ** 2
    u_sum = np.sum(u) + 0.5 * m[0] * omega[0] ** 2 * x[0] ** 2

    return u_sum


N = 2
mass = np.array([1.0, 1.0])
force_constants = np.array([0.0, 0.5, 0.0])

kb = 1.0
templ = 0
tempr = 0
gammal = 1.0
gammar = 1.0
# Teff = 1.
# Teff = (templ + tempr) / 2.
Teff = (templ * gammal + tempr * gammar) / (gammal + gammar)

tEnd = 1
deltat = 0.1
t_array = np.arange(0, tEnd, deltat)
t_size = t_array.size
x_t = np.zeros((N, t_size))
v_t = np.zeros((N, t_size))
# v_t[:, 0] = np.sqrt(kb*Teff/mass)*np.random.normal(0, 1, N)
v_t[0, 0] = 0.1
v_t[1, 0] = 0.2
a_t = np.zeros((N, t_size))
ql_t_instant = np.zeros(t_size)
qr_t_instant = np.zeros(t_size)
ql_t_accu = np.zeros(t_size)
qr_t_accu = np.zeros(t_size)
k_t_tot = np.zeros(t_size)
u_t_tot = np.zeros(t_size)

# x_t[0, 0] = 0.1
ql_sum = 0.0
qr_sum = 0.0

for i in range(t_size - 1):
    xil = random.gauss(0, 1)
    xir = random.gauss(0, 1)
    x_t[:, i + 1], v_t[:, i + 1], a_t[:, i + 1] = (vv_update(x_t[:, i],
                                                             v_t[:, i], a_t[:, i], mass, force_constants, deltat))
    k_t_tot[i] = kinetic_update(mass, v_t[:, i])
    # u_t_tot[i] = potential_update(mass, x_t[:, i], v_t[:, i])
    u_t_tot[i] = potential_update(mass, x_t[:, i], force_constants)
    ql_t_instant[i + 1], qr_t_instant[i + 1] = heat_update(v_t[:, i], v_t[:, i + 1], deltat)
    if i > (t_size//2):
        ql_sum += ql_t_instant[i + 1]
        ql_t_accu[i + 1] = ql_sum
        qr_sum += qr_t_instant[i + 1]
        qr_t_accu[i + 1] = qr_sum

    print(i)
print(ql_sum)
print(qr_sum)
# print(np.sum(ql_t_accu[-1]))
# print(np.sum(qr_t_accu[-1]))

# plt.figure(1)
# plt.plot(t_array, x_t[0, :])
# plt.plot(t_array, x_t[1, :])
# plt.plot(t_array, x_t[0, :] + x_t[1, :])
# plt.plot(t_array, v_t[0, :])
# plt.plot(t_array, v_t[1, :])
# plt.plot(t_array, v_t[0, :] + v_t[1, :])
#
# plt.figure(2)
# plt.plot(t_array, ql_t_accu)
# plt.plot(t_array, qr_t_accu)
#
# plt.figure(3)
# plt.plot(t_array, k_t_tot + u_t_tot)
#
# plt.figure(4)
# plt.plot(t_array, u_t_tot)

# plt.show()
