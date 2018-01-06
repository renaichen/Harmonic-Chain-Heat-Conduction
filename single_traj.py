
import numpy as np
import matplotlib.pyplot as plt
import random


def accel_update(x, v, m, omega, dt):
    """This is to update the acceleration for every atom at specific time,
    some parameters (ie. kb templ) are not passed but just as global
    defined in the main function"""

    x_forward = np.roll(x, 1)  # push the position list a step forward, now it can be treated as x_{n-1}
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    acc = -omega[:-1]**2*(x-x_forward) - omega[1:]**2*(x-x_backward)
    acc[0] = -omega[1]**2*(x[0]-x[1]) - omega[0]**2*x[0] - gammal*v[0] + \
             np.sqrt(2*kb*templ*gammal/(m[0]*dt))*xil
    acc[-1] = -omega[-2]**2*(x[-1]-x[-2]) - omega[-1]**2*x[-1] - gammar*v[-1] + \
               np.sqrt(2*kb*tempr*gammar/(m[-1]*dt))*xir

    return acc


def vv_update(xprev, vprev, accprev, m, omega, dt):
    x = xprev + vprev * dt + 0.5 * accprev * dt**2
    a = accel_update(x, vprev, m, omega, dt)
    v = vprev + 0.5 * (accprev + a) * dt

    return x, v, a


def heat_update():
    pass


N = 2
mass = np.array([1.0, 1.0])
freq = np.array([0.0, 1.0, 0.0])

tEnd = 10
deltat = 0.01
t_array = np.arange(0, tEnd, deltat)
t_size = t_array.size
x_t = np.zeros((N, t_size))
v_t = np.zeros((N, t_size))
a_t = np.zeros((N, t_size))

kb = 1.0
templ = 1.0
tempr = 3.0
gammar = 1.0
gammal = 1.0

x_t[0, 0] = 0.1

for i in range(t_size-1):

    xil = random.gauss(0, 1)
    xir = random.gauss(0, 1)
    x_t[:, i+1], v_t[:, i+1], a_t[:, i+1] = vv_update(x_t[:, i],
                                             v_t[:, i], a_t[:, i], mass, freq, deltat)

    print(i)


# x = [0,1,2,4]
# y = [0, 3, 5, 7]
plt.plot(t_array, v_t[1, :])
# plt.plot(x, y)
plt.show()
