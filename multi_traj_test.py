
import numpy as np
import matplotlib.pyplot as plt
import random
from multiprocessing import Pool
import time 


def accel_update(x, v, m, omega, dt, xil, xir):
    """This is to update the acceleration for every atom at specific time,
    some parameters (ie. kb templ) are not passed but just as global
    defined in the main function"""

    x_forward = np.roll(x, 1)  # push the position list a step forward, now it can be treated as x_{n-1}
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    acc = -omega[:-1]**2*(x-x_forward) - omega[1:]**2*(x-x_backward)
    acc[0] = (-omega[1]**2*(x[0]-x[1]) - omega[0]**2*x[0] - gammal*v[0] +
              np.sqrt(2*kb*templ*gammal/(m[0]*dt))*xil)
    acc[-1] = (-omega[-2]**2*(x[-1]-x[-2]) - omega[-1]**2*x[-1] - gammar*v[-1] +
               np.sqrt(2*kb*tempr*gammar/(m[-1]*dt))*xir)

    return acc


def vv_update(xprev, vprev, accprev, m, omega, dt, xil, xir):
    x = xprev + vprev * dt + 0.5 * accprev * dt**2
    a = accel_update(x, vprev, m, omega, dt, xil, xir)
    v = vprev + 0.5 * (accprev + a) * dt

    return x, v, a


def heat_update(vprev, vnow, dt, xil, xir):
    ql = (-gammal * 0.5 * (vprev[0] + vnow[0]) * dt +
          np.sqrt(2 * kb * templ * gammal * (mass[0] * dt)) * xil * 0.5 * (vprev[0] + vnow[0]))
    qr = (-gammar * 0.5 * (vprev[-1] + vnow[-1]) * dt +
          np.sqrt(2 * kb * tempr * gammar * (mass[0] * dt)) * xir * 0.5 * (vprev[-1] + vnow[-1]))

    return ql, qr


def kinetic_update(m, v):
    k = 0.5 * m * v**2

    return np.sum(k)


def potential_update(m, x, omega):
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    u = 0.5 * m * omega[1:]**2 * (x_backward - x)**2
    u_sum = np.sum(u) + 0.5 * m[0] * omega[0]**2 * x[0]**2

    return u_sum


def single_trajectory(n):

    x_t = np.zeros((N, t_size))
    v_t = np.zeros((N, t_size))
    a_t = np.zeros((N, t_size))
    ql_t_instant = np.zeros(t_size)
    qr_t_instant = np.zeros(t_size)
    ql_t_accu = np.zeros(t_size)
    qr_t_accu = np.zeros(t_size)
    k_t_tot = np.zeros(t_size)
    u_t_tot = np.zeros(t_size)

    x_t[0, 0] = 0.1
    ql_sum = 0.0
    qr_sum = 0.0

    for i in range(t_size-1):

        xil = random.gauss(0, 1)
        xir = random.gauss(0, 1)
        x_t[:, i + 1], v_t[:, i + 1], a_t[:, i + 1] = (vv_update(x_t[:, i],
                                                                 v_t[:, i], a_t[:, i], mass, freq, deltat, xil, xir))
        k_t_tot[i] = kinetic_update(mass, v_t[:, i])
        u_t_tot[i] = potential_update(mass, x_t[:, i], v_t[:, i])
        ql_t_instant[i + 1], qr_t_instant[i + 1] = heat_update(v_t[:, i], v_t[:, i + 1], deltat, xil, xir)
        ql_sum += ql_t_instant[i + 1]
        ql_t_accu[i + 1] = ql_sum
        qr_sum += qr_t_instant[i + 1]
        qr_t_accu[i + 1] = qr_sum

        print(i)

    return x_t, v_t, k_t_tot, u_t_tot, ql_t_instant, qr_t_instant, ql_t_accu, qr_t_accu


if __name__ == '__main__':

    start_time = time.time()
    
    N = 2
    mass = np.array([1.0, 1.0])
    freq = np.array([0.0, 1.0, 0.0])
    kb = 1.0
    templ = 10.0
    tempr = 30.0
    gammar = 100.0
    gammal = 100.0

    traj = 1000
    t_end = 100
    deltat = 0.01
    t_array = np.arange(0, t_end, deltat)
    t_size = t_array.size

    x_t = np.zeros((N, t_size))
    v_t = np.zeros((N, t_size))
    u_t_tot = np.zeros(t_size)
    k_t_tot = np.zeros(t_size)
    ql_t_instant = np.zeros(t_size)
    qr_t_instant = np.zeros(t_size)
    ql_t_accu = np.zeros(t_size)
    qr_t_accu = np.zeros(t_size)

    p = Pool()
    for x, v, k_tot, u_tot, ql_instant, qr_instant, ql_accu, qr_accu in p.map(single_trajectory, range(traj)):
        x_t += x
        v_t += v
        u_t_tot += u_tot
        k_t_tot += k_tot
        ql_t_instant += ql_instant
        qr_t_instant += qr_instant
        ql_t_accu += ql_accu
        qr_t_accu += qr_accu

    p.close()
    p.join()
    
    print("This run uses %s seconds " % (time.time() - start_time))

    # plt.figure(1)
    # plt.plot(t_array, v_t[0, :])
    #
    # plt.figure(2)
    # plt.plot(t_array, ql_t_accu)
    #
    # plt.figure(3)
    # plt.plot(t_array, k_t_tot)

    #plt.figure(5)
    #plt.plot(t_array, u_t_tot/traj)


    #plt.show()
