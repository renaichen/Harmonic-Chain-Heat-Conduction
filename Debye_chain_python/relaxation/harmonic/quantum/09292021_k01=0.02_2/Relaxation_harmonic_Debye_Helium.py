import numpy as np
import os
import sys
import random
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import time


def single_trajectory(n):
    
    x0 = np.zeros(t_size)
    x1 = np.zeros(t_size)
    x2 = np.zeros(t_size)
    v0 = np.zeros(t_size)
    v1 = np.zeros(t_size)
    v2 = np.zeros(t_size)
    K0 = np.zeros(t_size)
    U0 = np.zeros(t_size)

    tstep = 0
    # v0[0] = random.gauss(0,1)
    # v1[0] = random.gauss(0,1)
    # v2[0] = random.gauss(0,1)
    v0[0] = 2

    x0old = x0new = 0
    x1old = x1new = 0
    x2old = x2new = 0
    v0old = v0[0]
    v0new = 0
    v1old = v1[0]
    v1new = 0
    v2old = v2[0]
    v2new = 0
    a0old = a1old = a2old = 0
    a0new = a1new = a2new = 0

    while tstep < t_size - 1:
        Rn = 2 * np.sqrt(gamma * Teff_list *deltaomega / np.pi) * np.cos(omegan * dt *tstep + phin)
        K0[tstep] += 0.5 * v0old ** 2
        U0[tstep] += 0.5 * k0 * x0old ** 2

        F01old = - k01 * (x0old - x1old)
        a0old = - k0 * x0old + F01old
        a1old = - k1 * (x1old - x2old) - F01old
        a2old = - k1 * (x2old - x1old) - k2 * x2old
        x0new = x0old + v0old * dt + 0.5 * a0old * dt ** 2
        x1new = x1old + v1old * dt + 0.5 * a1old * dt ** 2
        x2new = x2old + v2old * dt + 0.5 * a2old * dt ** 2
        F01new = - k01 * (x0new - x1new)
        a0new = -k0 * x0new + F01new
        a1new = - k1 * (x1new - x2new) - F01new
        a2new = - k1 * (x2new - x1new) - k2 * x2new
        v0new = v0old + 0.5 * (a0old + a0new) * dt
        v1new = v1old + 0.5 * (a1old + a1new) * dt
        # Teff effect
        v2new = v2old + 0.5 * (a2old + a2new) * dt - gamma * v2old * dt + np.sum(Rn) * dt
        x0[tstep + 1] = x0new
        x1[tstep + 1] = x1new
        x2[tstep + 1] = x2new
        v0[tstep + 1] = v0new
        v1[tstep + 1] = v1new
        v2[tstep + 1] = v2new

        x0old = x0new
        x1old = x1new
        x2old = x2new
        v0old = v0new
        v1old = v1new
        v2old = v2new

        # if tstep % 10000 == 0:
        #     print(tstep)
        tstep += 1

    return x0, v0, K0, U0 



if __name__ == "__main__":

    start = time.time()

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    # k1 = 0.1
    # k2 = 0.1
    # gamma = 1
    k0 = 0.1
    k1 = 0.151388
    k2 = 0.131483
    gamma = 0.868517
    # gamma = 0.
    omegaD = (k1**2*k2**2)**(1/8)
    kB = 1
    # T = 0.1
    T = float(sys.argv[1])
    k01 = 0.02

    t_start = 1
    t_end = 1e2
    dt = 1e-3
    # dt = 1 / (np.pi * 100 * omegaD)
    t_array = np.arange(t_start, t_end, dt)

    t_size = len(t_array)

    omegaN = 6 * omegaD
    Nslice = int(5e4)
    deltaomega = omegaN / Nslice
    omegan = np.linspace(1e-5, omegaN, Nslice)
    phin = 2 * np.pi * np.random.randn(Nslice)
    Teff_list = omegan / (np.exp(omegan / T) -1)

    traj = 2400

    Kin0traj = np.zeros(t_size)
    Pot0traj = np.zeros(t_size)

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.
    i = 0
    for x_t, v_t, Kin, Pot in p.map(single_trajectory, range(traj)):
        Kin0traj += Kin
        Pot0traj += Pot

    p.close()
    p.join()
    

    Kin0traj /= traj
    Pot0traj /= traj

    np.savetxt("energy_relax_T" + str(T) + ".txt", np.c_[t_array[::100], Kin0traj[::100], Pot0traj[::100]])

    print("time used: ", (time.time() - start) / 60)

    # plt.figure(1)
    # plt.plot(x0)
    # plt.plot(x1)
    # plt.plot(x2)
    # plt.plot(K0)
    # plt.figure(2)
    # plt.plot(v0)
    # plt.plot(v1)
    # plt.plot(v2)
    plt.plot(Kin0traj+Pot0traj)
    plt.savefig("energy_T" + str(T) + ".pdf")




