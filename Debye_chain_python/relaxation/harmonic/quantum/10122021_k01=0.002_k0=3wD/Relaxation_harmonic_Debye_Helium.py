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
    
    tstep = 0
    # v0[0] = random.gauss(0,1)
    # v1[0] = random.gauss(0,1)
    # v2[0] = random.gauss(0,1)
    ph = random.gauss(0,1) * 2 * np.pi
    v0old = 2 * np.cos(ph)
    x0old = 2 * np.sin(ph) / np.sqrt(k0)

    x0new = x0old
    x1old = x1new = 0
    x2old = x2new = 0
    v0new = 0
    v1old = 0
    v1new = 0
    v2old = 0
    v2new = 0
    a0old = a1old = a2old = 0
    a0new = a1new = a2new = 0
    K0 = np.array([[0, 0.5 * v0old ** 2]])
    U0 = np.array([[0, 0.5 * k0 * x0old ** 2]])

    phin = np.pi * np.random.randn(Nslice)

    while tstep < t_size:
        Rn = 2 * np.sqrt(gamma * Teff_list *deltaomega / np.pi) * np.cos(omegan * dt *tstep + phin)

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

        if tstep > 0 and tstep % 1000 == 0:
            K0 = np.append(K0, [[tstep * dt, 0.5 * v0new ** 2]], axis=0)
            U0 = np.append(U0, [[tstep * dt, 0.5 * k0 * x0new ** 2]], axis=0)

        x0old = x0new
        x1old = x1new
        x2old = x2new
        v0old = v0new
        v1old = v1new
        v2old = v2new

        tstep += 1

    # return x0, v0, K0, U0 
    return K0, U0 



if __name__ == "__main__":

    start = time.time()

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    # k1 = 0.1
    # k2 = 0.1
    # gamma = 1
    # k0 = 0.1
    k0 = 1.2697
    k1 = 0.151388
    k2 = 0.131483
    gamma = 0.868517
    # gamma = 0.
    omegaD = (k1**2*k2**2)**(1/8)
    kB = 1
    # T = 0.1
    T = float(sys.argv[1])
    k01 = 0.002

    t_start = 0
    t_end = 5e4
    dt = 1e-2
    # dt = 1 / (np.pi * 100 * omegaD)
    t_array = np.arange(t_start, t_end, dt)

    t_size = len(t_array)

    omegaN = 100 * omegaD
    Nslice = int(5e4)
    deltaomega = omegaN / Nslice
    omegan = np.linspace(1e-5, omegaN, Nslice)
    Teff_list = omegan / (np.exp(omegan / T) -1)

    traj = 96

    # Kin0traj = np.zeros(t_size)
    # Pot0traj = np.zeros(t_size)
    Kin0traj = np.zeros((t_size // 1000, 2))
    Pot0traj = np.zeros((t_size // 1000, 2))

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.
    i = 0
    for Kin, Pot in p.map(single_trajectory, range(traj)):
        Kin0traj += Kin
        Pot0traj += Pot

    p.close()
    p.join()
    

    Kin0traj /= traj
    Pot0traj /= traj

    np.savetxt("energy_relax_T" + str(T) + ".txt", np.c_[Kin0traj[:,0],Kin0traj[:,1], Pot0traj[:,1]])

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
    plt.plot(Kin0traj[:,1]+Pot0traj[:,1])
    plt.savefig("energy_T" + str(T) + ".pdf")




