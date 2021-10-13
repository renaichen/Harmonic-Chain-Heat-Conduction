import numpy as np
import random
import matplotlib.pyplot as plt
import time



if __name__ == "__main__":

    start = time.time()

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
    T = 0.1
    # T = 0.
    k01 = 0.02

    t_start = 1
    t_end = 1e2
    dt = 1e-3
    # dt = 1 / (np.pi * 100 * omegaD)
    t_array = np.arange(t_start, t_end, dt)

    t_size = len(t_array)

    N = 2
    traj = 0
    sample_interval = 10
    # V1 = np.zeros((N, t_size // (2 * sample_interval)))

    K0 = np.zeros(t_size)
    U0 = np.zeros(t_size)
    
    while traj < N:

        x0 = np.zeros(t_size)
        x1 = np.zeros(t_size)
        x2 = np.zeros(t_size)
        v0 = np.zeros(t_size)
        v1 = np.zeros(t_size)
        v2 = np.zeros(t_size)

        tstep = 0
        # v0[0] = random.gauss(0,1)
        # v1[0] = random.gauss(0,1)
        # v2[0] = random.gauss(0,1)
        v0[0] = 1

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
            xi = random.gauss(0, 1)
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
            v2new = v2old + 0.5 * (a2old + a2new) * dt - gamma * v2old * dt + np.sqrt(2 * gamma * kB * T * dt) * xi
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

            if tstep % 10000 == 0:
                print(tstep)

            tstep += 1

        # V1[traj,:] = v1[t_size // 2::sample_interval]
        # V1[traj,:] = x1[t_size // 2::sample_interval]

        # print(K0[0])
        print("traj = ", traj)
        traj += 1
    K0 /= N
    U0 /= N





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
    plt.plot(K0+U0)

    np.savetxt("energy_relax.txt", np.c_[t_array[::100], K0[::100], U0[::100]])

    plt.savefig("energy.pdf")




