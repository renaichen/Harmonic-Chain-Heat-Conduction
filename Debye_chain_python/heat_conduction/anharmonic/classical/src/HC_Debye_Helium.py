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
    ph = random.gauss(0,1) * 2 * np.pi
    v0old = 2 * np.cos(ph)
    x0old = 2 * np.sin(ph) / np.sqrt(k0)

    x0new = x0old
    x1old = x1new = 0
    x2old = x2new = 0
    ## let "L" represent left side, while the "regular" the right side
    Lx1old = Lx1new = 0
    Lx2old = Lx2new = 0
    v0new = 0
    v1old = 0
    v1new = 0
    v2old = 0
    v2new = 0
    Lv1old = Lv1new = Lv2old = Lv2new = 0
    a0old = a1old = a2old = 0
    a0new = a1new = a2new = 0
    La1old = La1new = La2old = La2new = 0

    J01_total = LJ01_total = 0
    j01_inter = np.array([[0, 0]])
    Lj01_inter = np.array([[0, 0]])

    while tstep < t_size:
        xiL = random.gauss(0, 1)
        xiR = random.gauss(0, 1)

        LFmorseold = 2 * alphaL * DL * np.exp(-alphaL * (Lx1old - x0old)) * (np.exp(-alphaL * (Lx1old - x0old)) - 1)
        Fmorseold = 2 * alphaR * DR * np.exp(-alphaR * (x1old - x0old)) * (np.exp(-alphaR * (x1old - x0old)) - 1)
        a0old = - k0 * x0old + Fmorseold + LFmorseold
        a1old = - k1 * (x1old - x2old) - Fmorseold
        a2old = - k1 * (x2old - x1old) - k2 * x2old
        La1old = - Lk1 * (Lx1old - Lx2old) - LFmorseold
        La2old = - Lk1 * (Lx2old - Lx1old) - Lk2 * Lx2old
        x0new = x0old + v0old * dt + 0.5 * a0old * dt ** 2
        x1new = x1old + v1old * dt + 0.5 * a1old * dt ** 2
        x2new = x2old + v2old * dt + 0.5 * a2old * dt ** 2
        Lx1new = Lx1old + Lv1old * dt + 0.5 * La1old * dt ** 2
        Lx2new = Lx2old + Lv2old * dt + 0.5 * La2old * dt ** 2
        LFmorsenew = 2 * alphaL * DL * np.exp(-alphaL * (Lx1new - x0new)) * (np.exp(-alphaL * (Lx1new - x0new)) - 1)
        Fmorsenew = 2 * alphaR * DR * np.exp(-alphaR * (x1new - x0new)) * (np.exp(-alphaR * (x1new - x0new)) - 1)
        a0new = -k0 * x0new + Fmorsenew + LFmorsenew
        a1new = - k1 * (x1new - x2new) - Fmorsenew
        a2new = - k1 * (x2new - x1new) - k2 * x2new
        La1new = - Lk1 * (Lx1new - Lx2new) - LFmorsenew
        La2new = - Lk1 * (Lx2new - Lx1new) - Lk2 * Lx2new
        v0new = v0old + 0.5 * (a0old + a0new) * dt
        v1new = v1old + 0.5 * (a1old + a1new) * dt
        v2new = v2old + 0.5 * (a2old + a2new) * dt - gammaR * v2old * dt + np.sqrt(2 * gammaR * kB * TR * dt) * xiR
        Lv1new = Lv1old + 0.5 * (La1old + La1new) * dt
        Lv2new = Lv2old + 0.5 * (La2old + La2new) * dt - gammaL * Lv2old * dt + np.sqrt(2 * gammaL * kB * TL * dt) * xiL

        if tstep > (t_size // 2):
            j01inter = 0.5 * Fmorsenew * (v0new + v1new)
            Lj01inter = 0.5 * LFmorsenew * (v0new + Lv1new)
            J01_total += j01inter
            LJ01_total += Lj01inter
            if tstep % 1000 == 0:
                j01_inter = np.append(j01_inter, [[tstep * dt, j01inter]], axis=0)
                Lj01_inter = np.append(Lj01_inter, [[tstep * dt, Lj01inter]], axis=0)

        x0old = x0new
        x1old = x1new
        x2old = x2new
        Lx1old = Lx1new
        Lx2old = Lx2new
        v0old = v0new
        v1old = v1new
        v2old = v2new
        Lv1old = Lv1new
        Lv2old = Lv2new

        tstep += 1

    return J01_total * 2 / t_size, LJ01_total * 2 / t_size, j01_inter, Lj01_inter 



if __name__ == "__main__":

    start = time.time()

    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler

    k0 = 0.1
    k1 = Lk1 = 0.151388
    k2 = Lk2 = 0.131483
    gammaL = 0.868517
    gammaR = 0.868517
    # gamma = 0.
    omegaD = (k1**2*k2**2)**(1/8)
    kB = 1
    TL = 0.01
    # TR = float(sys.argv[1])
    TR = 0.05
    # k01 = 0.002
    # Lk01 = 0.002
    alphaL = 0.1 
    DL = 0.1 # 2*alpha^2*D=0.002
    alphaR = 0.1 
    DR = 0.1 # 2*alpha^2*D=0.002

    t_start = 0
    t_end = 5e3
    dt = 1e-2
    # dt = 1 / (np.pi * 100 * omegaD)
    t_array = np.arange(t_start, t_end, dt)

    t_size = len(t_array)

    traj = 48

    Kin0traj = np.zeros((t_size // 1000, 2))
    Pot0traj = np.zeros((t_size // 1000, 2))
    J01traj = LJ01traj = 0
    j01traj = np.zeros((t_size // 2000, 2)) # 2000 means half tsize
    Lj01traj = np.zeros((t_size // 2000, 2)) # 2000 means half tsize

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.
    # p = Pool(processes=2)# pass the number of core to the Pool so that I know how many cores I can use.
    i = 0
    for J01, LJ01, j01, Lj01 in p.map(single_trajectory, range(traj)):
        J01traj += J01
        LJ01traj += LJ01
        j01traj += j01
        Lj01traj += Lj01

    p.close()
    p.join()
    

    J01traj /= traj
    LJ01traj /= traj
    j01traj /= traj
    Lj01traj /= traj

    np.savetxt("hc_T" + str(TR) + ".txt", np.c_[j01traj[1:,0],j01traj[1:,1], Lj01traj[1:,1]])

    print("left and rigth currents are: ", LJ01traj, J01traj, "for T = ", TR)
    print("time used: ", (time.time() - start) / 60)

    plt.plot(j01traj[1:,0], j01traj[1:,1])
    plt.plot(Lj01traj[1:,0],Lj01traj[1:,1])
    plt.savefig("currents_T" + str(TR) + ".pdf")




