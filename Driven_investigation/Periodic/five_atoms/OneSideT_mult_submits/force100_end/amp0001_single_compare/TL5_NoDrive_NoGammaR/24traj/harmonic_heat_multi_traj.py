import os
import sys
import numpy as np
import random
from multiprocessing import Pool
import time


def accel_determine_update(x, m, K, dt):
    x_forward = np.roll(x, 1)  # push the position list a step forward, now it can be treated as x_{n-1}
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    acc = -(K[:-1] / m) * (x - x_forward) - (K[1:] / m) * (x - x_backward)
    acc[0] = -(K[1] / m[0]) * (x[0] - x[1]) - (K[0] / m[0]) * x[0]
    acc[-1] = -(K[-2] / m[-1]) * (x[-1] - x[-2]) - (K[-1] / m[-1]) * x[-1]
    return acc


def accel_damp_update(v):
    """return the damping part of acceleration to the left and right most atoms"""
    adamp = np.zeros(len(v))
    adamp[0] = - gammal * v[0]
    adamp[-1] = - gammar * v[-1]
    return adamp


def accel_fluct_update(m, dt, xil, xir):
    fluctuation = np.zeros(len(m))
    fluctuation[0] = np.sqrt(2 * kb * templ * gammal / (m[0] * dt)) * xil
    fluctuation[-1] = np.sqrt(2 * kb * tempr * gammar / (m[-1] * dt)) * xir

    return fluctuation


def apply_periodic(x, n, A, omega, t):
    """need the inputs for the amplitude, frequency and t for the oscillation"""
    x[n] = A * np.cos(omega * t)
    return x


def vv_update(xprev, vprev, acc_determine_prev, m, K, dt, xil, xir, tstep):
    x = xprev + vprev * dt + 0.5 * acc_determine_prev * dt ** 2
    # x[4] = amp * np.cos(omega_drive * tstep * dt)
    a_determine_new = accel_determine_update(x, m, K, dt)
    a_damp = accel_damp_update(vprev)
    a_fluct = accel_fluct_update(m, dt, xil, xir)
    v = vprev + 0.5 * (acc_determine_prev + a_determine_new) * dt \
        + a_damp * dt + a_fluct * dt
    return x, v, a_determine_new


def heat_update(vprev, vnow, m, dt, xil, xir):
    aver_vleft = 0.5 * (vprev[0] + vnow[0])
    aver_vright = 0.5 * (vprev[-1] + vnow[-1])
    ql = -gammal * m[0] *vprev[0] * aver_vleft  \
        + np.sqrt(2 * kb * templ * gammal * m[0] / dt) * xil * aver_vleft
    qr = -gammar * m[-1] * vprev[-1] * aver_vright  \
        + np.sqrt(2 * kb * tempr * gammar * m[-1] / dt) * xir * aver_vright
    return ql, qr


def heat_inter_update(x, K, v):
    """return an array of 1D interatomic heat currents"""
    x_forward = np.roll(x, 1)  # push the position list a step forward, now it can be treated as x_{n-1}
    fij = - K[:-1] * (x - x_forward)
    # the first element K_10 is not a part of the formula to be left out
    fij[0] = 0
    v_forward = np.roll(v, 1)
    return 0.5 * fij * (v + v_forward)  # again the first element needs not taken into account


def kinetic_update(m, v):
    k = m * v ** 2
    return np.sum(k)


def temperature_update(m, v):
    return m * v ** 2 / kb

def potential_update(m, x, omega):
    x_backward = np.roll(x, -1)  # push the position list a step backward, now it can be treated as x_{n+1}
    u = 0.5 *  omega[1:] * (x_backward - x) ** 2
    u_sum = np.sum(u) + 0.5 *  omega[0] * x[0] ** 2
    return u_sum


def single_trajectory(n):
    x_t = np.zeros((1,N))
    v_t = np.zeros((1,N))
    t_dis = np.array([])

    x_n_old = np.zeros(N)
    v_n_old = np.zeros(N)
    a_n_old = np.zeros(N)
    x_n_new = np.zeros(N)
    v_n_new = np.zeros(N)
    # v_n_new = np.sqrt(kb*Teff/mass)*np.random.normal(0, 1, N)
    a_n_new = np.zeros(N)
    powerL = 0.0
    powerR = 0.0
    power12 = 0.0
    power23 = 0.0
    powerNN_1 = 0.0
    Kin = 0.0
    Pot = 0.0
    temperature = np.zeros(N)

    tstep = 0
    while tstep < (t_size - 1):

        if tstep % sample_interval == 0 and tstep > half_tsize:
            x2d = x_n_new.reshape(1,N)
            x_t = np.append(x_t, x2d, axis=0)
            v2d = v_n_new.reshape(1,N)
            v_t = np.append(v_t, v2d, axis=0)
            t_dis = np.append(t_dis, deltat * tstep)


        xil = random.gauss(0, 1)
        xir = random.gauss(0, 1)

        x_n_old = x_n_new
        v_n_old = v_n_new
        a_n_old = a_n_new
        # x_n_new, v_n_new, a_n_new = (vv_update(x_n_old,v_n_old, a_n_old, mass,
        #                                         force_constants, deltat,xil, xir))
        x_n_new, v_n_new, a_n_new = (vv_update(x_n_old,v_n_old, a_n_old, mass,
                                                force_constants, deltat,xil,
                                                xir, tstep))
        # x_n_new = apply_periodic(x_n_new, 2, amp, omega_drive, tstep*deltat)
        # #No need for this apply_periodic function, since hard-coded into the
        # vv_update function

        if tstep > half_tsize:
            powerLtemp, powerRtemp = heat_update(v_n_old, v_n_new, mass, deltat, xil, xir)
            powerL += powerLtemp
            powerR += powerRtemp
            power12 += heat_inter_update(x_n_new, force_constants, v_n_new)[1]
            power23 += heat_inter_update(x_n_new, force_constants, v_n_new)[2]
            powerNN_1 += heat_inter_update(x_n_new, force_constants, v_n_new)[-1]
            Kin += kinetic_update(mass, v_n_old)
            Pot += potential_update(mass, x_n_old, force_constants)
            temperature += temperature_update(mass, v_n_old)
        tstep += 1

    return Kin, Pot, powerL, power12, power23,  powerNN_1, powerR, temperature, t_dis, x_t, v_t



if __name__ == '__main__':

    traj = 24
    SLOTS = int(os.getenv('NSLOTS')) # Get the NSLOTS environment variable provided by the scheduler
    start_time = time.time()

    # N = 2
    # mass = np.array([1.0, 1.0])
    # force_constants = np.array([0.0, 1., 0.0])
    N = 5
    mass = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    # force_constants = np.array([0.0, 1., 1., 1., 1., 0.0])
    force_constants = np.array([0.0, 100., 100., 100., 100., 0.0])

    amp = 0.0
    omega_drive = 1.
    # omega_drive = float(sys.argv[1])

    kb = 1.0
    # templ = 0.001
    # tempr = 0.005
    templ = 5.0
    tempr = 0.0
    # tempr = float(sys.argv[1])
    gammal = 1.0
    gammar = 0.0
    # Teff = 1.
    # Teff = (templ + tempr) / 2.
    Teff = (templ * gammal + tempr * gammar) / (gammal + gammar)
    tEnd = 1e6
    deltat = 0.001
    t_array = np.arange(0, tEnd, deltat)
    sample_interval = 1000
    t_size = t_array.size
    half_tsize = int(t_size // 2)

    Kintraj = np.zeros(traj)
    Pottraj = np.zeros(traj)
    powerLtraj = np.zeros(traj)
    powerRtraj = np.zeros(traj)
    power12traj = np.zeros(traj)
    power23traj = np.zeros(traj)
    powerNN_1traj = np.zeros(traj)
    temperaturetraj = np.zeros((N, traj))

    # # ##--only when necessary to investigate the trajectories
    # dirc = './atom_number' + str(N) + '/'
    # os.makedirs(os.path.dirname(dirc))

    p = Pool(processes=SLOTS)# pass the number of core to the Pool so that I know how many cores I can use.
    i = 0
    for Kin, Pot, powerL, power12, power23, powerNN_1, powerR, temperature, t_dis, x_t, v_t in p.map(single_trajectory, range(traj)):
        Kintraj[i] = Kin / half_tsize
        Pottraj[i] = Pot / half_tsize
        powerLtraj[i] = powerL / half_tsize
        power12traj[i] = power12 / half_tsize
        power23traj[i] = power23 / half_tsize
        powerNN_1traj[i] = powerNN_1 / half_tsize
        powerRtraj[i] = powerR / half_tsize
        temperaturetraj[:, i] = temperature / half_tsize

        time_dis = t_dis
        pos = x_t
        vel = v_t

        # ##--only when necessary to investigate the trajectories
        # filename = str(tempr) + str(N)  + 'position-' + str(i) + '.txt'
        # np.savetxt(filename, x_t)
        # filename = str(tempr) + str(N)  + 'velocity-' + str(i) + '.txt'
        # np.savetxt(filename, v_t)

        i += 1

    ##--only when necessary to investigate the trajectories
    # filename = str(tempr) + str(N)  + 'time_points' + str(i) + '.txt'
    # np.savetxt(filename, np.c_[time_dis])

    p.close()
    p.join()


    # ##----------
    Taver = np.mean(temperaturetraj, axis=1)
    Tstd = np.std(temperaturetraj, axis=1)
    print('T = ', Taver, Tstd)

    JLaver = np.mean(powerLtraj)
    JLstd = np.std(powerLtraj)
    print('heatL = ', JLaver, JLstd)
    J12aver = np.mean(power12traj)
    J12std = np.std(power12traj)
    print('heat12 = ', J12aver, J12std)
    J23aver = np.mean(power23traj)
    J23std = np.std(power23traj)
    print('heat23 = ', J23aver, J23std)
    JNN_1aver = np.mean(powerNN_1traj)
    JNN_1std = np.std(powerNN_1traj)
    print('heatNN_1 = ', JNN_1aver, JNN_1std)
    JRaver = np.mean(powerRtraj)
    JRstd = np.std(powerRtraj)
    print('heatR = ', JRaver, JRstd)
    #----------------

    with open ("TL-" +str(templ) + time.strftime('-%m%d-%H%M%S.txt'),"w") as outputfile:
        outputfile.write("This run uses mins {} \n".format((time.time() - \
            start_time)/60))
        outputfile.write("amplitude = {} \n".format(amp))
        outputfile.write("omega drive = {} \n".format(omega_drive))
        outputfile.write("gammal = {} \n".format(gammal))
        outputfile.write("gammar = {} \n".format(gammar))
        outputfile.write("time step = {} \n".format(deltat))
        outputfile.write("T = {}, Tstd = {} \n".format(Taver, Tstd))
        outputfile.write("heatL = {}, heatLstd = {} \n".format(JLaver, JLstd))
        outputfile.write("heat12 = {}, heat12std = {} \n".format(J12aver, J12std))
        outputfile.write("heat23 = {}, heat23std = {} \n".format(J23aver, J23std))
        outputfile.write("heatNN_1 = {}, heatNN_1std = {} \n".format(JNN_1aver, JNN_1std))
        outputfile.write("heatR = {}, heatRstd = {} \n".format(JRaver, JRstd))
    # ##-------------------------------




