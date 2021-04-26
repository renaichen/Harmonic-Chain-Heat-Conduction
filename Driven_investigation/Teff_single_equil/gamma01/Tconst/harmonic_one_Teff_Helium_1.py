import os
import sys
import numpy as np
import random
import time


##------------define global variables here
## for whatever reason, __name__=='__main__'
## doesn't work in local PC (though on the cluster----------
start_time = time.time()
N = 1
Nslice = int(5e3)

gammal = 0
gammar = 0.1
templ = 0.0
tempr = float(sys.argv[1])
tEnd = 1e4
deltat = 1e-3
t_array = np.arange(0, tEnd, deltat)
sample_interval = 1
t_size = t_array.size
half_tsize = int(t_size // 2)
# omegaN = 2 * np.pi / deltat
omegaN = 2 * np.pi * 10
deltaomega = omegaN / Nslice
omegan = np.arange(1e-15, omegaN, deltaomega)
# omegan = np.linspace(1e-15, omegaN, Nslice)
phin = np.pi * np.random.randn(Nslice)
# Teff_list = omegan / (np.exp(omegan /  tempr) - 1)
Teff_list = tempr
##----------------------------------------------

def accel_fluct_update(dt):
    # fluctuation = np.zeros(len(m))
    fluctuation = np.sqrt(2 * templ * gammal / (dt)) * xil
    fluctuation += np.sqrt(2 * tempr * gammar / (dt)) * xir
    return fluctuation

# def Teff(omegan, T):
#     return omegan / (np.exp(omegan /  T) - 1)

# # def Rt(omegan, T, t):
#     # Teff_list = Teff(omegan, T)
#     Rn = np.sqrt(2 * gammar * Teff_list * deltaomega / np.pi) * np.cos(omegan * t + phin)
#     # print(np.cos(omegan * t + phin))
#     return np.sum(Rn)


def vv_update(xprev, vprev, acc_determine_prev, dt, tstep):
    x = xprev + vprev * dt + 0.5 * acc_determine_prev * dt ** 2
    a_determine_new = - x
    a_damp = - (gammal + gammar) * vprev
    # a_fluct = accel_fluct_update(dt)
    # a_fluct = Rt(omegan, tempr, dt * tstep)
    Rn = 2 * np.sqrt(gammar * Teff_list * deltaomega / np.pi) * np.cos(omegan * dt * tstep + phin)
    a_fluct = np.sum(Rn)
    v = vprev + 0.5 * (acc_determine_prev + a_determine_new) * dt \
        + a_damp * dt + a_fluct * dt
    return x, v, a_determine_new

def heat_update(vprev, vnow, dt):
    aver_vleft = 0.5 * (vprev + vnow)
    aver_vright = 0.5 * (vprev + vnow)
    # ql = -gammal * vprev * aver_vleft  
    # qr = -gammar * vprev * aver_vright  
    ql = -gammal * vprev * aver_vleft  \
        + np.sqrt(2 * templ * gammal / dt) * xil * aver_vleft
    qr = -gammar * vprev * aver_vright  \
        + np.sqrt(2 * tempr * gammar / dt) * xir * aver_vright
    return ql, qr



if __name__ == '__main__':
    tstart = time.time()

    # x_t = np.zeros((1,N))
    x_t = np.zeros((t_size, N))
    # v_t = np.zeros((1,N))
    vx_t = np.zeros((t_size, N))
    # v_t = np.zeros((t_size, N))
    t_dis = np.array([])
    E_t = np.zeros(t_size)

    x_n_old = np.zeros(N)
    vx_n_old = np.zeros(N)
    ax_n_old = np.zeros(N)
    x_n_new = 0.1 * np.random.normal(0, 1, N)
    vx_n_new = np.random.normal(0, 1, N)
    ax_n_new = np.zeros(N)
    powerL = 0.0
    powerR = 0.0
    power12 = 0.0
    Kin = 0.0
    Kinhalf = 0.0
    Pot = 0.0
    Pothalf = 0.0
    temperature = np.zeros(N)

    tstep = 0
    while tstep < (t_size):
        ##----------not for production, for tests---------
        x2d = x_n_new.reshape(1, N)
        v2d = vx_n_new.reshape(1, N)
        x_t[tstep, :] = x2d
        vx_t[tstep, :] = x2d
        # vx2d = vx_n_new.reshape(1, N)
        # vx_t[tstep, :] = vx2d
        ##------------------------------------------------

        if tstep % 100000 == 0:
        # if tstep % 1000000 == 0:
            print(tstep)


        xil = random.gauss(0, 1)
        xir = random.gauss(0, 1)

        x_n_old = x_n_new
        vx_n_old = vx_n_new
        ax_n_old = ax_n_new
        x_n_new, vx_n_new, ax_n_new = (vv_update(x_n_old,vx_n_old, ax_n_old,
            deltat, tstep))

        if tstep > half_tsize:
            powerLtempx, powerRtempx = heat_update(vx_n_old, vx_n_new, deltat)
            powerL += powerLtempx 
            powerR += powerRtempx 
            Kinhalf += 0.5 * vx_n_new**2
            Pothalf += 0.5 * x_n_new**2
            E_t[tstep-1] = (Kinhalf + Pothalf) / (tstep + 1 - half_tsize)
        # Kin += 0.5 * vx_n_new**2
        # Pot += 0.5 * x_n_new**2
        tstep += 1

        # E_t[tstep-1] = 0.5 * (vx_n_new**2 + x_n_new**2)
        # E_t[tstep-1] = (Kin + Pot) / tstep
        
    # print(powerL / half_tsize, powerR / half_tsize)
    # print(powerL / half_tsize, powerR / half_tsize)
    Eeq = (Kinhalf + Pothalf) / half_tsize
    # print(np.mean(Teff(omegan, tempr)), np.mean(np.sqrt(Teff(omegan,
    #     tempr)/np.pi)))
    # print(0.5 * np.sum(np.sqrt(Teff(omegan,tempr)*deltaomega*2/np.pi)))
    # print(1/(np.exp(1/tempr)-1))
    # print(np.sum(Teff_list)/Nslice)



    # with open (str(traj) + "traj-" +str(tEnd) + "EndTime-" + time.strftime('-%m%d-%H%M%S.txt'),"w") as outputfile:
    # with open ("omega-" +str(omega_drive) + time.strftime('-%m%d-%H%M%S.txt'),"w") as outputfile:
    with open ("TR-" +str(tempr) + time.strftime('-%m%d-%H%M%S.txt'),"w") as outputfile:
        outputfile.write("This run uses mins {} \n".format((time.time() - \
            start_time)/60))
        outputfile.write("gammal = {} \n".format(gammal))
        outputfile.write("gammar = {} \n".format(gammar))
        outputfile.write("time step = {} \n".format(deltat))
        outputfile.write("mode number = {} \n".format(Nslice))
        outputfile.write("omegaN = {} \n".format(omegaN))
        outputfile.write("Eeq = {} \n".format(Eeq))
    # ##-------------------------------


