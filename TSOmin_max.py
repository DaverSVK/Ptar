import math
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

import control.matlab as matlab
K = 0.0883
# časová konštanta 
Tmin = 0.019 
Tmax = 0.068
#  obmedzenie akčnej veličiny 
M = 50
#  žiadaná hodnota 
w = 1
#  parameter riadenia 
k = 100
Tp = 0.001
bp = 2 
D=M/2
M1=M+D
M2 = M-D 
vmax = K*M
vlim = 0.5*vmax
caseTSO=1
alpha = 1/(Tmax*(1-(M2/M1)*math.log(1+(M1/M2))))
print(alpha)


a = -1
b = K    

def fcn_difRovniceMax(x, t, u):

    dotx = (a*x + u*b)/Tmax
    return dotx 

def fcn_difRovniceMin(x, t, u):

    dotx = (a*x + u*b)/Tmin
    return dotx 

def fcn_difRovnice1(x, t, u):

    dotx = u
    return dotx 



def fcn_simSch_max(t_start, T_s, finalIndex, sig_u_ext):

    #-----------------------------------------
    t_log = np.zeros([finalIndex, 1])
    t_log[0,:] = t_start

    #-----------------------------------------
    x_0 = 0

    x_log = np.zeros([finalIndex, 1])
    e_log = np.zeros([finalIndex, 1])
    x_log[0,:] = x_0
    x_logNorm = np.zeros([finalIndex, 1])
    x_logNorm[0,:] = x_0

    #-----------------------------------------


    #-----------------------------------------
    timespan = np.zeros(2)
    for idx in range(1, int(finalIndex)):

        timespan[0] = t_log[idx-1,:]
        timespan[1] = t_log[idx-1,:] + T_s
        e_log[idx-1,:] = (sig_u_ext[idx-1,:]-x_logNorm[idx-1,:])
        signum = ((sig_u_ext[idx-1,:]-x_logNorm[idx-1,:])*alpha)+(x_log[idx-1,:]*-1)
        if signum > 0:
            u=1*M
        elif signum < 0:
            u=-1*M
        else:
            u=0
        odeOut = odeint(fcn_difRovniceMax,
                        x_log[idx-1,:],
                        timespan,
                        args=(u,)
                        )
        x_log[idx,:] = odeOut[-1,:]
        t_log[idx,:] = timespan[-1]

        odeOut2 = odeint(fcn_difRovnice1,
                        x_logNorm[idx-1,:],
                        timespan,
                        args=(x_log[idx-1,:],)
                        )

        x_logNorm[idx,:] = odeOut2[-1,:]

    return [t_log, x_log,x_logNorm,e_log]                


def fcn_difRovnice(x, t, u):

    dotx = (a*x + u*b)/Tmax
    return dotx 
def fcn_difRovnice1(x, t, u):

    dotx = u
    return dotx 


def fcn_simSch_min(t_start, T_s, finalIndex, sig_u_ext):

    #-----------------------------------------
    t_log = np.zeros([finalIndex, 1])
    t_log[0,:] = t_start

    #-----------------------------------------
    x_0 = 0

    x_log = np.zeros([finalIndex, 1])
    e_log = np.zeros([finalIndex, 1])
    x_log[0,:] = x_0
    x_logNorm = np.zeros([finalIndex, 1])
    x_logNorm[0,:] = x_0

    #-----------------------------------------
    timespan = np.zeros(2)
    for idx in range(1, int(finalIndex)):

        timespan[0] = t_log[idx-1,:]
        timespan[1] = t_log[idx-1,:] + T_s
        e_log[idx-1,:] = (sig_u_ext[idx-1,:]-x_logNorm[idx-1,:])
        signum = ((sig_u_ext[idx-1,:]-x_logNorm[idx-1,:])*alpha)+(x_log[idx-1,:]*-1)
        if signum > 0:
            u=1*M
        elif signum < 0:
            u=-1*M
        else:
            u=0
        odeOut = odeint(fcn_difRovniceMin,
                        x_log[idx-1,:],
                        timespan,
                        args=(u,)
                        )
        x_log[idx,:] = odeOut[-1,:]
        t_log[idx,:] = timespan[-1]

        odeOut2 = odeint(fcn_difRovnice1,
                        x_logNorm[idx-1,:],
                        timespan,
                        args=(x_log[idx-1,:],)
                        )

        x_logNorm[idx,:] = odeOut2[-1,:]

    return [t_log, x_log,x_logNorm,e_log]                

sim_t_start = 0
sim_t_final = 1
sim_T_s = 0.0001
sim_finalIndex = int(((sim_t_final - sim_t_start)/sim_T_s) + 1) ### cellE c03 ###

# ---------------------------------------------------------------------------
                                                                ### cellB c04 ###

sig_delt_u = np.zeros([sim_finalIndex, 1])
for idx in range(sig_delt_u.shape[0]):
    sig_delt_u[idx] = 1


sig_u_ext = sig_delt_u

t_log, x_log, x_logNorm,e_log = fcn_simSch_max(
                    sim_t_start,
                    sim_T_s,
                    sim_finalIndex,
                    sig_u_ext,
                    )       
t_log2, x_log2, x_logNorm2,e_log2 = fcn_simSch_min(
                    sim_t_start,
                    sim_T_s,
                    sim_finalIndex,
                    sig_u_ext,
                    )

plt.figure()
plt.plot(t_log, x_log)
plt.plot(t_log2, x_log2)
plt.xlabel('Time')
plt.ylabel('X')
plt.title('Plot of X vs Time')
plt.grid()

plt.figure()
plt.plot(t_log, x_logNorm)
plt.plot(t_log2, x_logNorm2)
plt.xlabel('Time')
plt.ylabel('X')
plt.title('Plot of X vs Time')
plt.grid()

plt.figure()
plt.plot(e_log,x_log*-1)
plt.plot(e_log2,x_log2*-1)
plt.xlabel('Time')
plt.ylabel('X')
plt.title('Plot of X vs Time')
plt.grid()
# Show the plot
plt.show()

