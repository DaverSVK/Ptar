import numpy as np
import math
from scipy.integrate import odeint
import genSine
# zosilnenie
K = 0.0883
# časová konštanta 
Tmin = 0.019 
Tmax = 0.068
#  obmedzenie akčnej veličiny 
M = 50
#  žiadaná hodnota 
w = 1
#  parameter riadenia 
k = 1000
Tp = 0.001
bp = 2 
D=M/2
M1=M+D
M2 = M-D 
vmax = K*M
vlim = 0.5*vmax
alpha = 1/(Tmax*(1-(M2/M1)*math.log(1+(M1/M2))))


def fcn_difRovniceMax(x, t, u):
    A = np.array([[0, 1], [-0, -1/Tmax]])
    b = np.array([[0], [K/Tmax]])    
    dotx = np.dot(A,x) + np.dot(b,u)
    return dotx 

def fcn_difRovniceMin(x, t, u):
    A = np.array([[0, 1], [-0, -1/Tmin]])
    b = np.array([[0], [K/Tmin]])    
    dotx = np.dot(A,x) + np.dot(b,u)
    return dotx 

def fcn_simSch_max(t_start, T_s, finalIndex, sig_u_ext):

    t_log = np.zeros([finalIndex, 1])
    t_log[0,:] = t_start

    #-----------------------------------------
    x_0 = np.array([0, 0])

    x_log = np.zeros([finalIndex, len(x_0)])
    e_log = np.zeros([finalIndex, 1])
    x_log[0,:] = x_0
    u_log = np.zeros([finalIndex, 1])
    sine_sig = genSine.genSine(D,10,0)
    #-----------------------------------------
    timespan = np.zeros(2)
    for idx in range(1, int(finalIndex)):

        timespan[0] = t_log[idx-1,:]
        timespan[1] = t_log[idx-1,:] + T_s
        e_log[idx-1,:] = (sig_u_ext[idx-1,:]-x_log[idx-1,0])
        e_base = (sig_u_ext[idx-1,:]-x_log[idx-1,0])*k
        e_derivated = x_log[idx-1,1]*-1
        alpha_e_and_ed = (e_base + e_derivated)*alpha*Tmax
        if alpha_e_and_ed > k*Tmax*vlim:
            alpha_e_and_ed=1*(k*Tmax*vlim)
        elif alpha_e_and_ed < -(k*Tmax*vlim):
            alpha_e_and_ed=-1*(k*Tmax*vlim)
        else:
            alpha_e_and_ed=alpha_e_and_ed

        aplpha_e_derivated = e_derivated*(k*Tmax-1)
        saturation = (alpha_e_and_ed+aplpha_e_derivated)/K;

        if saturation > M:
            u_log[idx-1,:]=1*M
        elif saturation < -M:
            u_log[idx-1,:]=-1*M
        else:
            u_log[idx-1,:]=saturation
        u_log[idx-1,:]=u_log[idx-1,:]-sine_sig[idx-1,:]
        odeOut = odeint(fcn_difRovniceMax,
                        x_log[idx-1,:],
                        timespan,
                        args=(u_log[idx-1,:],)
                        )
        x_log[idx,:] = odeOut[-1,:]
        t_log[idx,:] = timespan[-1]

    return [t_log, x_log,e_log,u_log]   

def fcn_simSch_min(t_start, T_s, finalIndex, sig_u_ext):

    #-----------------------------------------
    t_log = np.zeros([finalIndex, 1])
    t_log[0,:] = t_start

    #-----------------------------------------
    x_0 = np.array([0, 0])

    x_log = np.zeros([finalIndex, len(x_0)])
    e_log = np.zeros([finalIndex, 1])
    x_log[0,:] = x_0
    u_log = np.zeros([finalIndex, 1])
    sine_sig = genSine.genSine(D,10,0)
    #-----------------------------------------
    timespan = np.zeros(2)
    for idx in range(1, int(finalIndex)):

        timespan[0] = t_log[idx-1,:]
        timespan[1] = t_log[idx-1,:] + T_s
        e_log[idx-1,:] = (sig_u_ext[idx-1,:]-x_log[idx-1,0])
        e_base = (sig_u_ext[idx-1,:]-x_log[idx-1,0])*k
        e_derivated = x_log[idx-1,1]*-1
        alpha_e_and_ed = (e_base + e_derivated)*alpha*Tmax
        if alpha_e_and_ed > k*Tmax*vlim:
            alpha_e_and_ed=1*(k*Tmax*vlim)
        elif alpha_e_and_ed < -(k*Tmax*vlim):
            alpha_e_and_ed=-1*(k*Tmax*vlim)
        else:
            alpha_e_and_ed=alpha_e_and_ed

        aplpha_e_derivated = e_derivated*(k*Tmax-1)
        saturation = (alpha_e_and_ed+aplpha_e_derivated)/K;
        
        if saturation > M:
            u_log[idx-1,:]=1*M
        elif saturation < -M:
            u_log[idx-1,:]=-1*M
        else:
            u_log[idx-1,:]=saturation
        u_log[idx-1,:]=u_log[idx-1,:]-sine_sig[idx-1,:]
        odeOut = odeint(fcn_difRovniceMin,
                        x_log[idx-1,:],
                        timespan,
                        args=(u_log[idx-1,:],)
                        )
        x_log[idx,:] = odeOut[-1,:]
        t_log[idx,:] = timespan[-1]

    return [t_log, x_log,e_log,u_log]