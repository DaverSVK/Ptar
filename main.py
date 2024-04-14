import math
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import tso
import etso
import etsoLimited
import tsoLimited

import control.matlab as matlab

def tsoNoLimitations(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext):
    t_log, x_log,e_log,u_log= tso.fcn_simSch_max(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )       
    t_log2, x_log2,e_log2,u_log2 = tso.fcn_simSch_min(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )
    printGraphs(t_log,t_log2,x_log,x_log2,u_log,u_log2,e_log,e_log2, 'TSO bez obmedzenia rýchlosti','tso_nolim')

def tsoSpeedLimitations(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext):
    t_log, x_log,e_log,u_log= tsoLimited.fcn_simSch_max(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )       
    t_log2, x_log2,e_log2,u_log2 = tsoLimited.fcn_simSch_min(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )
    printGraphs(t_log,t_log2,x_log,x_log2,u_log,u_log2,e_log,e_log2, 'TSO','tso_lim')

def etsoNoLim(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext):
    t_log, x_log,e_log,u_log= etso.fcn_simSch_max(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )       
    t_log2, x_log2,e_log2,u_log2 = etso.fcn_simSch_min(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )
    printGraphs(t_log,t_log2,x_log,x_log2,u_log,u_log2,e_log,e_log2, 'ETSO','etso_no_lim')

def etsoLim(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext):
    t_log, x_log,e_log,u_log= etsoLimited.fcn_simSch_max(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )       
    t_log2, x_log2,e_log2,u_log2 = etsoLimited.fcn_simSch_min(
                        sim_t_start,
                        sim_T_s,
                        sim_finalIndex,
                        sig_u_ext,
                        )
    printGraphs(t_log,t_log2,x_log,x_log2,u_log,u_log2,e_log,e_log2, 'ETSO Lim','etso_lim')

def printGraphs(t_log,t_log2,x_log,x_log2,u_log,u_log2,e_log,e_log2,addTit,saveName):
    plt.figure(figsize=(8, 4))
    plt.plot(t_log, x_log[:,0])
    plt.plot(t_log2, x_log2[:,0])
    plt.xlabel('čas')
    plt.ylabel('x [rad]')
    plt.legend(['Tmax','Tmin'])
    plt.title('Poloha systému '+ addTit)
    plt.savefig('poloha_' + saveName + '.png')
    plt.grid()

    plt.figure(figsize=(8, 4))
    plt.plot(t_log, u_log[:,0])
    plt.plot(t_log2, u_log2[:,0])
    plt.xlabel('čas')
    plt.ylabel('u [Nm]')
    plt.legend(['Tmax','Tmin'])
    plt.title('Akčný zásah '+ addTit)
    plt.savefig('akcnyZ_' + saveName + '.png')
    plt.grid()

    plt.figure(figsize=(8, 4))
    plt.plot(t_log, x_log[:,1])
    plt.plot(t_log2, x_log2[:,1])
    plt.xlabel('čas')
    plt.ylabel('v [rad/s]')
    plt.legend(['Tmax','Tmin'])
    plt.title('Rýchlosť systému '+ addTit)
    plt.savefig('rychlost_' + saveName + '.png')
    plt.grid()

    plt.figure(figsize=(8, 4))
    plt.plot(e_log,x_log[:,1]*-1)
    plt.plot(e_log2,x_log2[:,1]*-1)
    plt.xlabel('e [rad]')
    plt.ylabel('e\' [rad/s]')
    plt.legend(['Tmax','Tmin'])
    plt.title('Fázový graf '+ addTit)
    plt.savefig('faza_' + saveName + '.png')
    plt.grid()


#---------------------------------------------------------------------------

#  Main 

#---------------------------------------------------------------------------

sim_t_start = 0
sim_t_final = 1
sim_T_s = 0.0001
sim_finalIndex = int(((sim_t_final - sim_t_start)/sim_T_s) + 1) ### cellE c03 ###

# ---------------------------------------------------------------------------
sig_delt_u = np.zeros([sim_finalIndex, 1])
for idx in range(sig_delt_u.shape[0]):
    sig_delt_u[idx] = 1


sig_u_ext = sig_delt_u

tsoSpeedLimitations(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext)
tsoNoLimitations(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext)
etsoNoLim(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext)
etsoLim(sim_t_start,sim_T_s,sim_finalIndex,sig_u_ext)
# # Show the plot
plt.show()
