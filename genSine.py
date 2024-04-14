
 # Define sine wave parameters
# Amplitude of the sine wave
# Frequency of the sine wave (you can adjust this value)
# Phase shift of the sine wave (optional, you can adjust this value)
import numpy as np
def genSine(amplitude,frequency,phase):
    sim_t_start = 0
    sim_t_final = 1
    sim_T_s = 0.0001
    sim_finalIndex = int(((sim_t_final - sim_t_start) / sim_T_s) + 1)

    sig_vysl1 = np.zeros([sim_finalIndex, 1])

    for idx in range(sim_finalIndex):
        # Calculate the time corresponding to the current index
        t = sim_t_start + idx * sim_T_s
        # Generate the corresponding value of the sine wave at time t
        sig_vysl1[idx] = amplitude * np.sin(frequency * t + phase)

    return sig_vysl1