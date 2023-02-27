import numpy as np
from scipy import signal
from scipy import io as spio

NUM_TAPS = 900
M = 100
POLY_TAPS = NUM_TAPS//M

h_aa = signal.firwin(NUM_TAPS, 1/M).astype(np.float32)

h = np.zeros_like(h_aa)


for k in range(M):
    h[k*POLY_TAPS:(k+1)*POLY_TAPS] = h_aa[k::M]
    
spio.mmwrite("h_out",h[:,None])
