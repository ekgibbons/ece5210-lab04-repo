import numpy as np
from scipy import signal
from scipy import io as spio

def naive_decimation(x, h, M):
    y = np.convolve(h, x)[:x.size]
    
    return y[::M]

NUM_TAPS = 900
M = 100

print("your code is being tested with M=%i and NUM_TAPS=%i" %
      (M, NUM_TAPS))

h_aa = signal.firwin(NUM_TAPS, 1/M).astype(np.float32)

x = np.random.randint(-15000,15000,size=M*100).astype(np.float32)

y = naive_decimation(x, h_aa, M)

spio.mmwrite("x_out",x[:,None])
spio.mmwrite("y_out",y[:,None])


