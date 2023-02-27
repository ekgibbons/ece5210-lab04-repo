import numpy as np
from scipy import signal
from scipy import io as spio

NUM_TAPS = 900
L = 100

print("your code is being tested with L=%i and NUM_TAPS=%i" %
      (L, NUM_TAPS))

h_aa = signal.firwin(NUM_TAPS, 1/L).astype(np.float32)

x = np.random.uniform(-15000,15000,size=300).astype(np.float32)

x_e = np.zeros(x.size*L,
               dtype=np.float32)
x_e[::L] = x
y = L*np.convolve(x_e,h_aa)[:x_e.size]

spio.mmwrite("x_out",x[:,None])
spio.mmwrite("y_out",y[:,None])


