##################################################################################
# CFD 2 Example scripts 1
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Sine Function Initial Condition
L = 1
a = 1
dx = L/40
x = np.linspace(0,1,40)
time = np.linspace(0,20,40)
u_o = (1+np.sin((2*np.pi*x)/L))/2
c = [0.1, 0.5, 0.9]
# Exact solution:
#u = u_o(x-a*t)
# Finite difference method using Forward Time, First Order Upwind FD
# i = Timestep
# n = location (node). In this case the x axis.
# Iterate
u = np.zeros((40,41))
u[:,0] = u_o
for i, counting in enumerate(time):
    for n, counting2 in enumerate(x):
        u[i,n+1] = u[i,n] - c[2]*(u[i,n] - u[i,n-1])

plt.figure(1)
for i, counter in enumerate(time):
    plt.plot(x,u[:,i])

plt.show()

test = "Test"











