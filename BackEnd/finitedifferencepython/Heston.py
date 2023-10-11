import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import spdiags, csc_array
from scipy.sparse.linalg import inv

from scipy.sparse.linalg import spsolve

np.set_printoptions(precision=6, suppress=True)

dt = 1/4000
T = 1
S = 10.0
K = 5.0
I = 20
J = 50
V = 1.0
r = 0.15
sigma = 0.2
rho = 0.8
kappa = 2.0
eta = 0.2
theta = 0.5

nt = T / dt
# dt = T/nt

sbound = S * 1.2
vbound = V * 1.3
ds = S / I  # step length of s
dv = V / J  # step length of v
svalue = np.arange(0, sbound, ds)
vvalue = np.arange(0, vbound, dv)
ns = svalue.size - 1
nv = vvalue.size - 1
iline = np.linspace(0, ns, ns + 1, dtype = np.float64)
jline = np.linspace(0, nv, nv + 1, dtype = np.float64)


def setM(M, i0, j0, i1, j1, v):
    M[i0*(nv+1)+j0, i1*(nv+1)+j1] = v

def getM(M, i0, j0, i1, j1):
    return M[i0*(nv+1)+j0, i1*(nv+1)+j1]

M0 = np.zeros(((ns+1)*(nv+1), (ns+1)*(nv+1)))
for i in range(1, ns):
    for j in range(1, nv):
        v = rho * sigma * iline[i] * jline[j] / 4
        setM(M0, i, j, i-1, j-1, v)
        setM(M0, i, j, i-1, j+1, -v)
        setM(M0, i, j, i+1, j-1, -v)
        setM(M0, i, j, i+1, j+1, v)

#CSC format
M0 = csc_array(M0)

M1 = np.zeros(((ns+1)*(nv+1), (ns+1)*(nv+1)))
for i in range(1, ns):
    for j in range(1, nv):
        setM(M1, i, j, i-1, j, iline[i] * iline[i] * vvalue[j] / 2 - r * iline[i] / 2)
        setM(M1, i, j, i, j, -iline[i] * iline[i] * vvalue[j] - r/2)
        setM(M1, i, j, i+1, j, iline[i] * iline[i] * vvalue[j] / 2 + r * iline[i] / 2)
M1 = csc_array(M1)

M2 = np.zeros(((ns+1)*(nv+1), (ns+1)*(nv+1)))
for i in range(1, ns):
    for j in range(1,nv):
        setM(M2, i, j, i, j-1, sigma * sigma * jline[j] / (2 * dv) - kappa * (eta - vvalue[j]) / (2 * dv))
        setM(M2, i, j, i, j, -sigma * sigma * jline[j] / dv - r / 2)
        setM(M2, i, j, i, j+1, sigma * sigma * jline[j] / (2 * dv) + kappa * (eta - vvalue[j]) / (2 * dv))
M2 = csc_array(M2)

M = M0 + M1 + M2

M3 = np.zeros(((ns+1)*(nv+1), (ns+1)*(nv+1)))
for i in range(1, ns):
    setM(M3, i, 0, i+1, 0, r * i)
    setM(M3, i, 0, i, 0, -r * i - kappa * eta / dv - r)
    setM(M3, i, 0, i, 1, kappa * eta / dv)

M3 = csc_array(M3)


MI = np.identity((ns+1)*(nv+1))
MI = csc_array(MI)

def DR(U, dt):
    # Y = spsolve(MI - dt * theta * M1, (MI + dt * M0 + (1-theta)*dt*M1 + dt*M2) @ U)
    # return spsolve(MI - theta * dt*M2, Y - theta * dt*M2 @ U)
    Y1 = spsolve(MI - dt * theta * M1, (MI + dt * M0 + (1-theta)*dt*M1 + dt*M2) @ U)
    Y2 = spsolve(MI - theta * dt*M2, Y1 - theta * dt*M2 @ U)
    Y3 = spsolve(MI - theta * dt*M1, (MI-theta*dt*M1)@Y1+0.5*dt*M0@(Y2-U))
    return spsolve(MI - theta * dt*M2, Y3 - theta * dt*M2 @ U)


def run(dt):
    temp = np.maximum(0, svalue - K)
    uinitial = np.outer(temp, np.ones(nv + 1, dtype = np.float64))
    U = uinitial
    for n in range(int(1/dt)):
        if n % (int(0.1/dt)) == 0:
            print('n = ', n)
        U = np.reshape(U, (ns+1)*(nv+1))
        # U = U + M @ U
        U = DR(U, dt)
        U = np.reshape(U, (ns+1, nv+1))

        # S = 0: U = 0
        U[0, :] = np.zeros(nv + 1, dtype = np.float64)

        # V = \infty: U = S
        U[:, nv] = svalue

        U = np.reshape(U, (ns+1)*(nv+1))
        # V = 0 : TODO: type equation
        U = U + dt * M3 @ U
        U = np.reshape(U, (ns+1, nv+1))

        # S = \infty: dU/dS = 1
        U[ns, 1:nv] = (2 * ds + 4 * U[ns - 1, 1:nv] - U[ns - 2, 1:nv]) / 3


    return U[0:I + 1, 0:J + 1]

U = run(1/4000)
s = np.linspace(0, S, I + 1)
v = np.linspace(0, V, J + 1)
[Sval, Vval] = np.meshgrid(v, s)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
surf = ax.plot_surface(Sval, Vval, U)
ax.set_xlabel('Volatility (V)')
ax.set_ylabel('Underlying Asset Price (S)')
ax.set_zlabel('Option Price (U)')
ax.set_title('Option Price Surface at time t = 0')
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

