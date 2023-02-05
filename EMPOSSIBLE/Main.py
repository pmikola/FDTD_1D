# Main FDTD LOOP
import numpy as np
import math as m
import matplotlib.pyplot as plt
from time import sleep

from matplotlib.offsetbox import AnchoredText
from scipy import constants as const
from IPython.core.display import HTML
from celluloid import Camera
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.animation as animation
from PIL import Image, ImageDraw
from C import C
from Source import source, source_sin
import matplotlib

matplotlib.use('Qt5Agg')
size = 500
Z = size - 1
cc = C()
nAir = cc.nAir
nGe = cc.nGe
nSi = cc.nSi
epsR = []
muR = []
frequency = 10*10**1
z1 = round(size / 2 - size / 10)
z2 = round(size / 2 + size / 10)
for i in range(0, size):
    if i < z1:
        epsR.append(m.sqrt(nAir))
        muR.append(1)
    elif z1 <= i <= z2:
        epsR.append(m.sqrt(nSi))
        muR.append(1)
    else:
        epsR.append(m.sqrt(nAir))
        muR.append(1)

epsilon = cc.epislon_r(nGe)
# eps = np.ones(size)
eps = epsR

# Source
source_width = 25. * np.sqrt(epsilon)
# source_width = size*np.sqrt(epsilon)
fmax = (1 / m.pi * source_width) * 0.1
if fmax >= frequency:
    pass
else:
    fmax = frequency
deltalambda = cc.gridRes(fmax)
Nt = 6
Nlambda = 60
Nd = 6
dmin = 4
lambdaMin = cc.c0 / fmax * nGe
delay = 10. * source_width
source_z = round(size * 0.1)  # Source position

# MODEL
steps = 10000  # int((size + delay) * np.sqrt(epsilon))  # Time stepping
frame_interval = 20
all_steps = np.linspace(0, size - 1, size)
VoidImp = cc.Impedance(1, 1)
Ey = np.zeros(size)
Hx = np.zeros(size)
mEy = np.zeros(size)
mHx = np.zeros(size)
delta1 = cc.deltat(cc.c0, 1, 1, nSi)
delta2 = (1 / m.pi * fmax) / Nt
delta3 = deltalambda
dt = min([delta1, delta2, delta3])
# dz1 = lambdaMin / Nlambda
# dz2 = dmin / Nd
# dz = min([dz1, dz2])
dzprim = min([lambdaMin / Nlambda, dmin / Nd])
dz = dzprim * 5
while dz > cc.c0 * dt * 10 / nAir:
    print(dz)
    dz -= dz * 0.001
    # print(dz)
delta4 = nAir * dz / 2 * cc.c0

if delta4 > dt:
    pass
else:
    dt = delta4

mHx[:] = cc.c0 * dt
mEy[:] = cc.c0 * dt
fig = plt.figure()
grid = plt.GridSpec(1, 1, wspace=6, hspace=0.6)
ax = fig.add_subplot(grid[:, :])
ims = []
for time in range(0, steps):
    # plt.clf()
    Hx[-1] = Hx[-2]  # simple ABC for Hx
    Hx[0] = Hx[1]  # simple ABC for Hx
    for z in range(0, Z):
        # Hx field from derivation of Ey field with normalization VoidImp and Coefficients
        Hx[z] = Hx[z] + (mHx[z] / muR[z] * (Ey[z + 1] - Ey[z]) / dz) / VoidImp / muR[z]
    #Hx[source_z - 1] -= source("Hx", time, delay, source_width, 100.) / VoidImp
    # Hx[Z] = Hx[Z] + (mHx[Z] * (0 - Ey[Z]) / dz) / VoidImp
    # Ey[0] = Ey[0] + mEy[0] * ((Hx[0] - 0) / dz) * VoidImp
    Ey[0] = Ey[1]  # simple ABC for Ey
    Ey[-1] = Ey[-2]  # simple ABC for Ey
    for z in range(1, Z):
        # Ey field from derivation of Hx field with normalization VoidImp and Coefficients
        Ey[z] = Ey[z] + (mEy[z] / epsR[z] * (Hx[z] - Hx[z - 1]) / dz) * VoidImp / epsR[z]
    Ey[source_z] += source("Ey", time, delay, source_width, 100.)
    #Ey[source_z] += source_sin(frequency, source_width, dt)
    # Hx[source_z] += cc.Impedance(epsilon, 1) * source(time, delay, source_width) * dz
    if time % frame_interval == 0:  # snapshots
    # ax.title("EH after time t=%i" % time)
        anchored_text = AnchoredText("EH after time t=" + '{:<.3e}'.format(time), loc=1, prop=dict(size=8))
        title = ax.add_artist(anchored_text)
        ims1, = ax.plot(all_steps, Ey, "r")
        ims2, = ax.plot(all_steps, Hx * cc.Impedance(1, 1), "b-")
        # ax.set_ylim(-max(Ey), max(Ey))
        # plt.grid()
        # ims, = [plt.plot(all_steps, Ey, all_steps, Hx * cc.Impedance(1, 1))]
        # plt.pause(0.001)
        #print(time)
        ims.append([ims1, ims2, title])
        # plt.clf()

plt.axvline(0, color='black')
plt.axvline(size, color='black')
plt.axvspan(z1, z2, facecolor='green', alpha=0.3)
# figManager = plt.get_current_fig_manager()
# figManager.window.resize(700, 700)
file_name = "1d_fdtd"
# file_name = "./" + file_name + '.gif'
file_name = "./" + file_name + '.gif'
ani = animation.ArtistAnimation(fig, ims, interval=10, blit=True)
#ani.save(file_name, writer='pillow', fps=25, dpi=100)
# ani.save(file_name, writer="imagemagick", fps=30)
print("OK")
plt.show()