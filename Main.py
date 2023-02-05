import math as M
import numpy as np
from matplotlib import pyplot as plt, animation
from matplotlib.offsetbox import AnchoredText

from C import C

cc = C()
KE = 500
freq = 1E12

freqp = np.zeros(3)
freqp[0] = 100.E6
freqp[1] = 200.E6
freqp[2] = 500.E6

highest_er = 50  # ~biological tissue er for 500 MHz
ddx = ((cc.c0 / M.sqrt(highest_er)) / freq) / 10  # Cells Size
dt = ddx / 6E8  # Time Steps
T = 0
nsteps = 5000

arg = np.zeros(3)
arg[0] = 2 * M.pi * freqp[0] * dt
arg[1] = 2 * M.pi * freqp[1] * dt
arg[2] = 2 * M.pi * freqp[2] * dt

all_steps = np.linspace(0, KE - 1, KE)
ex_low_m2, ex_low_m1, ex_high_m2, ex_high_m1 = 0, 0, 0, 0

ga = np.ones(KE)
gb = np.zeros(KE)
dx = np.zeros(KE)
ix = np.zeros(KE)
ex = np.zeros(KE)
hy = np.zeros(KE)
cb = np.zeros(KE)
ca = np.zeros(KE)

# Fourier transform arguments
m = 5
mag = np.zeros(KE)
real_pt = np.zeros([m, KE])
imag_pt = np.zeros([m, KE])
amp = np.zeros([m, KE])
phase = np.zeros([m, KE])
real_in = np.zeros(m)
imag_in = np.zeros(m)
amp_in = np.zeros(m)
phase_in = np.zeros(m)

kc = round(KE / 3)
kp = round(KE / 1.5)
spread = 12
epsz = 8.8E-12
t0 = 40.

# Lossy Field Eq
# eaf = dt * cc.sigmaSi / (2 * epsz * cc.epislon_r(cc.nSi))
# ca[0:kc] = 1.
# cb[0:kc] = 0.5 / cc.epislon_r(cc.nAir)
# ca[kc:kp] = (1. - eaf) / (1 + eaf)
# cb[kc:kp] = 0.5 / (cc.epislon_r(cc.nSi) * (1 + eaf))
# ca[kp:KE] = 1.
# cb[kp:KE] = 0.5 / cc.epislon_r(cc.nAir)

# Lossy field Eq with Flix density
ga[kc:kp] = 1. / (cc.epislon_r(cc.nSi) + (cc.sigmaSi * dt / epsz))
gb[kc:kp] = cc.sigmaSi * dt / epsz

fig = plt.figure()
grid = plt.GridSpec(2, 1, wspace=6, hspace=0.6)
ax = fig.add_subplot(grid[:, :])
#ay = fig.add_subplot(grid[1, 0])
# Cyclic Number of image snapping
frame_interval = 10
ims = []
for n in range(1, nsteps):
    T = T + 1
    # Calculate dx field
    # Calculate ex field
    for k in range(1, KE):
        dx[k] = dx[k] + 0.5 * (hy[k - 1] - hy[k])  # DX field
        # ex[k] = ca[k] * ex[k] + cb[k] * (hy[k - 1] - hy[k]) #EX field

    # pulse generation gauss and sin
    pulse = M.exp(-.5 * (pow((t0 - T) / spread, 2.0)))
    # pulse = m.sin(2 * m.pi * freq * dt * T)
    dx[5] = dx[5] + pulse
    # ex[5] = ex[5] + pulse

    # Calc ex field from dx field
    for k in range(0, KE - 1):
        ex[k] = ga[k] * (dx[k] - ix[k])
        ix[k] = ix[k] + gb[k] * ex[k]

    # Calc Fourier transform of EX
    for k in range(0, KE):
        mag[k] = mag[k] + pow(ex[k], 2)
        for m in range(0, m - 2):
            real_pt[m][k] = real_pt[m][k] + M.cos(arg[m] * T) * ex[k]
            imag_pt[m][k] = imag_pt[m][k] - M.sin(arg[m] * T) * ex[k]

    # Fourier Transform of the input pulse signal
    if T < 100:
        for m in range(0, 2):
            real_in[m] = real_in[m] + M.cos(arg[m] * T) * ex[10]
            imag_in[m] = imag_in[m] - M.sin(arg[m] * T) * ex[10]

    # Absorbing Boundary Conditions - left end
    ex[0] = ex_low_m2
    ex_low_m2 = ex_low_m1
    ex_low_m1 = ex[1]
    # right end
    ex[-1] = ex_high_m2
    ex_high_m2 = ex_high_m1
    ex_high_m1 = ex[-2]

    # calculate hy field
    for k in range(0, KE - 1):
        hy[k] = hy[k] + 0.5 * (ex[k] - ex[k + 1])  # HY field

    # Calc amp and phase for each freq
    for m in range(0, 2):
        amp_in[m] = M.sqrt(pow(imag_in[m], 2)) + pow(real_in[m], 2)
        phase_in[m] = M.atan2(imag_in[m], real_in[m])
        for k in range(1, KE):
            try:
                amp[m][k] = (1. / amp_in[m]) * M.sqrt(pow(real_pt[m][k], 2)) + pow(imag_pt[m][k], 2)
            except ZeroDivisionError:
                amp[m][k] = 0.
            phase[m][k] = M.atan2(imag_pt[m][k], real_pt[m][k]) - phase_in[m]

    # Drawing of the EM and FT plots
    if T % frame_interval == 0:
        anchored_text = AnchoredText("EX and EM field after time t=" + '{:<.3e}'.format(T), loc=1, prop=dict(size=8))
        title = ax.add_artist(anchored_text)
        ims1, = ax.plot(all_steps, ex, "r")
        ims2, = ax.plot(all_steps, hy, "b--")
        #ims3, = ay.plot(mag[:],"g")
        ims.append([ims1, ims2, title])

ax.axvspan(kc, kp, facecolor='green', alpha=0.3)
file_name = "1d_fdtd_lossy2"
# file_name = "./" + file_name + '.gif'
file_name = "./" + file_name + '.gif'
ani = animation.ArtistAnimation(fig, ims, interval=10, blit=True)
#ani.save(file_name, writer='pillow', fps=25, dpi=200)
# ani.save(file_name, writer="imagemagick", fps=30)
print("OK")
plt.show()
