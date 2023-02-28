import time
from pyevtk.hl import *
import numpy as np
import taichi as ti
import initialization

# ti.init(arch=ti.cpu, default_fp=ti.f64, cpu_max_num_threads=1, advanced_optimization=False, fast_math=False)
ti.init(arch=ti.cpu, default_fp=ti.f64, advanced_optimization=False, fast_math=False)

DT = 0.0005
END_T = 1.0
OPT_FQC = 100
it = int(END_T / DT)
tn0, tlmd = 0.0, 0.0
PCL_DST = 0.02
MIN_X = (initialization.MINX - PCL_DST * 3)
MIN_Y = (initialization.MINY - PCL_DST * 3)
MIN_Z = (initialization.MINZ - PCL_DST * 3)
MAX_X = (initialization.MAXX + PCL_DST * 3)
MAX_Y = (initialization.MAXY + PCL_DST * 3)
MAX_Z = (initialization.MAXZ + PCL_DST * 30)
GST, FLD, WLL = -1, 0, 1
SND = 22.0
KNM_VSC = 1e-6
DIM = 3
CRT_NUM = 0.1
COL_RAT = 0.2
DST_LMT_RAT = 0.9
G_X, G_Y, G_Z = 0.0, 0.0, -9.8
t = 0
r = PCL_DST * 2.1
r2 = r * r
DB = r * (1.0 + CRT_NUM)
DB2 = DB * DB
DBinv = 1.0 / DB
DNS = 1000
invDNS = 1.0 / DNS
nBx = int(((MAX_X - MIN_X) * DBinv) + 3)
nBy = int(((MAX_Y - MIN_Y) * DBinv) + 3)
nBz = int(((MAX_Z - MIN_Z) * DBinv) + 3)
nBxy = (nBx * nBy)
nBxyz = (nBx * nBy * nBz)
nP = initialization.k
num = np.zeros(nP, dtype=int)
px_p = np.zeros(nP, dtype=float)
py_p = np.zeros(nP, dtype=float)
pz_p = np.zeros(nP, dtype=float)
typ_p = -1 * np.ones(nP, dtype=int)
data = np.loadtxt('int3dpy.txt', skiprows=1)
num[:] = data[:, 0]
typ_p[:] = data[:, 1]
px_p[:] = data[:, 2]
py_p[:] = data[:, 3]
pz_p[:] = data[:, 4]
typ = ti.field(ti.i64, nP)
bfst = ti.field(ti.i64, nBxyz)
blst = ti.field(ti.i64, nBxyz)
nxt = ti.field(ti.i64, nP)
px = ti.field(ti.f64, nP)
py = ti.field(ti.f64, nP)
pz = ti.field(ti.f64, nP)
vx = ti.field(ti.f64, nP)
vy = ti.field(ti.f64, nP)
vz = ti.field(ti.f64, nP)
ax = ti.field(ti.f64, nP)
ay = ti.field(ti.f64, nP)
az = ti.field(ti.f64, nP)
prs = ti.field(ti.f64, nP)
pav = ti.field(ti.f64, nP)

typ.from_numpy(typ_p)
px.from_numpy(px_p)
py.from_numpy(py_p)
pz.from_numpy(pz_p)

for ix in range(-4, 5):
    for iy in range(-4, 5):
        for iz in range(-4, 5):
            x = PCL_DST * ix
            y = PCL_DST * iy
            z = PCL_DST * iz
            dist2 = x * x + y * y + z * z
            if dist2 <= r2:
                if dist2 == 0:
                    continue
                dist = np.sqrt(dist2)
                wei = ((r / dist) - 1.0)
                tn0 += wei
                tlmd += dist2 * wei
n0 = tn0
lmd = tlmd / tn0
A1 = 2.0 * KNM_VSC * DIM / n0 / lmd
A2 = SND * SND / n0
A3 = -DIM / n0
rlim = PCL_DST * DST_LMT_RAT
rlim2 = rlim ** 2
COL = 1.0 + COL_RAT


#
# def Output():
#     positionx = px.to_numpy()
#     positiony = py.to_numpy()
#     positionz = pz.to_numpy()
#     vel_x = vx.to_numpy()
#     vel_y = vy.to_numpy()
#     vel_z = vz.to_numpy()
#     pressave = pav.to_numpy() / OPT_FQC
#     type_out = typ.to_numpy()
#     name = str(format(t))
#     npa = str(nP)
#     Velocity = np.sqrt(vel_x * vel_x + vel_y * vel_y + vel_z * vel_z)
#     pxa = positionx[np.where(type_out != WLL)]
#     pya = positiony[np.where(type_out != WLL)]
#     pza = positionz[np.where(type_out != WLL)]
#     pa = pressave[np.where(type_out != WLL)]
#     uxa = vel_x[np.where(type_out != WLL)]
#     uya = vel_y[np.where(type_out != WLL)]
#     uza = vel_z[np.where(type_out != WLL)]
#     ua = Velocity[np.where(type_out != WLL)]
#     point_data = {"Pressure": pa, "Velocity": ua}
#     pointsToVTK("vtu/./emps" + name, pxa, pya, pza, data=point_data)
#
#     np.savetxt(r"txt/emps" + name + ".prof", np.column_stack((pxa, pya, pza, uxa, uya, uza, pa)),
#                fmt=' %.6f %.6f %.6f %.6f %.6f %.6f %.6f ', header=npa, comments='')
#     pav.from_numpy(pz_p)
def Output():
    positionx = px.to_numpy()
    positiony = py.to_numpy()
    positionz = pz.to_numpy()
    vel_x = vx.to_numpy()
    vel_y = vy.to_numpy()
    vel_z = vz.to_numpy()
    prssure_out = prs.to_numpy()
    pressave = pav.to_numpy() / OPT_FQC
    type_out = typ.to_numpy()
    name = str(format(t))
    npa = str(nP)
    Velocity = np.sqrt(vel_x * vel_x + vel_y * vel_y + vel_z * vel_z)
    point_data = {"Pressure": pressave, "Type": type_out, "Velocity": Velocity}
    pointsToVTK("vtu/./emps" + name, positionx, positiony, positionz, data=point_data)
    np.savetxt(r"txt/emps" + name + ".prof", np.column_stack((num[0:nP], type_out, positionx,
                                                              positiony, positionz,
                                                              vel_x, vel_y,
                                                              vel_z, prssure_out, pressave)),
               fmt='%d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f', header=npa, comments='')
    pav.from_numpy(pz_p)


@ti.func
def Check_position(i):
    if px[i] > MAX_X or px[i] < MIN_X or py[i] > MAX_Y or py[i] < MIN_Y or pz[i] > MAX_Z or pz[i] < MIN_Z:
        typ[i] = GST
        prs[i] = vx[i] = vy[i] = vz[i] = 0.0


@ti.kernel
def Calculate_num():
    p_num = 0
    for x in range(nP):
        if typ[x] != GST:
            p_num += 1
    print("p_num:", p_num)


@ti.kernel
def Reset_bucket():
    for i in nxt:
        nxt[i] = -1
    for i in blst:
        blst[i] = bfst[i] = -1


@ti.kernel
def Assigh_bucket():
    for ii in range(1):
        # ti.loop_config(serialize=True)
        for i in range(nP):
            if typ[i] == GST:
                continue
            ix = int(((px[i] - MIN_X) * DBinv) + 1)
            iy = int(((py[i] - MIN_Y) * DBinv) + 1)
            iz = int(((pz[i] - MIN_Z) * DBinv) + 1)
            ib = int(iz * nBxy + iy * nBx + ix)
            j = int(blst[ib])
            blst[ib] = int(i)
            if j == -1:
                bfst[ib] = int(i)
            else:
                nxt[j] = int(i)


@ti.kernel
def Calculate_a():
    for i in range(nP):
        if typ[i] == FLD:
            Acc_x, Acc_y, Acc_z = 0.0, 0.0, 0.0
            pix, piy, piz = px[i], py[i], pz[i]
            vix, viy, viz = vx[i], vy[i], vz[i]
            ix = int(((pix - MIN_X) * DBinv) + 1)
            iy = int(((piy - MIN_Y) * DBinv) + 1)
            iz = int(((piz - MIN_Z) * DBinv) + 1)
            for jz in range(iz - 1, iz + 2):
                for jy in range(iy - 1, iy + 2):
                    for jx in range(ix - 1, ix + 2):
                        jb = int(jz * nBxy + jy * nBx + jx)
                        j = int(bfst[jb])
                        if j == -1:
                            continue
                        while True:
                            v0 = px[j] - pix
                            v1 = py[j] - piy
                            v2 = pz[j] - piz
                            dist2 = v0 ** 2 + v1 ** 2 + v2 ** 2
                            if typ[j] != GST and dist2 < r2 and j != i:
                                dist = ti.sqrt(dist2)
                                wei = ((r / dist) - 1.0)
                                Acc_x += (vx[j] - vix) * wei
                                Acc_y += (vy[j] - viy) * wei
                                Acc_z += (vz[j] - viz) * wei
                            j = int(nxt[j])
                            if j == -1:
                                break
            ax[i] = Acc_x * A1 + G_X
            ay[i] = Acc_y * A1 + G_Y
            az[i] = Acc_z * A1 + G_Z


@ti.kernel
def Move():
    for i in range(nP):
        if typ[i] == FLD:
            vx[i] += ax[i] * DT
            vy[i] += ay[i] * DT
            vz[i] += az[i] * DT
            px[i] += vx[i] * DT
            py[i] += vy[i] * DT
            pz[i] += vz[i] * DT
            ax[i] = ay[i] = az[i] = 0
            Check_position(i)


@ti.kernel
def Move2():
    for i in range(nP):
        if typ[i] == FLD:
            vx[i] += ax[i] * DT
            vy[i] += ay[i] * DT
            vz[i] += az[i] * DT
            px[i] += ax[i] * DT * DT
            py[i] += ay[i] * DT * DT
            pz[i] += az[i] * DT * DT
            ax[i] = ay[i] = az[i] = 0
            Check_position(i)


@ti.kernel
def Collision():
    for i in range(nP):
        if typ[i] == FLD:
            pix, piy, piz = px[i], py[i], pz[i]
            vix, viy, viz = vx[i], vy[i], vz[i]
            vix2, viy2, viz2 = vx[i], vy[i], vz[i]
            ix = int(((pix - MIN_X) * DBinv) + 1)
            iy = int(((piy - MIN_Y) * DBinv) + 1)
            iz = int(((piz - MIN_Z) * DBinv) + 1)
            for jz in range(iz - 1, iz + 2):
                for jy in range(iy - 1, iy + 2):
                    for jx in range(ix - 1, ix + 2):
                        jb = int(jz * nBxy + jy * nBx + jx)
                        j = int(bfst[jb])
                        if j == -1:
                            continue
                        while True:
                            v0 = px[j] - pix
                            v1 = py[j] - piy
                            v2 = pz[j] - piz
                            dist2 = v0 ** 2 + v1 ** 2 + v2 ** 2
                            if typ[j] != GST and dist2 < rlim2 and j != i:
                                fDT = (vix - vx[j]) * v0 + (viy - vy[j]) * v1 + (viz - vz[j]) * v2
                                if fDT > 0.0:
                                    fDT *= COL * DNS / (DNS + DNS) / dist2
                                    vix2 -= v0 * fDT
                                    viy2 -= v1 * fDT
                                    viz2 -= v2 * fDT
                            j = int(nxt[j])
                            if j == -1:
                                break
            ax[i] = vix2
            ay[i] = viy2
            az[i] = viz2


@ti.kernel
def updav():
    for i in range(nP):
        vx[i] = ax[i]
        vy[i] = ay[i]
        vz[i] = az[i]


@ti.kernel
def Calculate_p():
    for i in range(nP):
        if typ[i] == FLD or typ[i] == WLL:
            pix, piy, piz = px[i], py[i], pz[i]
            ix = int(((pix - MIN_X) * DBinv) + 1)
            iy = int(((piy - MIN_Y) * DBinv) + 1)
            iz = int(((piz - MIN_Z) * DBinv) + 1)
            ni = 0.0
            for jz in range(iz - 1, iz + 2):
                for jy in range(iy - 1, iy + 2):
                    for jx in range(ix - 1, ix + 2):
                        jb = int(jz * nBxy + jy * nBx + jx)
                        j = int(bfst[jb])
                        if j == -1:
                            continue
                        while True:
                            v0 = px[j] - pix
                            v1 = py[j] - piy
                            v2 = pz[j] - piz
                            dist2 = v0 ** 2 + v1 ** 2 + v2 ** 2
                            if typ[j] != GST and dist2 < r2 and j != i:
                                dist = ti.sqrt(dist2)
                                wei = ((r / dist) - 1.0)
                                ni += wei
                            j = int(nxt[j])
                            if j == -1:
                                break
            pressure = (ni > n0) * (ni - n0) * A2 * DNS
            prs[i] = pressure


@ti.kernel
def Pressure_gradient():
    for i in range(nP):
        if typ[i] == FLD:
            Acc_x, Acc_y, Acc_z = 0.0, 0.0, 0.0
            pix, piy, piz = px[i], py[i], pz[i]
            pre_min = prs[i]
            ix = int(((pix - MIN_X) * DBinv) + 1)
            iy = int(((piy - MIN_Y) * DBinv) + 1)
            iz = int(((piz - MIN_Z) * DBinv) + 1)
            for jz in range(iz - 1, iz + 2):
                for jy in range(iy - 1, iy + 2):
                    for jx in range(ix - 1, ix + 2):
                        jb = int(jz * nBxy + jy * nBx + jx)
                        j = int(bfst[jb])
                        if j == -1:
                            continue
                        while True:
                            v0 = px[j] - pix
                            v1 = py[j] - piy
                            v2 = pz[j] - piz
                            dist2 = v0 ** 2 + v1 ** 2 + v2 ** 2
                            if typ[j] != GST and dist2 < r2 and j != i and pre_min > prs[j]:
                                pre_min = prs[j]
                            j = int(nxt[j])
                            if j == -1:
                                break
            for jz in range(iz - 1, iz + 2):
                for jy in range(iy - 1, iy + 2):
                    for jx in range(ix - 1, ix + 2):
                        jb = int(jz * nBxy + jy * nBx + jx)
                        j = int(bfst[jb])
                        if j == -1:
                            continue
                        while True:
                            v0 = px[j] - pix
                            v1 = py[j] - piy
                            v2 = pz[j] - piz
                            dist2 = v0 ** 2 + v1 ** 2 + v2 ** 2
                            if typ[j] != GST and dist2 < r2 and j != i:
                                dist = ti.sqrt(dist2)
                                wei = ((r / dist) - 1.0)
                                wei *= (prs[j] - pre_min) / dist2
                                Acc_x += v0 * wei
                                Acc_y += v1 * wei
                                Acc_z += v2 * wei
                            j = int(nxt[j])
                            if j == -1:
                                break
            ax[i] = Acc_x * invDNS * A3
            ay[i] = Acc_y * invDNS * A3
            az[i] = Acc_z * invDNS * A3


@ti.kernel
def Calculate_pav():
    for i in range(nP):
        pav[i] += prs[i]


T1 = time.perf_counter()
tim = 0.0
print("Loop start\n", "nP:", nP, ",nBxyz:", nBxyz)
Output()
for i in range(it + 1):
    if t % OPT_FQC == 0 and t != 0:
        print(t, "steps")
        Calculate_num()
        Output()
        if tim >= END_T:
            break
    Reset_bucket()
    Assigh_bucket()
    Calculate_a()
    Move()
    Collision()
    updav()
    Calculate_p()
    Pressure_gradient()
    Move2()
    Calculate_p()
    Calculate_pav()
    t += 1
    tim += DT
T2 = time.perf_counter()
print("Total loop time:", T2 - T1)
