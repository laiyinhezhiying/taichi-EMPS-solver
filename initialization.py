from pyevtk.hl import *
import numpy as np

MAXX = 1.0
MINX = 0.0
MAXY = 0.2
MINY = 0.0
MAXZ = 0.6
MINZ = 0.0
GHOST = -1
FLUID = 0
WALL = 1

WAVE_HEIGHT = 0.5
WAVE_WIDTH = 0.25

PARTICLE_DISTANCE = 0.02
nx = int(((MAXX - MINX) / PARTICLE_DISTANCE) + 6)
ny = int(((MAXY - MINY) / PARTICLE_DISTANCE) + 6)
nz = int(((MAXZ - MINZ) / PARTICLE_DISTANCE) + 6)
nxy = int(nx * ny)
nxyz = int(nxy * nz)
Type = -1 * np.ones((nxyz))
Pos = np.zeros((3 * nxyz))
Vel = np.zeros((3 * nxyz))
Prs = np.zeros((nxyz))
Num = 0
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            ip = int(iz * nxy + iy * nx + ix)
            Type[ip] = GHOST
            Pos[ip] = MINX + PARTICLE_DISTANCE * (ix - 3)
            Pos[ip + nxyz] = MINY + PARTICLE_DISTANCE * (iy - 3)
            Pos[ip + 2 * nxyz] = MINZ + PARTICLE_DISTANCE * (iz - 3)
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            ip = int(iz * nxy + iy * nx + ix)
            if (ix < 3 or ix >= nx - 3 or iy < 3 or iy >= ny - 3 or iz < 3):
                Type[ip] = WALL
                Num += 1
            elif (Pos[ip + 2 * nxyz] <= WAVE_HEIGHT and Pos[ip] <= WAVE_WIDTH):
                Type[ip] = FLUID
                Num += 1
            # elif (Pos[ip + 2 * nxyz] <= WAVE_HEIGHT and Pos[ip] >= MAXX - WAVE_WIDTH):
            #     Type[ip] = FLUID
            #     Num += 1

Num = str(Num)

def Output(Pos2, Vel2, Prs2, Type2):
    px = Pos2[0:k]
    py = Pos2[k:2 * k]
    pz = Pos2[2 * k:3 * k]
    vx = Vel2[0:k]
    vy = Vel2[k:2 * k]
    vz = Vel2[2 * k:3 * k]
    Velocity = np.sqrt(vx * vx + vy * vy + vz * vz)
    point_data = {"Pressure": Prs2, "Type": Type2, "Velocity": Velocity}
    pointsToVTK("int", px, py, pz, data=point_data)


k = 0
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            ip = int(iz * nxy + iy * nx + ix)
            if (Type[ip] == GHOST): continue
            k += 1
n = 0
Type2 = -1 * np.ones((k))
Pos2 = np.zeros((3 * k))
Vel2 = np.zeros((3 * k))
Prs2 = np.zeros((k))
Pav2 = np.zeros((k))
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            ip = int(iz * nxy + iy * nx + ix)
            if (Type[ip] == GHOST): continue
            Type2[n] = Type[ip]
            Pos2[n] = Pos[ip]
            Pos2[n + k] = Pos[ip + nxyz]
            Pos2[n + 2 * k] = Pos[ip + 2 * nxyz]
            Vel2[n] = Vel[ip]
            Vel2[n + 2 * k] = Vel[ip + 2 * nxyz]
            Prs2[n] = Prs[ip]
            n += 1
b = np.arange(0, k, 1)
nP = str(k)
np.savetxt("int3dpy.txt",
           np.column_stack((b, Type2, Pos2[0:k], Pos2[k:2 * k], Pos2[2 * k:3 * k], Vel2[0:k], Vel2[k:2 * k],
                            Vel2[2 * k:3 * k], Prs2)),
           fmt='%d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f ', header=nP, comments='')
Output(Pos2, Vel2, Prs2, Type2)
print(k)
