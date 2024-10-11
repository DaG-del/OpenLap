import math
import helper_lib
import matplotlib.pyplot as mlt
import numpy
from scipy.interpolate import interp1d


f = open("motor.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(float(f[ef].split(",")[2]))


factor_power = line[0]
n_thermal = line[1]
fuel_LHV = line[2] # Energy per unit mass for battery also?
n_primary = line[3]
n_final = line[4]
n_gearbox = line[5] # replace with one efficiency factor?
ratio_primary = line[6]
ratio_final = line[7]
ratio_gearbox = line[8] # replace with one ratio?
max_power = line[9]*1000
max_torque = line[10]
max_rpm = line[11]
nog = 1


f = open("aero.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(float(f[ef].split(",")[2]))

Cl = line[0]
Cd = line[1]
factor_Cl = line[2]
factor_Cd = line[3]
da = line[4]/100
A = line[5]
rho = line[6]


f = open("mass_dimensions.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(float(f[ef].split(",")[2]))

M = line[0]
df = line[1]/100
L = line[2]/1000
rack = line[3]


f = open("tires.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(f[ef].split(",")[2])

factor_grip = float(line[0])
tyre_radius = float(line[1])/1000
Cr = float(line[2])
mu_x = float(line[3])
mu_x_M = float(line[4])
sens_x = float(line[5])
mu_y = float(line[6])
mu_y_M = float(line[7])
sens_y = float(line[8])
CF = float(line[9])
CR = float(line[10])

f = open("brake.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(f[ef].split(",")[2])

br_disc_d = float(line[0])/1000
br_pad_h = float(line[1])/1000
br_pad_mu = float(line[2])
br_nop = float(line[3])
br_pist_d = float(line[4])
br_mast_d = float(line[5])
br_ped_r = float(line[6])

br_pist_a = math.pow(br_nop*math.pi*(br_pist_d/1000), 2)/4
br_mast_a = math.pow(math.pi*(br_mast_d/1000), 2)/4
beta = tyre_radius/(br_disc_d/2 - br_pad_h/2)/br_pist_a/br_pad_mu/4
phi = br_mast_a/br_ped_r*2

a = (1 - df) * L
b = -df * L
C = [[2*CF, 2*(CF+CR)], [2*CF*a, 2*(CF*a + CR*b)]]

[torque, speed] = helper_lib.tsc(max_power, max_torque, max_rpm)
mlt.plot(speed, torque)
#mlt.show()

f = open("curve.csv", "w")
f = open("curve.csv", "a")
for i in range(len(speed)):
    f.write(str(speed[i]) + "," + str(torque[i]) + "\n")

power = []
wheel_speed_gear = []
vehicle_speed_gear = []
wheel_torque_gear = []
for s in range(len(speed)):
    power.append(torque[s] * speed[s] * 2 * math.pi/60)
    wheel_speed_gear.append(speed[s]/ratio_primary/ratio_gearbox/ratio_final)
    vehicle_speed_gear.append(wheel_speed_gear[s]*2*math.pi/60*tyre_radius)
    wheel_torque_gear.append(torque[s]*ratio_primary*n_primary*ratio_gearbox*n_gearbox*ratio_final*n_final)


v_min = min(vehicle_speed_gear)
v_max = max(vehicle_speed_gear)

dv = 0.5/3.6
vehicle_speed: list = numpy.linspace(v_min, v_max, math.floor((v_max-v_min)/dv)).tolist()

temp = []
for wtg in wheel_torque_gear:
    temp.append(wtg/tyre_radius)
fx_engine = []
gear = []
for i in range(len(vehicle_speed)):
    fx_engine.append(numpy.interp(vehicle_speed[i], vehicle_speed_gear, temp))
    gear.append(1)

vehicle_speed = [0] + vehicle_speed
gear.append(1)
temp = [fx_engine[0]]
temp.extend(fx_engine)
fx_engine = temp

engine_speed = []
for vs in vehicle_speed:
    engine_speed.append(ratio_final * ratio_gearbox * ratio_primary * vs / tyre_radius * 60 / 2 / math.pi)

wheel_torque = []
for fe in fx_engine:
    wheel_torque.append(fe*tyre_radius)

engine_torque = []
for wt in wheel_torque:
    engine_torque.append(wt / ratio_final / ratio_gearbox / ratio_primary / n_primary / n_gearbox / n_final)

# engine_power = engine_torque.*engine_speed*2*pi/60;
engine_power = []
for x in range(len(engine_torque)):
    engine_power.append(engine_torque[x] * engine_speed[x] * 2 * math.pi / 60)

g = 9.81
factor_drive = (1-df)
factor_aero = (1-da)
driven_wheels = 2

fz_mass = -M*g
fz_aero = []
fz_total = []
for vs in vehicle_speed:
    fz_aero.append(1/2*rho*factor_Cl*Cl*A*vs*vs)
    fz_total.append(1/2*rho*factor_Cl*Cl*A*vs*vs + fz_mass)

fz_tyre = []
for fa in fz_aero:
    fz_tyre.append((factor_drive*fz_mass+factor_aero*fa)/driven_wheels)

fx_aero = []
for vs in vehicle_speed:
    fx_aero.append(1/2*rho*factor_Cd*Cd*A*vs*vs)

fx_roll = []
for ft in fz_total:
    fx_roll.append(Cr * abs(ft))

fx_tyre = []
for ft in fz_tyre:
    fx_tyre.append(driven_wheels*(mu_x+sens_x*(mu_x_M*g-abs(ft)))*abs(ft))

bank = 0
incl = 0

dmy = factor_grip*sens_y
muy = factor_grip*mu_y
Ny = mu_y_M*g

dmx = factor_grip*sens_x
mux = factor_grip*mu_x
Nx = mu_x_M*g

Wz = M*g
Wy = 0
Wx = 0

dv = 2
v = helper_lib.vector(0, dv, v_max)
if v[-1] != v_max:
    v = v + [v_max]

N = 45
GGV = numpy.zeros((len(v), 2*N-1, 3))
fx_engine_time_factor_power = [fx_engine[i] * factor_power for i in range(len(fx_engine))]

for i in range(len(v)):
    Aero_Df = 0.5 * rho * factor_Cl * Cl * A * v[i]**2
    Aero_Dr = 0.5 * rho * factor_Cd * Cd * A * v[i]**2
    Roll_Dr = Cr * abs(-Aero_Df + Wz)
    Wd = (factor_drive * Wz + (-factor_aero * Aero_Df)) / driven_wheels
    ax_drag = (Aero_Dr + Roll_Dr + Wx) / M
    ay_max = (1 / M) * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
    ax_tyre_max_acc = (1 / M) * (mux + dmx * (Nx - Wd)) * Wd * driven_wheels
    ax_tyre_max_dec = -(1 / M) * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)
    interp_fx_engine = interp1d(vehicle_speed, fx_engine_time_factor_power, fill_value="extrapolate")
    ax_power_limit = (1 / M) * interp_fx_engine(v[i])
    ax_power_limit = ax_power_limit * numpy.ones(N)
    ay = ay_max * numpy.cos(numpy.radians(numpy.linspace(0, 180, N)))
    ax_tyre_acc = ax_tyre_max_acc * numpy.sqrt(1 - (ay / ay_max) ** 2)
    ax_acc = numpy.minimum(ax_tyre_acc, ax_power_limit) + ax_drag
    ax_dec = ax_tyre_max_dec * numpy.sqrt(1 - (ay / ay_max) ** 2) + ax_drag  # Friction ellipse
    GGV[i, :, 0] = numpy.concatenate([ax_acc, ax_dec[1:]])
    GGV[i, :, 1] = numpy.concatenate([ay, numpy.flipud(ay[1:])])
    GGV[i, :, 2] = v[i] * numpy.ones(2 * N - 1)
print("debug")