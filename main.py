import math


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


f = open("torque_speed_curve.csv", "r")
f = f.readlines()
speed = []
torque = []

for ef in range(len(f)):
    speed.append(float(f[ef].split(",")[0]))
    torque.append(float(f[ef].split(",")[1]))

power = [torque[i] * speed[i] * 2 * math.pi / 60 for i in range(len(speed))]

wheel_speed_gear = [speed[i] /ratio_primary/ratio_gearbox/ratio_final for i in range(len(speed))]
vehicle_speed_gear = []
vehicle_torque_gear = []