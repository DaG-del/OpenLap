import math


f = open("tires.csv", "r")
f = f.readlines()
line = []

for ef in range(len(f)):
    line.append(f[ef].split(",")[2])

M = float(line[0].split(",")[2])
df = float(line[1].split(",")[2])
L = float(line[2].split(",")[2])
rack = float(line[3].split(",")[2])


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