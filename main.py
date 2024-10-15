import logging
import math

import scipy
import scipy.interpolate as interp
import helper_lib
import matplotlib.pyplot as mlt
import numpy
from scipy.interpolate import interp1d


def center_points(l):
    sum = 0
    cumulative = []
    for el in l:
        sum += el
        cumulative.append(sum - el / 2)

    return cumulative


def end_points(l):
    sum = 0
    cumulative = []
    for el in l:
        sum += el
        cumulative.append(sum)

    return cumulative


def nos(R):
    n = 0
    for r in R:
        if r == float("inf"):
            n += 1

    return n


def fine(total_length):
    tl = total_length
    x_fine = []
    sum = 0
    while True:
        x_fine.append(sum)
        sum += _mesh_size
        if x_fine[-1] > math.floor(tl):
            x_fine.pop()
            x_fine.append(total_length)
            return x_fine


def min_custom(dh):
    dh_abs = []
    for d in dh:
        dh_abs.append(abs(d))
    return min(dh_abs)


_KAPPA = 10
_mesh_size = 0.1
f = open("track.txt", "r")
track = []
ip = f.readlines()

for l in ip:
    l = l[:-1]
    tup = tuple(l.split(","))
    track.append(tup)

track = tuple(track)

R = [float(t[2]) for t in track]
for ar in range(len(R)):
    if R[ar] == 0:
        R[ar] = float("inf")

l = [float(t[1]) for t in track]
type = [-int(t[0]) for t in track]
total_length = sum(l)

track = [(type[i], l[i], R[i]) for i in range(len(track))]

angle_seg = [(l[i] / R[i]) * 180 / math.pi for i in range(len(track))]

R_injected = R
l_injected = l
type_injected = type

j = 0
for i in range(len(l)):
    if angle_seg[i] > _KAPPA:
        distance_of_injected_point_from_corner_entry_and_exit = min(l_injected[j] / 3, (_KAPPA * math.pi / 180) * R[i])
        temp = l_injected[0:j]
        temp.append(distance_of_injected_point_from_corner_entry_and_exit)
        temp.append(l_injected[j] - 2 * distance_of_injected_point_from_corner_entry_and_exit)
        temp.append(distance_of_injected_point_from_corner_entry_and_exit)
        temp.extend(l_injected[j + 1:])
        l_injected = temp

        temp = R_injected[0:j]
        temp.extend([R_injected[j]] * 3)
        temp.extend(R_injected[j + 1:])

        R_injected = temp

        temp = type_injected[0:j]
        temp.extend([type_injected[j]] * 3)
        temp.extend(type_injected[j + 1:])

        type_injected = temp

        del temp

        j += 3
    else:
        j += 1

R = R_injected
l = l_injected
type = type_injected

'''
for i=1:length(l)-1
    j = 1 ;
    while true
        if type(i+j)==0 && type(i)==0 && l(i)~=-1
            l(i) = l(i)+l(i+j) ;
            l(i+j) = -1 ;
        else
            break
        end
        j = j+1 ;
    end
end
R(l==-1) = [] ;
type(l==-1) = [] ;
l(l==-1) = [] ;
'''

for i in range(len(l) - 1):
    j = 1
    while True:
        if type[i + j] == 0 and type[i] == 0 and l[i] != -1:
            l[i] = l[i + j] + l[i]
            l[i + j] = -1
        else:
            break
        j += 1

temp = len(l) - 1
el = 0
while el < temp:
    if l[el] == -1:
        l.pop(el)
        R.pop(el)
        type.pop(el)
        temp -= 1
    el += 1

segment_end_point = end_points(l)

segment_center_point = center_points(l)

no_of_straights = nos(R)
x_coarse = [0] * (len(segment_end_point) + no_of_straights)
r: list = [0] * len(x_coarse)

j = 0
for i in range(len(segment_center_point)):
    if R[i] == float("inf"):
        x_coarse[j] = segment_end_point[i] - l[i]
        x_coarse[j + 1] = segment_end_point[i]
        j += 2
    else:
        x_coarse[j] = segment_center_point[i]
        r[j] = type[i] / R[i]
        j += 1

x = fine(total_length)

distance_step_vector = numpy.diff(x)
distance_step_vector = distance_step_vector.tolist()
distance_step_vector.append(distance_step_vector[-1])
number_of_mesh_points = len(x)

r = interp.pchip_interpolate(x_coarse, r, x)
n = len(x)

X = [0] * n
Y = [0] * n

angle_seg = []
for d in range(len(distance_step_vector)):
    angle_seg.append((distance_step_vector[d] * r[d]) * (180 / math.pi))

angle_heading = numpy.cumsum(angle_seg)
angle_heading = list(angle_heading.tolist())

dh = []

temp = angle_heading[-1] % (360 * (abs(angle_heading[-1]) / angle_heading[-1]))
dh.append(temp)
temp = angle_heading[-1] - (360 * (abs(angle_heading[-1]) / angle_heading[-1]))
dh.append(temp)

dh = min(dh)

for ah in range(len(angle_heading)):
    angle_heading[ah] -= x[ah] / total_length * dh

angle_seg.clear()
angle_seg.append(angle_heading[0])
diff_as = numpy.diff(angle_heading)
diff_as = list(diff_as.tolist())
angle_seg.extend(diff_as)

temp = angle_heading[0]

for ah in range(len(angle_heading)):
    angle_heading[ah] -= temp

for i in range(1, len(x)):
    p = [X[i - 1], Y[i - 1]]

    xyz = helper_lib.rotz(angle_heading[i - 1], distance_step_vector[i - 1], p)

    X[i] = xyz[0]
    Y[i] = xyz[1]

r_abs = helper_lib.absol(r)
r_apex_indices = scipy.signal.find_peaks(r_abs)
r_apex_indices = r_apex_indices[0]
r_apex_indices = r_apex_indices.tolist()

r_apex = helper_lib.find(r_apex_indices, r)

X = helper_lib.calc_d(x, total_length, X)
Y = helper_lib.calc_d(x, total_length, Y)
Z = [0] * len(X)

########################################################################################################
########################################################################################################
########################################################################################################


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

######################################################################################################
######################################################################################################
######################################################################################################

def vehicle_model_lat(i, r):

    r = r[p]
    # Speed solution
    if r == 0:  # Straight road (limited by engine speed limit or drag)
        v = v_max
        tps = 1  # Full throttle
        bps = 0  # No brake
    else:  # Cornering (may be limited by engine, drag, or cornering ability)
        # Initial speed solution
        D = -0.5 * rho * factor_Cl * Cl * A  # Downforce coefficient
        dmy = factor_grip * sens_y
        muy = factor_grip * mu_y
        Ny = mu_y_M * g

        dmx = factor_grip * sens_x
        mux = factor_grip * mu_x
        Nx = mu_x_M * g

        # Polynomial coefficients for quadratic equation a*x^2 + b*x + c = 0
        a = -numpy.sign(r) * dmy / 4 * D ** 2
        b = numpy.sign(r) * (muy * D + (dmy / 4) * (Ny * 4) * D - 2 * (dmy / 4) * Wz * D) - M * r
        c = numpy.sign(r) * (muy * Wz + (dmy / 4) * (Ny * 4) * Wz - (dmy / 4) * Wz ** 2) + Wy

        # Solving the quadratic equation
        if a == 0:
            v = numpy.sqrt(-c / b)
        elif b ** 2 - 4 * a * c >= 0:
            root1 = (-b + numpy.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            root2 = (-b - numpy.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
            if root1 >= 0:
                v = numpy.sqrt(root1)
            elif root2 >= 0:
                v = numpy.sqrt(root2)
            else:
                raise ValueError(f'No real roots at point index: {p}')
        else:
            raise ValueError(f'Discriminant < 0 at point index: {p}')

        # Checking for engine speed limit
        v = min(v, v_max)

        # Adjusting speed for drag force compensation
        adjust_speed = True
        while adjust_speed:
            # Aero forces
            Aero_Df = 0.5 * rho * factor_Cl* Cl * A * v ** 2
            Aero_Dr = 0.5 * rho * factor_Cd * Cd * A * v ** 2

            # Rolling resistance
            Roll_Dr = Cr * (-Aero_Df + Wz)

            # Normal load on driven wheels
            Wd = (factor_drive * Wz + (-factor_aero * Aero_Df)) / driven_wheels

            # Drag acceleration
            ax_drag = (Aero_Dr + Roll_Dr + Wx) / M

            # Maximum lateral acceleration available from tyres
            ay_max = numpy.sign(r) / M * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)

            # Needed lateral acceleration to make the turn
            ay_needed = v ** 2 / r + g * numpy.sin(numpy.radians(bank))

            # Calculating driver inputs
            if ax_drag <= 0:  # Need throttle to compensate for drag
                ax_tyre_max_acc = (1 / M) * (mux + dmx * (Nx - Wd)) * Wd * driven_wheels
                temp = [fx_engine[i] * factor_power for i in range(len(fx_engine))]
                ax_power_limit = (1 / M) * interp1d(vehicle_speed, temp)(v)

                # Available combined lateral acceleration when ax_net == 0
                ay = ay_max * numpy.sqrt(1 - (ax_drag / ax_tyre_max_acc) ** 2)

                # Available combined longitudinal acceleration at ay_needed
                ax_acc = ax_tyre_max_acc * numpy.sqrt(1 - (ay_needed / ay_max) ** 2)

                # Throttle input (tps)
                scale = min([-ax_drag, ax_acc]) / ax_power_limit
                tps = max(min(1, scale), 0)
                bps = 0  # No brake pressure
            else:  # Need brake to compensate for drag
                ax_tyre_max_dec = -(1 / M) * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df)

                # Available combined lateral acceleration when ax_net == 0
                ay = ay_max * numpy.sqrt(1 - (ax_drag / ax_tyre_max_dec) ** 2)

                # Available combined longitudinal deceleration at ay_needed
                ax_dec = ax_tyre_max_dec * numpy.sqrt(1 - (ay_needed / ay_max) ** 2)

                # Brake input (bps)
                fx_tyre = max([ax_drag, -ax_dec]) * M
                bps = max([fx_tyre, 0]) * beta
                tps = 0  # No throttle

            # Check if tyres can produce the available combined lateral acceleration
            if ay / ay_needed < 1:  # Not enough grip
                v = numpy.sqrt((ay - g * numpy.sin(numpy.radians(bank))) / r) - 1E-3  # Adjust speed for convergence
            else:  # Enough grip
                adjust_speed = False

    return v, tps, bps

def other_points(i, i_max):
    i_rest = numpy.arange(1, i_max + 1)
    i_rest = numpy.delete(i_rest, i - 1)  # Remove the i-th element (adjusting for 0-based index)
    return i_rest

def next_point(j, j_max, mode, tr_config):
    if mode == 1:  # Acceleration
        if tr_config == 'Closed':
            if j == j_max - 1:
                j = j_max
                j_next = 1
            elif j == j_max:
                j = 1
                j_next = j + 1
            else:
                j = j + 1
                j_next = j + 1
        elif tr_config == 'Open':
            j = j + 1
            j_next = j + 1

    elif mode == -1:  # Deceleration
        if tr_config == 'Closed':
            if j == 2:
                j = 1
                j_next = j_max
            elif j == 1:
                j = j_max
                j_next = j - 1
            else:
                j = j - 1
                j_next = j - 1
        elif tr_config == 'Open':
            j = j - 1
            j_next = j - 1

    return j_next, j

def vehicle_model_comb(veh, tr, v, v_max_next, j, mode, dx, r):

    overshoot = False  # assuming no overshoot initially

    # Getting track data
    dx = dx[j]
    r = r[j]
    factor_grip = factor_grip[j] * veh['factor_grip']
    g = 9.81

    # Getting vehicle data
    if mode == 1:
        factor_drive = veh['factor_drive']
        factor_aero = veh['factor_aero']
        driven_wheels = veh['driven_wheels']
    else:
        factor_drive = 1
        factor_aero = 1
        driven_wheels = 4

    # External forces
    M = veh['M']
    Wz = M * g * numpy.cos(numpy.radians(bank)) * numpy.cos(numpy.radians(incl))  # normal load on all wheels
    Wy = -M * g * numpy.sin(numpy.radians(bank))  # induced weight from banking
    Wx = M * g * numpy.sin(numpy.radians(incl))  # induced weight from inclination

    Aero_Df = 0.5 * veh['rho'] * veh['factor_Cl'] * veh['Cl'] * veh['A'] * v ** 2  # downforce
    Aero_Dr = 0.5 * veh['rho'] * veh['factor_Cd'] * veh['Cd'] * veh['A'] * v ** 2  # drag force
    Roll_Dr = veh['Cr'] * (-Aero_Df + Wz)  # rolling resistance
    Wd = (factor_drive * Wz + (-factor_aero * Aero_Df)) / driven_wheels  # normal load on driven wheels

    # Overshoot acceleration
    ax_max = mode * (v_max_next ** 2 - v ** 2) / (
                2 * dx)  # maximum allowed longitudinal acceleration to avoid overshoot
    ax_drag = (Aero_Dr + Roll_Dr + Wx) / M  # drag acceleration
    ax_needed = ax_max - ax_drag  # required acceleration to avoid overshoot

    # Current lateral acceleration
    ay = v ** 2 / r + g * numpy.sin(numpy.radians(bank))

    # Tyre forces
    dmy = factor_grip * veh['sens_y']
    muy = factor_grip * veh['mu_y']
    Ny = veh['mu_y_M'] * g

    dmx = factor_grip * veh['sens_x']
    mux = factor_grip * veh['mu_x']
    Nx = veh['mu_x_M'] * g

    # Friction ellipse multiplier
    if numpy.sign(ay) != 0:  # In corner or compensating for banking
        ay_max = 1 / M * (numpy.sign(ay) * (muy + dmy * (Ny - (Wz - Aero_Df) / 4)) * (Wz - Aero_Df) + Wy)
        if abs(ay / ay_max) > 1:  # Checking if vehicle overshot
            ellipse_multi = 0  # Overshot condition
        else:
            ellipse_multi = numpy.sqrt(1 - (ay / ay_max) ** 2)  # Friction ellipse
    else:  # In straight or no compensation needed
        ellipse_multi = 1

    # Calculating driver inumpyuts
    if ax_needed >= 0:  # Throttle required
        ax_tyre_max = 1 / M * (mux + dmx * (Nx - Wd)) * Wd * driven_wheels  # Max pure longitudinal acceleration
        ax_tyre = ax_tyre_max * ellipse_multi  # Max combined longitudinal acceleration

        # Engine power limit
        fx_engine_interp = interp1d(veh['vehicle_speed'], veh['factor_power'] * veh['fx_engine'],
                                    fill_value="extrapolate")
        ax_power_limit = 1 / M * fx_engine_interp(v)

        # Throttle position
        scale = min([ax_tyre, ax_needed] / ax_power_limit)
        tps = max(min(1, scale), 0)  # Making sure it's between 0 and 1
        bps = 0  # No braking
        ax_com = tps * ax_power_limit  # Final longitudinal acceleration command
    else:  # Braking required
        ax_tyre_max = -1 / M * (mux + dmx * (Nx - (Wz - Aero_Df) / 4)) * (
                    Wz - Aero_Df)  # Max pure longitudinal deceleration
        ax_tyre = ax_tyre_max * ellipse_multi  # Max combined longitudinal deceleration

        fx_tyre = min(-[ax_tyre, ax_needed]) * M  # Tyre braking force
        bps = max(fx_tyre, 0) * veh['beta']  # Brake pressure
        tps = 0  # No throttle
        ax_com = -min(-[ax_tyre, ax_needed])  # Final longitudinal deceleration command

    # Final results
    ax = ax_com + ax_drag  # Total longitudinal acceleration
    v_next = numpy.sqrt(v ** 2 + 2 * mode * ax * dx)  # Next speed value

    # Correcting throttle for full throttle when at v_max on straights
    if tps > 0 and v / veh['v_max'] >= 0.999:
        tps = 1

    # Checking for overshoot
    if v_next / v_max_next > 1:
        overshoot = True
        v_next = numpy.inf  # Reset speed for overshoot
        ax = 0
        ay = 0
        tps = -1  # Special value for throttle
        bps = -1  # Special value for brake
        return v_next, ax, ay, tps, bps, overshoot

    return v_next, ax, ay, tps, bps, overshoot

def flag_update(flag, j, k, prg_size, logid, prg_pos):
    """
    Updates the flag matrix and checks if the progress bar needs updating.
    """
    # Current flag state
    p = numpy.sum(flag) / (flag.shape[0] * flag.shape[1])  # Calculate percentage of flags set to True
    n_old = int(p * prg_size)  # Old number of progress steps

    # New flag state
    flag[j, k] = True  # Update the flag at position (j, k)
    p = numpy.sum(flag) / (flag.shape[0] * flag.shape[1])  # Recalculate percentage after update
    n = int(p * prg_size)  # New number of progress steps

    return flag

def simulate(veh, tr, logid):
    v_max = numpy.zeros(tr.n, dtype=numpy.float32)
    bps_v_max = numpy.zeros(tr.n, dtype=numpy.float32)
    tps_v_max = numpy.zeros(tr.n, dtype=numpy.float32)
    for i in range(tr.n):
        v_max[i], tps_v_max[i], bps_v_max[i] = vehicle_model_lat(i, r)

    v_apex, apex = scipy.findpeaks(-v_max)
    v_apex = -v_apex
    if tr.info['config'] == 'Open':
        if apex[0] != 0:
            apex = numpy.insert(apex, 0, 0)
            v_apex = numpy.insert(v_apex, 0, 0)
        else:
            v_apex[0] = 0

    if len(apex) == 0:
        v_apex, apex = numpy.min(v_max), numpy.argmin(v_max)

    apex_table = numpy.column_stack((v_apex, apex))
    apex_table = apex_table[numpy.argsort(apex_table[:, 0])]
    v_apex = apex_table[:, 0]
    apex = apex_table[:, 1].astype(int)

    tps_apex = tps_v_max[apex]
    bps_apex = bps_v_max[apex]

    N = len(apex)
    flag = numpy.zeros((tr.n, 2), dtype=bool)
    v = numpy.full((tr.n, N, 2), numpy.inf, dtype=numpy.float32)
    ax = numpy.zeros((tr.n, N, 2), dtype=numpy.float32)
    ay = numpy.zeros((tr.n, N, 2), dtype=numpy.float32)
    tps = numpy.zeros((tr.n, N, 2), dtype=numpy.float32)
    bps = numpy.zeros((tr.n, N, 2), dtype=numpy.float32)

    prg_size = 30
    prg_pos = 0

    for i in range(N):
        for k in range(1, 3):
            mode = 1 if k == 1 else -1
            k_rest = 2 if k == 1 else 1
            if not (tr.info['config'] == 'Open' and mode == -1 and i == 0):
                i_rest = other_points(i, N)
                if i_rest is None:
                    i_rest = i

                j = apex[i]
                v[j, i, k - 1] = v_apex[i]
                ay[j, i, k - 1] = v_apex[i] ** 2 * tr.r[j]
                tps[j, :, 0] = tps_apex[i]
                bps[j, :, 0] = bps_apex[i]
                tps[j, :, 1] = tps_apex[i]
                bps[j, :, 1] = bps_apex[i]
                flag[j, k - 1] = True

                _, j_next = next_point(j, tr.n, mode, tr.info['config'])
                if not (tr.info['config'] == 'Open' and mode == 1 and i == 0):
                    v[j_next, i, k - 1] = v[j, i, k - 1]
                    j_next, j = next_point(j, tr.n, mode, tr.info['config'])

                while True:
                    logging.info(f"{i}\t{j}\t{k}\t{tr.x[j]}\t{v[j, i, k - 1]}\t{v_max[j]}")
                    v[j_next, i, k - 1], ax[j, i, k - 1], ay[j, i, k - 1], tps[j, i, k - 1], bps[
                        j, i, k - 1], overshoot = vehicle_model_comb(veh, tr, v[j, i, k - 1], v_max[j_next], j, mode)

                    if overshoot:
                        break

                    if flag[j, k - 1] or flag[j, k_rest - 1]:
                        if max(v[j_next, i, k - 1] >= v[j_next, i_rest, k - 1]) or max(
                                v[j_next, i, k - 1] > v[j_next, i_rest, k_rest - 1]):
                            break

                    flag = flag_update(flag, j, k, prg_size, logid, prg_pos)
                    j_next, j = next_point(j, tr.n, mode, tr.info['config'])

                    if tr.info['config'] == 'Closed':
                        if j == apex[i]:
                            break
                    elif tr.info['config'] == 'Open':
                        if j == tr.n or j == 0:
                            break

    V, AX, AY, TPS, BPS = numpy.zeros(tr.n), numpy.zeros(tr.n), numpy.zeros(tr.n), numpy.zeros(tr.n), numpy.zeros(tr.n)
    for i in range(tr.n):
        IDX = len(v[i, :, 0])
        v = v.tolist()
        V[i], idx = min([v[i, :, 0], v[i, :, 1]])
        if idx < IDX:
            AX[i] = ax[i, idx, 0]
            AY[i] = ay[i, idx, 0]
            TPS[i] = tps[i, idx, 0]
            BPS[i] = bps[i, idx, 0]
        else:
            AX[i] = ax[i, idx - IDX, 1]
            AY[i] = ay[i, idx - IDX, 1]
            TPS[i] = tps[i, idx - IDX, 1]
            BPS[i] = bps[i, idx - IDX, 1]

    if tr.info['config'] == 'Open':
        time = numpy.cumsum([tr.dx[1] / V[1]] + list(tr.dx[1:] / V[1:]))
    else:
        time = numpy.cumsum(tr.dx / V)

    sector_time = numpy.zeros(numpy.max(tr.sector))
    for i in range(1, numpy.max(tr.sector) + 1):
        sector_time[i - 1] = numpy.max(time[tr.sector == i]) - numpy.min(time[tr.sector == i])

    laptime = time[-1]

    A = numpy.sqrt(AX ** 2 + AY ** 2)
    Fz_mass = -M * g * numpy.cos(numpy.radians(tr.bank)) * numpy.cos(numpy.radians(tr.incl))
    Fz_aero = 0.5 * veh.rho * veh.factor_Cl * veh.Cl * A * V ** 2
    Fz_total = Fz_mass + Fz_aero
    Fx_aero = 0.5 * veh.rho * veh.factor_Cd * veh.Cd * A * V ** 2
    Fx_roll = veh.Cr * numpy.abs(Fz_total)

    yaw_rate = V * tr.r
    delta = numpy.zeros(tr.n)
    beta = numpy.zeros(tr.n)
    for i in range(tr.n):
        B = numpy.array([M * V[i] ** 2 * tr.r[i] + M * g * numpy.sin(numpy.radians(tr.bank[i])), 0])
        sol = numpy.linalg.solve(veh.C, B)
        delta[i] = sol[0] + numpy.degrees(numpy.arctan(veh.L * tr.r[i]))
        beta[i] = sol[1]

    steer = delta * veh.rack

    wheel_torque = TPS * interp.interp1d(veh.enginerpm_axis, veh.enginetorque_axis)(V / veh.Re)
    P_loss = TPS * (interp.interp1d(veh.enginerpm_axis, veh.enginefr_axis)(V / veh.Re) + Fx_aero + Fx_roll)
    Power = P_loss + (TPS * wheel_torque * V / veh.Re)

    return laptime, sector_time, V, AX, AY, yaw_rate, steer, beta, wheel_torque, Power
    