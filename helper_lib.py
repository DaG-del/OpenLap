import math

def tsc(max_power, max_torque, max_rpm):
    rpms = list(range(100, int(max_rpm), 100))
    torque_limit = [max_torque] * len(rpms)
    power_limit = [max_power/(rpms[i] * 2 * math.pi / 60) for i in range(len(rpms))]
    torque = [min(torque_limit[i], power_limit[i]) for i in range(len(rpms))]
    return [torque, rpms]


def vector(min, dv, v_max):
    n = []
    curr = min
    while curr <= v_max:
        n.append(curr)
        curr += dv

    return n

def compare_floats(a, b, error):
    a_rounded = []
    for ay in a:
        a_rounded.append(round(ay, 2))
    b_rounded = []
    for by in b:
        b_rounded.append(round(by, 2))

    if len(b_rounded) != len(a_rounded):
        print("Unequal length!")
        return

    errors = {}
    for x in range(len(a_rounded)):
        if a_rounded[x] - b_rounded[x] > error:
            errors[x] = [a_rounded[x], b_rounded[x]]

    for e in errors:
        print(e, end = " : ")
        print(errors[e])


def read_csv(filename):
    f = open(filename, "r")
    ret = []
    ip = f.readlines()
    for i in ip:
        ret.append(float(i))

    return ret


def rotz(angle, dsv, p):
    angle *= math.pi/180
    x = dsv*math.cos(angle) + p[0]
    y = dsv*math.sin(angle) + p[1]

    return [x, y]


def absol(r: list):
    r_abs = []
    for aar in r:
        r_abs.append(abs(aar))
    return r_abs


def find(r_apex_indices, r):
    ret = []
    for aar in r_apex_indices:
        ret.append(r[aar])
    return ret


def calc_d(x, total_length, Y):
    '''
    % linear correction vectors
    DX = x/L*(X(1)-X(end)) ;
    DY = x/L*(Y(1)-Y(end)) ;
    DZ = x/L*(Z(1)-Z(end)) ;
    % adding correction
    X = X+DX ;
    Y = Y+DY ;
    Z = Z+DZ ;
    '''
    ret = []
    DX = []
    for el in range(len(x)):
        DX.append(x[el] / total_length * (Y[0] - Y[-1]))

    for el in range(len(DX)):
        ret.append(Y[el] + DX[el])

    return ret


def writecsv(r_apex, distance_step_vector, n, r, x, X, Y, Z):
    f = open("OpenTRACK_output.csv", "a")
    d = {}
    d["r_apex"] = r_apex
    d["dx"] = distance_step_vector
    d["n"] = n
    d["r"] = r
    d["x"] = x
    d["X"] = X
    d["Y"] = Y
    d["Z"] = Z
    f.write(str(d))