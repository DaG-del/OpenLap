import math

def tsc(max_power, max_torque, max_rpm):
    rpms = list(range(0, max_rpm, 100))
    torque_limit = [max_torque] * len(rpms)
    power_limit = [max_power/(rpms[i] * 2 * math.pi / 60) for i in range(len(rpms))]
    torque = [min(torque_limit[i], power_limit[i]) for i in range(len(rpms))]
    return [torque, rpms]