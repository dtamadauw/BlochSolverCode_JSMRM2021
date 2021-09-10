from psdk import *
import numpy as np

gamma = 42.57747892 # [MHz/T]
TR = 15.0e+3 # [us]
TE = 6.0e+3 # [us]
NR = 256 # Number of readout points
NPE1 = 256 # Number of 1st phase encoding
fov = [256.0, 256.0, 256.0] # [mm]
dwell_time = 10.0 # [us]
slice_width = 5.0 # [mm]
gx_value = 1e+6 / (dwell_time * gamma * fov[0]) # [mT/m]
gy_value = 2e+6 / (dwell_time * gamma * fov[1]) * NPE1 / NR # [mT/m]
gz_value = 1.25 / (slice_width * 1.0e-3) / gamma # [mT/m]
gx_rt = 300.0 # [us] gx rise time
gy_rt = 300.0 # [us] gy rise time
gz_rt = 300.0 # [us] gz rise time
ex_pulse_width = 3200.0 # [us]
ex_pulse_flip_angle = 30.0 # [degree]

def sinc_with_hamming(flip_angle, pulse_width, points, *, min = -2.0 * np.pi, max = 2.0 * np.pi):
    x0 = np.arange(min, max, (max - min) / points)
    x1 = x0 + (max - min) / points
    y = (np.sinc(x0 / np.pi) + np.sinc(x1 / np.pi)) * 0.5 * np.hamming(points)
    return flip_angle * y * points / (y.sum() * pulse_width * 360.0e-6 * gamma)

with Sequence('2D GradientEcho'):

    with Block('Excitation', ex_pulse_width + 2.0*gz_rt):
        GZ(0.0, gz_value, gz_rt)
        RF(gz_rt, sinc_with_hamming(ex_pulse_flip_angle, ex_pulse_width, 160), ex_pulse_width / 160)
        GZ(ex_pulse_width + gz_rt, 0.0, gz_rt)

    with Block('PhaseEncoding', NR // 2 * dwell_time + gx_rt * 2.5):
        GX(0.0, -gx_value, gx_rt)
        GY(0.0, ([gy_value * (i - NPE1 // 2) / NPE1 for i in range(NPE1)], ['PE1']), gy_rt)
        GY(NR // 2 * dwell_time, 0.0, gy_rt)
        GX(NR // 2 * dwell_time + gx_rt * 0.5, gx_value, gx_rt * 2.0)
        GZ(0.0, -gz_value, gz_rt * 0.6)
        GZ(ex_pulse_width * 0.5 + 150.0, 0.0, gz_rt * 0.6) 

    with Block('Readout', NR * dwell_time):
        AD(0.0, NR, dwell_time)

    with Block('Rewinding', NR // 2 * dwell_time + gx_rt*2):
        GY(0.0, ([gy_value * (NPE1 // 2 - i) / NPE1 for i in range(NPE1)], ['PE1']), gy_rt)
        GX(0.0, 0.0, gx_rt)
        GX(gx_rt, gx_value, gx_rt)
        GY(NR // 2 * dwell_time, 0.0, gy_rt)
        GX(gx_rt+NR // 2 * dwell_time -0.5*gx_rt, 0, gx_rt)

    with Main():
        with Loop('PE1', NPE1):
            BlockRef('Excitation')
            WaitUntil(TE + ex_pulse_width * 0.5 + gz_rt - NR // 2 * 2 * dwell_time - gx_rt * 2.5)
            BlockRef('PhaseEncoding')
            BlockRef('Readout')
            BlockRef('Rewinding')
            WaitUntil(TR)