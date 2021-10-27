# ----- Calculate the longitudinal dynamic stability of rigid configuration -----

# -------- Import libraries ---------
import pandas as pd
import dyn_functions as dynfn
import numpy as np
import csv
import aero_functions as aerofn
from datetime import date

# -------- Generic run script -------
rho = 1.225
# taken from the bird body specimen tested in AvInertia - chord was scaled to match this bird
m = 1.0154  # mass (kg)
W = m * 9.81  # weight (N)
S_max = 0.244657  # wings and body reference area
c_max = 0.2861011  # wing root chord
coef_data = pd.read_csv('/Volumes/GoogleDrive/My Drive/DoctoralThesis/Chapter3_DynamicStability'
                        '/2021_10_27_coefficients_sw0_di20.csv')

a_test = [-2, 0, 2, 4, 6]
elbow_test = np.arange(80, 170, 1)
manus_test = np.arange(100, 185, 1)
today = date.today()
date_adj = today.strftime("%Y_%m_%d")

with open((date_adj+'_LongDynStability_Rigid.csv'), 'w', newline='') as res_file:
    writer = csv.writer(res_file)
    writer.writerow(['date', 'alpha', 'Ve', 'elbow', 'manus', 'Iyy', 'theta0', 'eignum', 'zeta', 'omega_n',
                     'eig_real', 'eig_imag', 'mag1', 'mag2', 'mag3', 'mag4',
                     'phase1', 'phase2', 'phase3', 'phase4',
                     'CL', 'CD', 'CL_alp', 'CD_alp', 'Cm_alp', 'CL_q', 'Cm_q', 'trim_x', "trim_z", "trim_m"])

for e in elbow_test:
    for w in manus_test:

        # Step 1: Set morphological data
        elbow = e
        manus = w
        Iyy = dynfn.get_Iyy(elbow, manus, coef_data)

        # Step 2: Calculate the trim flight path angle and speed for steady gliding flight
        theta_rad, V_e, alpha_0 = aerofn.trim_aero(W, rho, S_max, elbow, manus, coef_data)

        if V_e == "NA":
            continue

        # Step 2: Solve the homogenous eigen problem
        A = dynfn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, alpha_0, V_e,
                               theta_rad, coef_data)

        # initialize to not compute time series unless specified
        compute = False

        # Only grab time series of the following configurations
        if e == 90 and w == 120:
            compute = True
        if e == 120 and w == 120:
            compute = True
        if e == 120 and w == 157:
                compute = True
        if compute:
            # Step 3: Extract the free response of the system wrt d_alp
            x0 = [0, 5 * np.pi / 180, 0, 0]
            out = dynfn.calc_res_dalp(x0, 60, 6000, A)

            filename = "outputdata/" + \
                        date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                        "_alpha" + str(alpha_0) + "_dalp.csv"
            np.savetxt(filename, out, delimiter=",")

            # Step 4: Extract the time response of the system wrt ramped up speed
            x0 = [0, 0, 0, 0]
            B = np.asarray([0, 0.05, 0, 0])   # time dependent portion of the forced response
            C = np.asarray([0, 0, 0, 0])  # constant in time portion of the forced response

            ramp_end = 30  # seconds
            length_t = 120
            out = dynfn.calc_res_uramp(x0, length_t, ramp_end*200, A, B, C, ramp_end)

            filename = "outputdata/" + \
                        date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                        "_alpha" + str(alpha_0) + "_uramp.csv"
            np.savetxt(filename, out, delimiter=",")

