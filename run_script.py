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

coef_q_data = pd.read_csv('/Volumes/GoogleDrive/My Drive/DoctoralThesis/Chapter3_DynamicStability/coefficients'
                                '/2021_11_03_coefficients_q.csv')

a_test = [-2, 0, 2, 4, 6]
dihedral_test = [0, 10, 20]
sweep_test = [-20, -10, 0, 10, 20]
elbow_test = np.arange(86, 166, 2)
manus_test = np.arange(106, 180, 2)

today = date.today()
date_adj = today.strftime("%Y_%m_%d")

with open((date_adj+'_LongDynStability_Rigid.csv'), 'w', newline='') as res_file:
    writer = csv.writer(res_file)
    writer.writerow(['date', 'alpha', 'U0', 'elbow', 'manus', 'sweep', 'dihedral', 'Iyy', 'gamma0',
                     'eignum', 'zeta', 'omega_n',
                     'eig_real', 'eig_imag', 'mag1', 'mag2', 'mag3', 'mag4',
                     'phase1', 'phase2', 'phase3', 'phase4',
                     'CL', 'CD', 'CL_alp', 'CD_alp', 'Cm_alp', 'CL_q', 'Cm_q', 'trim_force', "Cm"])
for sw in sweep_test:
    for di in dihedral_test:

        coef_data = pd.read_csv('/Volumes/GoogleDrive/My Drive/DoctoralThesis/Chapter3_DynamicStability/coefficients'
                                '/2021_11_03_coefficients_sw' + str(sw) + '_di' + str(di) + '.csv')

        for e in elbow_test:
            for w in manus_test:

                # Step 1: Set morphological data
                elbow = e
                manus = w
                Iyy = dynfn.get_Iyy(elbow, manus, coef_data)

                # Step 2: Calculate the trim flight path angle and speed for steady gliding flight
                gamma_rad, V_e, alpha_0 = aerofn.trim_aero(W, rho, S_max, elbow, manus, coef_data)

                if V_e == "NA":
                    continue

                # Step 2: Solve the homogenous eigen problem
                A = dynfn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, sw, di, alpha_0, V_e,
                                       gamma_rad, coef_data, coef_q_data)

                # initialize to not compute time series unless specified
                compute = False

                # Only grab time series of the following configurations
                if e == 90 and w == 120:
                    compute = True
                if e == 120 and w == 120:
                    compute = True
                if e == 120 and w == 150:
                    compute = True

                if compute:
                    # --------------------------------------------------------------------
                    # Step 3: Extract the free response of the system wrt d_alp
                    x0 = [0, 2 * np.pi / 180, 0, 0]     # initial condition
                    length_t = 60                       # total time saved
                    timesteps = 6000                    # amount of steps s.t. dt = timesteps/length_t
                    # solve for time response
                    out = dynfn.calc_res_dalp(x0, length_t, timesteps, A)
                    # save outputs
                    filename = "outputdata/" + \
                                date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                               "_sw" + str(sw) + "_di" + str(di) + "_dalp.csv"
                    np.savetxt(filename, out, delimiter=",")
                    # --------------------------------------------------------------------
                    # Step 4: Extract the time response of the system wrt ramped up speed
                    x0 = [0, 0, 0, 0]                   # initial condition
                    u_ramp = 0.01                       # change in the velocity divided by total velocity (1%)
                    B = np.asarray([u_ramp, 0, 0, 0])   # time dependent portion of the forced response
                    C = np.asarray([0, 0, 0, 0])        # constant in time portion of the forced response
                    ramp_end = 5                        # seconds until time dependent part stops
                    length_t = 60                       # total time saved
                    timesteps = ramp_end*1200           # amount of steps s.t. dt = timesteps/length_t
                    # solve for time response
                    out = dynfn.calc_res_uramp(x0, length_t, timesteps, A, B, C, ramp_end)
                    # save outputs
                    filename = "outputdata/" + \
                                date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                               "_sw" + str(sw) + "_di" + str(di) + "_uramp.csv"
                    np.savetxt(filename, out, delimiter=",")
                    # --------------------------------------------------------------------
