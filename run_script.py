# ----- Calculate the longitudinal dynamic stability of rigid configuration -----

# -------- Import libraries ---------
import pandas as pd
import dyn_functions as dynfn
import numpy as np
import csv
import aero_functions as aerofn
from datetime import date
import control
import matplotlib.pyplot as plt

# -------- Generic run script -------
rho = 1.225
# taken from the bird body specimen tested in AvInertia - chord was scaled to match this bird
m = 1.0154  # mass (kg)
W = m * 9.81  # weight (N)
S_max = 0.244657  # wings and body reference area
c_max = 0.2861011  # wing root chord

coef_data = pd.read_csv('/Volumes/GoogleDrive/My Drive/DoctoralThesis/Chapter3_DynamicStability/coefficients'
                                '/2021_11_21_coefficients_all.csv')

dihedral_test = [10, 15, 20]
sweep_test = [-20, -15, -10, -5, 0, 5, 10]
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

        for e in elbow_test:
            for w in manus_test:

                # Step 1: Set morphological data
                elbow = e
                manus = w
                Iyy = dynfn.get_Iyy(elbow, manus, sw, di, coef_data)

                # Step 2: Calculate the trim flight path angle and speed for steady gliding flight
                gamma_rad, V_e, alpha_0 = aerofn.trim_aero(W, rho, S_max, elbow, manus, sw, di, coef_data)

                if V_e == "NA":
                    continue

                # Step 2: Solve the homogenous eigen problem
                A, B_gust = dynfn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, sw, di, alpha_0, V_e,
                                               gamma_rad, coef_data)

                # initialize to not compute time series unless specified
                compute = False

                # Only grab time series of the following configurations
                if e == 126 and (w == 106 or w == 116 or w == 126 or w == 136 or w == 146 or w == 156):
                    compute = True

                if compute:

                    # --------------------------------------------------------------------
                    # Step 3: Extract the free response of the system wrt d_alp
                    B = np.transpose(np.asarray([[0, 0, 0, 0]]))
                    C = np.identity(4)
                    D = B
                    t = [np.arange(0, 100, step=0.001)]
                    # Define the state space model
                    system = control.ss(A, B, C, D)
                    x0 = np.transpose([[0, 2 * np.pi / 180, 0, 0]]) # increase of 2deg
                    # solve the system
                    out_step = control.initial_response(system, T=t, X0=x0)
                    # save the data
                    dat_step = np.transpose(np.insert(out_step[1], values=np.asarray(out_step[0]), obj=0, axis=0))
                    filename = "outputdata/" + \
                                date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                               "_sw" + str(sw) + "_di" + str(di) + "_dalp.csv"
                    np.savetxt(filename, dat_step, delimiter=",")
                    # --------------------------------------------------------------------
                    # Step 4: Extract the time response of the system wrt to a ramping incoming velocity
                    # allows speed control
                    system = control.ss(A, B_gust, C, D)
                    x0 = np.transpose([[0, 0, 0, 0]])  # initial condition
                    # ramps speed linearly with time up to 5 secs
                    u_ramp = np.append(0.01*t[0][0:5001], 0.01*t[0][5000]*np.ones(len(t[0]) - len(t[0][0:5001])))
                    # solve the system
                    out_force = control.forced_response(system, T=t, X0=x0, U=u_ramp)
                    # save the data
                    dat_force = np.transpose(np.insert(out_force[1], values=np.asarray(out_force[0]), obj=0, axis=0))
                    filename = "outputdata/" + \
                                date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                               "_sw" + str(sw) + "_di" + str(di) + "_uramp.csv"
                    np.savetxt(filename, dat_force, delimiter=",")

                    # Step 5: Extract the time response of the system wrt to a ramping incoming velocity
                    # NO speed control
                    A_no_u = np.delete(np.delete(A, obj=0, axis=0), obj=0, axis=1)
                    B_no_u = np.delete(B_gust, obj=0)
                    B_no_u.shape = (3, 1)
                    D_no_u = np.zeros(3)
                    D_no_u.shape = (3, 1)
                    x0_no_u = np.transpose([[0, 0, 0]])  # initial condition
                    system = control.ss(A_no_u, B_no_u, np.identity(3), D_no_u)
                    out_force_no_u = control.forced_response(system, T=t, X0=x0_no_u, U=u_ramp)
                    # save the data
                    dat_force_no_u = np.insert(out_force_no_u[1],
                                               values=np.asarray(u_ramp), obj=0, axis=0)
                    dat_force_no_u = np.transpose(np.insert(dat_force_no_u,
                                                            values=np.asarray(out_force_no_u[0]), obj=0, axis=0))
                    filename = "outputdata/" + \
                                date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                               "_sw" + str(sw) + "_di" + str(di) + "_uramp_no_u.csv"
                    np.savetxt(filename, dat_force_no_u, delimiter=",")
                    # --------------------------------------------------------------------
