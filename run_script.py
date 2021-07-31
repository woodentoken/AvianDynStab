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
coef_data = pd.read_csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability'
                        '/2021_07_29_coefficients.csv')
a_test = [-2, 0, 2, 4, 6]
elbow_test = np.arange(80, 170, 5)
manus_test = np.arange(100, 185, 5)

with open('LongDynStability_Rigid.csv', 'w', newline='') as res_file:
    writer = csv.writer(res_file)
    writer.writerow(['date', 'alpha', 'U0', 'elbow', 'manus', 'Iyy', 'theta0', 'eignum', 'zeta', 'omega_n',
                     'eig_real', 'eig_imag', 'mag1', 'mag2', 'mag3', 'mag4',
                     'phase1', 'phase2', 'phase3', 'phase4',
                     'CL', 'CD', 'CL_alp', 'CD_alp', 'Cm_alp', 'CL_q', 'Cm_q', 'del_x', 'del_z'])

for i in a_test:
    for e in elbow_test:
        for w in manus_test:
            # Step 1: Set flight conditions
            alpha_0 = i

            # Step 2: Set morphological data
            elbow = e
            manus = w
            Iyy = coef_data['elbow'][0] * elbow + coef_data['manus'][0] * manus + \
                  coef_data['elbow2'][0] * elbow ** 2 + coef_data['manus2'][0] * manus ** 2 + \
                  coef_data['elbow3'][0] * elbow ** 3 + coef_data['manus3'][0] * manus ** 3 + \
                  coef_data['elbowmanus'][0] * elbow * manus + coef_data['intercept'][0]

            # Step 1: Calculate the trim flight path angle and speed for steady gliding flight
            theta_0, U_0, del_x, del_z = aerofn.trim_aero(W, rho, S_max, c_max, elbow, manus, alpha_0, coef_data)

            if U_0 == 0:
                continue

            # Step 2: Solve the eigen problem
            A = dynfn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, alpha_0, U_0,
                                   theta_0, coef_data, del_x, del_z)

            # Step 3: Extract the time response of the system
            x0 = [0, 5*np.pi/180, 0, 0]
            t = np.linspace(0, 60, 5000)
            x = dynfn.solve_IVP(A, x0, t)

            # Step 4: arrange the data to be saved
            t.shape = (5000, 1)
            out = np.hstack((t, x))
            today = date.today()
            date_adj = today.strftime("%Y_%m_%d")

            filename = "outputdata/" + \
                       date_adj + "_elbow" + str(elbow) + "_manus" + str(manus) + \
                       "_alpha" + str(alpha_0) + ".csv"
            np.savetxt(filename, out, delimiter=",")


