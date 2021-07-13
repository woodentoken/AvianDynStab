# ----- Calculate the longitudinal dynamic stability of rigid configuration -----

# -------- Import libraries ---------
import pandas as pd
import aero_functions as aerofn
import numpy as np
import csv
from datetime import date

# -------- Generic run script -------
rho = 1.225
m = 1.0154
W = m * 9.81
S_max = 0.1074285
c_max = 0.2861011
coef_data = pd.read_csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3_DynamicStability'
                        '/2021_07_13_coefficients.csv')
a_test = [-2, 0, 2, 4, 6]
U_test = [10, 15, 20, 25]
elbow_test = np.arange(80, 155, 5)
manus_test = np.arange(100, 170, 5)

with open('LongDynStability_Rigid.csv', 'w', newline='') as res_file:
    writer = csv.writer(res_file)
    writer.writerow(['date', 'alpha', 'U0', 'elbow', 'manus', 'Iyy', 'theta0', 'eignum', 'zeta', 'omega_n',
              'eig_real', 'eig_imag', 'mag1', 'mag2', 'mag3', 'mag4',
              'phase1', 'phase2', 'phase3', 'phase4'])

for i in a_test:
    for j in U_test:
        for e in elbow_test:
            for w in manus_test:
                # Step 1: Set flight conditions
                alpha_0 = i
                U_0 = j

                # Step 2: Set morphological data
                elbow = e
                manus = w
                Iyy = coef_data['elbow'][0] * elbow + coef_data['manus'][0] * manus + \
                      coef_data['elbow2'][0] * elbow ** 2 + coef_data['manus2'][0] * manus ** 2 + \
                      coef_data['elbow3'][0] * elbow ** 3 + coef_data['manus3'][0] * manus ** 3 + \
                      coef_data['elbowmanus'][0] * elbow * manus + coef_data['intercept'][0]

                # Step 1: Calculate the trim position for steady gliding flight
                theta_0 = aerofn.trim_aero(W, elbow, manus, alpha_0, coef_data)

                # Step 2: Solve the eigen problem
                damp, freq, eig_val, mag, phase = aerofn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus,
                                                                      alpha_0, U_0, theta_0, coef_data)

                today = date.today()
                date_adj = today.strftime("%Y_%m_%d")

                with open('LongDynStability_Rigid.csv', 'a', newline="") as res_file:

                    writer = csv.writer(res_file)
                    eignum = 1
                    writer.writerow([date_adj, i, j, e, w, Iyy, theta_0, eignum, damp[0], freq[0],
                                     eig_val[0].real, eig_val[0].imag, mag[0], mag[1], mag[2], mag[3],
                                     phase[0], phase[1], phase[2], phase[3]])
                    eignum = 2
                    writer.writerow([date_adj, i, j, e, w, Iyy, theta_0, eignum, damp[0], freq[0],
                                     eig_val[1].real, eig_val[1].imag, mag[4], mag[5], mag[6], mag[7],
                                     phase[4], phase[5], phase[6], phase[7]])
                    eignum = 3
                    writer.writerow([date_adj, i, j, e, w, Iyy, theta_0, eignum, damp[1], freq[1],
                                     eig_val[2].real, eig_val[2].imag, mag[8], mag[9], mag[10], mag[11],
                                     phase[8], phase[9], phase[10], phase[11]])
                    eignum = 4
                    writer.writerow([date_adj, i, j, e, w, Iyy, theta_0, eignum, damp[1], freq[1],
                                     eig_val[3].real, eig_val[3].imag, mag[12], mag[13], mag[14], mag[15],
                                     phase[12], phase[13], phase[14], phase[15]])
