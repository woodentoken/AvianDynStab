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
W = m*9.81
S_max = 0.1074285
c_max = 0.2861011
coef_data = pd.read_csv('/Users/christinaharvey/Google Drive/DoctoralThesis/Chapter3 - DynamicStability'
                        '/2021_07_12_coefficients.csv')
a_test = [-2, 0, 2, 4, 6]
U_test = [10, 15, 20, 25]
elbow_test = np.arange(80, 155, 5)
manus_test = np.arange(100, 170, 5)

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
                damp,freq = aerofn.solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, alpha_0, U_0, theta_0, coef_data)

                today = date.today()
                date_adj = today.strftime("%Y_%m_%d")

                with open('LongDynStability_Rigid.csv', 'a', newline="") as res_file:
                    writer = csv.writer(res_file)
                    writer.writerow([date_adj, i, j, e, w, Iyy, theta_0, damp[0], damp[1], freq[0], freq[1]])
