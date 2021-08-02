import numpy as np
import cmath
from datetime import date
import csv
import aero_functions as aerofn
from scipy.integrate import odeint


def solve_linsys(m, Iyy, rho, S_max, c_max, elbow, manus, alpha_0, U_0, theta_0, aero_data, del_x, del_z):

    # Collect all necessary aerodynamic coefficients
    CL = aerofn.get_CL(aero_data, elbow, manus, alpha_0)
    CD = aerofn.get_CD(aero_data, elbow, manus, CL)
    # Collect all angle of attack derivatives
    CL_alp = aerofn.get_dCL_dalp(aero_data, elbow, manus, alpha_0)
    CD_alp = aerofn.get_dCD_dalp(aero_data, elbow, manus, CL, CL_alp)
    Cm_alp = aerofn.get_dCm_dalp(aero_data, elbow, manus, CL_alp)
    # Collect all pitch rate derivatives
    CL_q = aerofn.get_dCL_dq(aero_data, elbow, manus)
    Cm_q = aerofn.get_dCm_dq(aero_data, elbow, manus)
    # define common constants
    m1_inv = 1/(2*m/(rho*U_0*S_max))
    Iyy1_inv = 1/(2*Iyy/(rho*(U_0**2)*S_max*c_max))
    w1 = m*9.81/(0.5*rho*(U_0**2)*S_max)
    # define the linear system
    A = np.array([[m1_inv * (-2 * CD), m1_inv * (CL - CD_alp), 0, m1_inv * (-w1*np.cos(theta_0))],
                  [m1_inv * (-2 * CL), m1_inv * (-CD - CL_alp), 1+(-CL_q*m1_inv), m1_inv * (-w1*np.sin(theta_0))],
                  [0, Iyy1_inv * Cm_alp, Iyy1_inv * Cm_q, 0],
                  [0, 0, 1, 0]])
    # solve for the free response of the linear system
    eig_val, eig_vec = np.linalg.eig(A)

    # pre-define inputs
    damp = [0] * 2
    freq = [0] * 2
    mag = [0] * 16
    phase = [0] * 16

    # Calculate the frequency and damping for imaginary root
    if abs(eig_val[0].imag) > 0 and abs(eig_val[1].imag) > 0:
        freq[0] = np.sqrt(eig_val[0].real ** 2 + eig_val[0].imag ** 2)
        damp[0] = - eig_val[0].real / freq[0]

    if abs(eig_val[2].imag) > 0 and abs(eig_val[3].imag) > 0:
        freq[1] = np.sqrt(eig_val[2].real ** 2 + eig_val[2].imag ** 2)
        damp[1] = - eig_val[2].real / freq[1]

    # Calculate the phase and offset of each eigenvector
    count = 0
    for it1 in range(0, 4, 1):
        for it2 in range(0, 4, 1):
            phase[count] = cmath.phase(eig_vec[it2, it1])
            mag[count] = abs(eig_vec[it2, it1])
            count = count + 1

    # Pull in todays date
    today = date.today()
    date_adj = today.strftime("%Y_%m_%d")

    # Save data
    with open('LongDynStability_Rigid.csv', 'a', newline="") as res_file:

        writer = csv.writer(res_file)
        eignum = 1
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, Iyy, theta_0, eignum, damp[0], freq[0],
                         eig_val[0].real, eig_val[0].imag, mag[0], mag[1], mag[2], mag[3],
                         phase[0], phase[1], phase[2], phase[3],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, del_x, del_z])
        eignum = 2
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, Iyy, theta_0, eignum, damp[0], freq[0],
                         eig_val[1].real, eig_val[1].imag, mag[4], mag[5], mag[6], mag[7],
                         phase[4], phase[5], phase[6], phase[7],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, del_x, del_z])
        eignum = 3
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, Iyy, theta_0, eignum, damp[1], freq[1],
                         eig_val[2].real, eig_val[2].imag, mag[8], mag[9], mag[10], mag[11],
                         phase[8], phase[9], phase[10], phase[11],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, del_x, del_z])
        eignum = 4
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, Iyy, theta_0, eignum, damp[1], freq[1],
                         eig_val[3].real, eig_val[3].imag, mag[12], mag[13], mag[14], mag[15],
                         phase[12], phase[13], phase[14], phase[15],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, del_x, del_z])

    return A


def def_model(x, t, A):
    dxdt = A.dot(x)
    return dxdt


def solve_IVP(A, x0, t):
    x = odeint(def_model, x0, t, args=(A,))
    return x


def get_Iyy(elbow, manus, coef_data):
    Iyy = coef_data['elbow'][0] * elbow + coef_data['manus'][0] * manus + \
          coef_data['elbow2'][0] * elbow ** 2 + coef_data['manus2'][0] * manus ** 2 + \
          coef_data['elbow3'][0] * elbow ** 3 + coef_data['manus3'][0] * manus ** 3 + \
          coef_data['elbowmanus'][0] * elbow * manus + coef_data['intercept'][0]

    return Iyy


def U_ramp_par(G, A, t):

    # assumes a shape of a*t + b
    a = np.linalg.solve(A, -G)
    b = np.linalg.solve(A, a)

    x_p = a*t+b

    return x_p