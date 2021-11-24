import copy

import numpy as np
import cmath
from datetime import date
import csv
import aero_functions as aerofn
from scipy.integrate import odeint


def solve_linsys(m, Iyy, rho, S, c, elbow, manus, sw, di, alpha_0, U_0, gamma_rad, aero_data):
    # Caution: theta_rad is in rad and alpha_0 is in deg
    # Collect all necessary aerodynamic coefficients
    CL = aerofn.get_CL(aero_data, elbow, manus, sw, di, alpha_0)
    CD = aerofn.get_CD(aero_data, elbow, manus, alpha_0)

    # Speed derivatives
    # Drag could extract from the experimental tests CDu = 5.593e-04 - choose to neglect due to small magnitude
    # There was no significant effect of speed on the pitching moment or lift

    # Collect all angle of attack derivatives
    CL_alp = aerofn.get_dCL_dalp(aero_data, elbow, manus, sw, di, alpha_0)
    CD_alp = aerofn.get_dCD_dalp(aero_data, elbow, manus, alpha_0)
    Cm_alp = aerofn.get_dCm_dCL(aero_data, elbow, manus, sw, di)*CL_alp

    # Collect all pitch rate derivatives
    CL_q = aerofn.get_dCL_dq(aero_data, elbow, manus, sw, di)
    Cm_q = aerofn.get_dCm_dq(aero_data, elbow, manus, sw, di, CL_q)

    # define common constants
    m_til = rho*U_0*S/(2*m)
    I_til = rho*(U_0**2)*S*c/(2*Iyy)
    w_til = 9.81/U_0

    # define the linear system
    A = np.array([[2 * m_til * (-CD),
                   m_til * (CL - CD_alp),
                   0,
                   -w_til*np.cos(gamma_rad)],
                  [2 * m_til * (-CL),
                   m_til * (-CL_alp - CD),
                   (1-(m_til*CL_q)),
                   -w_til*np.sin(gamma_rad)],
                  [0, I_til * Cm_alp, I_til * Cm_q, 0],
                  [0, 0, 1, 0]])

    # solve for the free response of the linear system
    eig_val, eig_vec = np.linalg.eig(A)

    # need to make sure that the complex conjugate pairs are grouped together
    new_eig_val = copy.deepcopy(eig_val)
    new_eig_vec = copy.deepcopy(eig_vec)
    if abs(eig_val[1].imag) == abs(eig_val[2].imag) and abs(eig_val[1].real) == abs(eig_val[2].real):
        new_eig_val[2] = copy.deepcopy(eig_val[1])
        new_eig_val[3] = copy.deepcopy(eig_val[2])
        new_eig_val[0] = copy.deepcopy(eig_val[0])
        new_eig_val[1] = copy.deepcopy(eig_val[3])
        new_eig_vec[2] = copy.deepcopy(eig_vec[1])
        new_eig_vec[3] = copy.deepcopy(eig_vec[2])
        new_eig_vec[0] = copy.deepcopy(eig_vec[0])
        new_eig_vec[1] = copy.deepcopy(eig_vec[3])

        eig_val = copy.deepcopy(new_eig_val)
        eig_vec = copy.deepcopy(new_eig_vec)

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
    # just for saving purposes
    trim = aerofn.trim_force(np.array([U_0, gamma_rad]), m*9.81, rho, S, CL, CD)
    Cm = aerofn.get_Cm(aero_data, elbow, manus, sw, di, CL)

    # Save data
    with open((date_adj+'_LongDynStability_Rigid.csv'), 'a', newline="") as res_file:

        writer = csv.writer(res_file)
        eignum = 1
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, sw, di, Iyy, gamma_rad, eignum, damp[0], freq[0],
                         eig_val[0].real, eig_val[0].imag, mag[0], mag[1], mag[2], mag[3],
                         phase[0], phase[1], phase[2], phase[3],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 2
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, sw, di, Iyy, gamma_rad, eignum, damp[0], freq[0],
                         eig_val[1].real, eig_val[1].imag, mag[4], mag[5], mag[6], mag[7],
                         phase[4], phase[5], phase[6], phase[7],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 3
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, sw, di, Iyy, gamma_rad, eignum, damp[1], freq[1],
                         eig_val[2].real, eig_val[2].imag, mag[8], mag[9], mag[10], mag[11],
                         phase[8], phase[9], phase[10], phase[11],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 4
        writer.writerow([date_adj, alpha_0, U_0, elbow, manus, sw, di, Iyy, gamma_rad, eignum, damp[1], freq[1],
                         eig_val[3].real, eig_val[3].imag, mag[12], mag[13], mag[14], mag[15],
                         phase[12], phase[13], phase[14], phase[15],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])

    B = np.array([2 * m_til * (-CD), 2 * m_til * (-CL), 0, 0])  # for gust response only

    return A, B


def def_model_free(x, t, A):
    dxdt = A.dot(x)
    return dxdt


def def_model_forced(x, t, A, B, C, end):
    dxdt = A.dot(x) + B*t + C
    if t > end:
        dxdt = A.dot(x) + B*end + C
    return dxdt


def solve_IVP(A, x0, t, B=None, C=None, end=None, mod_type="free"):
    if mod_type == "free":
        x = odeint(def_model_free, x0, t, args=(A,))
    if mod_type == "forced":
        # assumes that u(t) = t
        x = odeint(def_model_forced, x0, t, args=(A, B, C, end))
    return x


def calc_res_dalp(x0, t_max, timesteps, A):

    t = np.linspace(0, t_max, timesteps)
    x = solve_IVP(A, x0, t, mod_type="free")

    # arrange the data to be saved
    t.shape = (timesteps, 1)
    out_dalp = np.hstack((t, x))

    return out_dalp


def calc_res_uramp(x0, t_max, timesteps, A, B, C, end):
    t = np.linspace(0, t_max, timesteps)
    x = solve_IVP(A, x0, t, B, C, end, mod_type="forced")

    # arrange the data to be saved
    t.shape = (timesteps, 1)
    out_du = np.hstack((t, x))

    return out_du


def get_Iyy(elbow, manus, sweep, dihedral, coef_data):
    Iyy = coef_data['intercept'][0] + \
          coef_data['elbow'][0] * elbow + \
          coef_data['manus'][0] * manus + \
          coef_data['sweep'][0] * sweep + \
          coef_data['dihedral'][0] * dihedral + \
          coef_data['elbow2'][0] * elbow ** 2 + \
          coef_data['manus2'][0] * manus ** 2 + \
          coef_data['sweep2'][0] * sweep ** 2 + \
          coef_data['dihedral2'][0] * dihedral ** 2 + \
          coef_data['elbowmanus'][0] * elbow * manus + \
          coef_data['elbowsweep'][0] * elbow * sweep + \
          coef_data['manussweep'][0] * manus * sweep + \
          coef_data['elbowdihedral'][0] * elbow * dihedral + \
          coef_data['manusdihedral'][0] * manus * dihedral + \
          coef_data['sweepdihedral'][0] * sweep * dihedral + \
          coef_data['elbowmanussweep'][0] * elbow * manus * sweep + \
          coef_data['elbowmanusdihedral'][0] * elbow * manus * dihedral + \
          coef_data['elbowsweepdihedral'][0] * elbow * sweep * dihedral + \
          coef_data['manussweepdihedral'][0] * manus * sweep * dihedral + \
          coef_data['elbowmanussweepdihedral'][0] * elbow * manus * sweep * dihedral

    return Iyy