import numpy as np
import cmath
from datetime import date
import csv
import aero_functions as aerofn
from scipy.integrate import odeint


def solve_linsys(m, Iyy, rho, S, c, elbow, manus, sw, di, alpha_0, V_e, theta_rad, aero_data, coef_q_data):
    # Caution: theta_rad is in rad and alpha_0 is in deg
    # Collect all necessary aerodynamic coefficients
    CL = aerofn.get_CL(aero_data, elbow, manus, alpha_0)
    CD = aerofn.get_CD(aero_data, elbow, manus, CL)

    # Collect all angle of attack derivatives
    CL_alp = aerofn.get_dCL_dalp(aero_data, elbow, manus, alpha_0)
    CD_alp = aerofn.get_dCD_dalp(aero_data, elbow, manus, CL, CL_alp)
    Cm_alp = aerofn.get_dCm_dCL(aero_data, elbow, manus)*CL_alp

    # Collect all pitch rate derivatives
    CL_q = aerofn.get_dCL_dq(coef_q_data, elbow, manus, sw, di)
    Cm_q = aerofn.get_dCm_dq(coef_q_data, elbow, manus, sw, di)

    # define common constants
    alpha_rad = np.deg2rad(alpha_0)
    m_til = rho*V_e*S/(2*m)
    I_til = rho*(V_e**2)*S*c/(2*Iyy)
    w_til = 9.81/V_e
    # just for saving purposes
    trim = aerofn.trim_force(np.array([V_e, theta_rad]), m*9.81, rho, S, CL, CD, alpha_rad)
    Cm = aerofn.get_Cm(aero_data, elbow, manus, CL)

    # define the linear system
    A = np.array([[2 * m_til * (CL*np.sin(alpha_rad)-CD*np.cos(alpha_rad)),
                   m_til * ((CL - CD_alp)*np.cos(alpha_rad) + (CD + CL_alp)*np.sin(alpha_rad)),
                   (m_til*CL_q+1)*np.sin(alpha_rad),
                   -w_til*np.cos(theta_rad)],
                  [2 * m_til * (-CL*np.cos(alpha_rad) - CD*np.sin(alpha_rad)),
                   m_til * ((CL - CD_alp)*np.sin(alpha_rad) - (CD + CL_alp)*np.cos(alpha_rad)),
                   (m_til*CL_q+1)*np.cos(alpha_rad),
                   -w_til*np.sin(theta_rad)],
                  [0, I_til * Cm_alp, I_til * Cm_q, 0],
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
    with open((date_adj+'_LongDynStability_Rigid.csv'), 'a', newline="") as res_file:

        writer = csv.writer(res_file)
        eignum = 1
        writer.writerow([date_adj, alpha_0, V_e, elbow, manus, sw, di, Iyy, theta_rad, eignum, damp[0], freq[0],
                         eig_val[0].real, eig_val[0].imag, mag[0], mag[1], mag[2], mag[3],
                         phase[0], phase[1], phase[2], phase[3],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 2
        writer.writerow([date_adj, alpha_0, V_e, elbow, manus, sw, di, Iyy, theta_rad, eignum, damp[0], freq[0],
                         eig_val[1].real, eig_val[1].imag, mag[4], mag[5], mag[6], mag[7],
                         phase[4], phase[5], phase[6], phase[7],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 3
        writer.writerow([date_adj, alpha_0, V_e, elbow, manus, sw, di, Iyy, theta_rad, eignum, damp[1], freq[1],
                         eig_val[2].real, eig_val[2].imag, mag[8], mag[9], mag[10], mag[11],
                         phase[8], phase[9], phase[10], phase[11],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])
        eignum = 4
        writer.writerow([date_adj, alpha_0, V_e, elbow, manus, sw, di, Iyy, theta_rad, eignum, damp[1], freq[1],
                         eig_val[3].real, eig_val[3].imag, mag[12], mag[13], mag[14], mag[15],
                         phase[12], phase[13], phase[14], phase[15],
                         CL, CD, CL_alp, CD_alp, Cm_alp, CL_q, Cm_q, trim, Cm])

    return A


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


def get_Iyy(elbow, manus, coef_data):
    Iyy = coef_data['elbow'][0] * elbow + coef_data['manus'][0] * manus + \
          coef_data['elbow2'][0] * elbow ** 2 + coef_data['manus2'][0] * manus ** 2 + \
          coef_data['elbow3'][0] * elbow ** 3 + coef_data['manus3'][0] * manus ** 3 + \
          coef_data['elbowmanus'][0] * elbow * manus + coef_data['intercept'][0]

    return Iyy


def calc_res_dalp(x0, t_max, dt, A):

    t = np.linspace(0, t_max, dt)
    x = solve_IVP(A, x0, t, mod_type="free")

    # arrange the data to be saved
    t.shape = (dt, 1)
    out_dalp = np.hstack((t, x))

    return out_dalp


def calc_res_uramp(x0, t_max, dt, A, B, C, end):
    t = np.linspace(0, t_max, dt)
    x = solve_IVP(A, x0, t, B, C, end, mod_type="forced")

    # arrange the data to be saved
    t.shape = (dt, 1)
    out_du = np.hstack((t, x))

    return out_du