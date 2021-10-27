# -- This file contains all functions that are required to properly calculate the trim position --
import numpy as np
import scipy.optimize as opt


def trim_aero(W, rho, S, elbow, manus, aero_data):
    # calculate the necessary flight path angle and speed for trim
    pitch_guess = np.array([0])
    pitch_sol = opt.fsolve(trim_pitch, pitch_guess,
                           args=(elbow, manus, aero_data))
    alpha_0 = pitch_sol[0]
    alpha_rad = np.deg2rad(alpha_0)

    CL = get_CL(aero_data, elbow, manus, alpha_0)
    CD = get_CD(aero_data, elbow, manus, CL)
    # solve the force equations for V and theta at this alpha
    force_guess = np.array([10, 0])
    force_sol = opt.fsolve(trim_force, force_guess,
                           args=(W, rho, S, CL, CD, alpha_rad))

    # save outputs note that these constraints are used in the solving which is why the numbers are saved like this
    V_e = abs(force_sol[0])
    theta_rad = -abs(force_sol[1])

    if theta_rad < -0.5*np.pi:
        theta_rad = -0.5*np.pi

    # make sure this result actually solved for an equilibrium
    out_forces = trim_force(np.array([V_e, theta_rad]), W, rho, S, CL, CD, alpha_rad)
    if abs(out_forces[0]) > 10e-4 or abs(out_forces[0]) > 10e-4:
        theta_rad = "NA"
        V_e = "NA"

    return theta_rad, V_e, alpha_0


def get_CL(aero_data, elbow, manus, alpha_0):
    # incoming alpha should be in degrees
    CL = aero_data['elbow'][2]*elbow + aero_data['manus'][2]*manus + \
         aero_data['elbow2'][2]*(elbow**2) + aero_data['manus2'][2]*(manus**2) + \
         aero_data['elbow3'][2]*(elbow**3) + aero_data['manus3'][2]*(manus**3) + \
         aero_data['elbowmanus'][2]*elbow*manus + aero_data['intercept'][2] + \
         aero_data['alpha'][2]*alpha_0 + aero_data['alpha2'][2]*(alpha_0**2) + \
         aero_data['alpha3'][2]*(alpha_0**3) + aero_data['elbowalpha'][2]*elbow*alpha_0 + \
         aero_data['manusalpha'][2]*manus*alpha_0 + aero_data['elbowmanusalpha'][2]*elbow*manus*alpha_0

    return CL


def get_dCL_dalp(aero_data, elbow, manus, alpha_0):
    # incoming alpha should be in degrees
    CL_alp = aero_data['alpha'][2] + 2 * aero_data['alpha2'][2] * alpha_0 + \
             3 * aero_data['alpha3'][2] * (alpha_0 ** 2) + aero_data['elbowalpha'][2] * elbow + \
             aero_data['manusalpha'][2] * manus + aero_data['elbowmanusalpha'][2] * elbow * manus

    # CAUTION THE INPUT ANGLE OF ATTACK IS IN DEGREES
    CL_alp = CL_alp/(np.pi/180)
    # output will be in radians
    return CL_alp


def get_CD(aero_data, elbow, manus, CL):
    CD = aero_data['elbow'][4]*elbow + aero_data['manus'][4]*manus + \
         aero_data['intercept'][4] + aero_data['CL'][4]*CL + \
         aero_data['CL2'][4]*(CL**2) + aero_data['CL3'][4]*(CL**3) + \
         aero_data['elbowCL'][4]*elbow*CL + aero_data['manusCL'][4]*manus*CL + \
         aero_data['elbowCL2'][4] * elbow * (CL**2) + aero_data['manusCL2'][4] * manus * (CL**2)

    return CD


def get_dCD_dalp(aero_data, elbow, manus, CL, CL_alp):
    cdcl = aero_data['CL'][4] + 2 * aero_data['CL2'][4] * CL + \
           3 * aero_data['CL3'][4] * (CL ** 2) + aero_data['elbowCL'][4] * elbow + \
           aero_data['manusCL'][4] * manus + aero_data['elbowmanusCL'][4] * elbow * manus + \
           2*aero_data['elbowCL2'][4] * elbow * CL + 2*aero_data['manusCL2'][4] * manus * CL

    CD_alp = cdcl * CL_alp
    # output will be in radians because of adjustment in get_dCL_dalp
    return CD_alp


def get_Cm(aero_data, elbow, manus, CL):
    Cm = aero_data['elbow'][3]*elbow + aero_data['manus'][3]*manus + \
         aero_data['elbow2'][3]*(elbow**2) + aero_data['manus2'][3]*(manus**2) + \
         aero_data['elbow3'][3]*(elbow**3) + aero_data['manus3'][3]*(manus**3) + \
         aero_data['elbowmanus'][3]*elbow*manus + aero_data['intercept'][3] + \
         aero_data['CL'][3]*CL + aero_data['CL2'][3]*(CL**2) + \
         aero_data['CL3'][3]*(CL**3) + aero_data['elbowCL'][3]*elbow*CL + \
         aero_data['manusCL'][3]*manus*CL + aero_data['elbowmanusCL'][3]*elbow*manus*CL

    return Cm


def get_dCm_dalp(aero_data, elbow, manus, CL_alp):
    cmcl = aero_data['elbow'][1] * elbow + aero_data['manus'][1] * manus + \
           aero_data['elbow2'][1] * (elbow ** 2) + aero_data['manus2'][1] * (manus ** 2) + \
           aero_data['elbow3'][1] * (elbow ** 3) + aero_data['manus3'][1] * (manus**3) +  \
           aero_data['elbowmanus'][1] * elbow * manus + aero_data['intercept'][1]

    Cm_alp = cmcl * CL_alp
    # output will be in radians because of adjustment in get_dCL_dalp
    return Cm_alp


def get_dCL_dq(aero_data, elbow, manus):
    CL_q = aero_data['elbow'][5]*elbow + aero_data['manus'][5]*manus + \
           aero_data['elbow2'][5] * elbow ** 2 + aero_data['manus2'][5] * manus ** 2 + \
           aero_data['elbowmanus'][5] * elbow * manus + aero_data['intercept'][5]
    # output will be in radians because incoming q is in radians
    return CL_q


def get_dCm_dq(aero_data, elbow, manus):
    Cm_q = aero_data['elbow'][6]*elbow + aero_data['manus'][6]*manus + \
           aero_data['elbow2'][6] * elbow ** 2 + aero_data['manus2'][6] * manus ** 2 + \
           aero_data['elbowmanus'][6] * elbow * manus + aero_data['intercept'][6]
    # output will be in radians because incoming q is in radians
    return Cm_q


def min_geo_angle(x, Cm, CL, CD, alpha, c_root):

    z = (c_root/(CL*np.sin(alpha)-CD*np.cos(alpha)))*(Cm - ((CL*np.cos(alpha)+CD*np.sin(alpha))*(x/c_root)))
    out = x**2 + z**2

    return out


def trim_force(x, W, rho, S, CL, CD, alpha_rad):
    V = abs(x[0])
    theta_rad = -abs(x[1])

    if theta_rad < -0.5*np.pi:
        theta_rad = -0.5*np.pi

    F = np.empty(2)
    # x direction forces
    F[0] = 0.5 * rho * V ** 2 * S * (CL * np.sin(alpha_rad) - CD * np.cos(alpha_rad)) - W * np.sin(theta_rad)
    # z direction forces
    F[1] = 0.5 * rho * V ** 2 * S * (-CL * np.cos(alpha_rad) - CD * np.sin(alpha_rad)) + W * np.cos(theta_rad)

    return F


def trim_pitch(x, elbow, manus, aero_data):
    alpha_0 = x
    CL = get_CL(aero_data, elbow, manus, alpha_0)
    Cm = get_Cm(aero_data, elbow, manus, CL)

    return Cm