# -- This file contains all functions that are required to properly calculate the trim position --
import numpy as np


def trim_aero(W, rho, S, elbow, manus, alpha_0, aero_data):
    U_0 = 0
    # use lift equation to solve for theta
    CL = get_CL(aero_data, elbow, manus, alpha_0)
    if CL < 0:
        neg_CL = 1
    CD = get_CD(aero_data, elbow, manus, CL)
    # calculate the necessary flight path angle for trim
    theta_0 = np.arctan(-CD/CL)
    # calculate the necessary speed for trim
    if CL > 0:
        U_0 = np.sqrt(2*W*np.cos(theta_0)/(rho*S*CL))

    return theta_0, U_0


def get_CL(aero_data, elbow, manus, alpha_0):
    CL = aero_data['elbow'][2]*elbow + aero_data['manus'][2]*manus + \
         aero_data['elbow2'][2]*elbow**2 + aero_data['manus2'][2]*manus**2 + \
         aero_data['manus3'][2]*manus**3 + \
         aero_data['elbowmanus'][2]*elbow*manus + aero_data['intercept'][2] + \
         aero_data['alpha'][2]*alpha_0 + aero_data['alpha2'][2]*alpha_0**2 + \
         aero_data['alpha3'][2]*alpha_0**3 + aero_data['elbowalpha'][2]*elbow*alpha_0 + \
         aero_data['manusalpha'][2]*manus*alpha_0 + aero_data['elbowmanusalpha'][2]*elbow*manus*alpha_0

    return CL


def get_dCL_dalp(aero_data, elbow, manus, alpha_0):
    CL_alp = aero_data['alpha'][2] + 2 * aero_data['alpha2'][2] * alpha_0 + \
             3 * aero_data['alpha3'][2] * alpha_0 ** 2 + aero_data['elbowalpha'][2] * elbow + \
             aero_data['manusalpha'][2] * manus + aero_data['elbowmanusalpha'][2] * elbow * manus

    # CAUTION THE INPUT ANGLE OF ATTACK IS IN DEGREES
    CL_alp = CL_alp/(np.pi/180)
    # output will be in radians
    return CL_alp


def get_CD(aero_data, elbow, manus, CL):
    CD = aero_data['elbow'][4]*elbow + aero_data['manus'][4]*manus + \
         aero_data['elbow2'][4]*elbow**2 + aero_data['manus2'][4]*manus**2 + \
         aero_data['manus3'][4]*manus**3 + \
         aero_data['elbowmanus'][4]*elbow*manus + aero_data['intercept'][4] + \
         aero_data['CL'][4]*CL + aero_data['CL2'][4]*CL**2 + \
         aero_data['CL3'][4]*CL**3 + aero_data['elbowCL'][4]*elbow*CL + \
         aero_data['manusCL'][4]*manus*CL + aero_data['elbowmanusCL'][4]*elbow*manus*CL

    return CD


def get_dCD_dalp(aero_data, elbow, manus, CL, CL_alp):
    cdcl = aero_data['CL'][4] + 2 * aero_data['CL2'][4] * CL + \
           3 * aero_data['CL3'][4] * CL ** 2 + aero_data['elbowCL'][4] * elbow + \
           aero_data['manusCL'][4] * manus + aero_data['elbowmanusCL'][4] * elbow * manus

    CD_alp = cdcl * CL_alp
    # output will be in radians
    return CD_alp


def get_dCm_dalp(aero_data, elbow, manus, CL_alp):
    cmcl = aero_data['elbow'][1] * elbow + aero_data['manus'][1] * manus + \
           aero_data['elbow2'][1] * elbow ** 2 + aero_data['manus2'][1] * manus ** 2 + \
           aero_data['elbow3'][1] * elbow ** 3 + \
           aero_data['elbowmanus'][1] * elbow * manus + aero_data['intercept'][1]

    Cm_alp = cmcl * CL_alp
    # output will be in radians
    return Cm_alp


def get_dCL_dq(aero_data, elbow, manus):
    CL_q = aero_data['elbow'][5]*elbow + aero_data['manus'][5]*manus + \
           aero_data['elbowmanus'][5] * elbow * manus + aero_data['intercept'][5]
    # output will be in radians
    return CL_q


def get_dCm_dq(aero_data, elbow, manus):
    Cm_q = aero_data['elbow'][6]*elbow + aero_data['manus'][6]*manus + \
           aero_data['elbowmanus'][6] * elbow * manus + aero_data['intercept'][6]
    # output will be in radians
    return Cm_q