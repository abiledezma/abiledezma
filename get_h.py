import math as m


def calc_discharge(b, h, m_bank, S, **kwargs):
    A=h*(b+h*m_bank)
    P=b+2*h*(m.sqrt(m_bank**2 +1))
    Rh=A/P
    for k in kwargs.items():
        if "nm" in k[0]:
            k_st =1/float(k[1])
            print("Using k_st= " + str (k_st ))
        if "D" in k[0]:
            k_st=26/(float(k[1])**(1/6))
            print("Using k_st= " + str (k_st ))
        elif "k" in k[0]:
            k_st=float(k[1])
            print("Using k_st= " + str (k_st ))

    Q=k_st* m.sqrt(S) * (Rh**(2/3))*A
    return Q

def interpolate_h(Q, b, m_bank, S, **kwargs):
    h = 1.0
    eps = 1.0
    count=1

    for k in kwargs.items():
        if "nm" in k[0]:
            k_st =1/float(k[1])
            print("Using k_st= " + str (k_st ))
        if "D" in k[0]:
            k_st=26/(float(k[1])**(1/6))
            print("Using k_st= " + str (k_st ))
        elif "k" in k[0]:
            k_st=float(k[1])
            print("Using k_st= " + str (k_st ))
    n_m =1/float(k_st)

    while eps > 10**-3 and count < 5000:

        A=h*(b+h*m_bank)
        P=b+2*h*(m.sqrt(m_bank**2 +1))
        Qk = A**(5/3)*m.sqrt(S_0)/(n_m*P**(2/3))
        eps = abs(Q-Qk)/Q
        dA_dh = b+2*m_bank*h
        dP_dh = 2*m.sqrt(m_bank**2+1)
        F = n_m*Q*P**(2/3)-A**(5/3)*m.sqrt(S_0)
        dF_dh = 2/3*n_m*Q*P**(-1/3)*dP_dh-5/3*A**(2/3)*m.sqrt(S_0)*dA_dh
        h = abs(h-F/dF_dh)
        count = count+1

    return [float(h),float(eps)]

if __name__ == '__main__':
    # input parameters
    Q = 15.5        # discharge in (m3/s)
    b = 5.1         # bottom channel width (m)
    m_bank = 2.5    # bank slope
    k_st = 20       # Strickler value
    n_m = 1 / k_st  # Manning's n
    S_0 = 0.005     # channel slope


    # call the solver with user-defined channel geometry and discharge
    Qo = calc_discharge(b, 1, m_bank, S_0, D=4)
    print("The first calculation of Q is %0.4f m3/s" % (Qo))
    h_n = interpolate_h(Q, b, m_bank, S_0,D=4)
    print("The calculated h after interpolation is %0.4f m and the error is %0.4f " % (h_n[0],h_n[1]))
    myq = calc_discharge(b, h_n[0], m_bank, S_0, D=4)
    print("The calculated flux with the interpolated h is %0.4f m3/s" % (myq))
