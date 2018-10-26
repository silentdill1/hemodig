# TODO: improve polynomial fit


def get_hb_concentration(t):  # returns host cell hb concentration in mM
    if t <= 32:
        a = [7.06231, -0.362935, 0.4128, -0.174946, 0.0377169, -0.00473409, 0.0003715, -0.0000189184,
             6.33731 * 10 ** -7, -1.38266 * 10 ** -8, 1.88175 * 10 ** -10, -1.4294 * 10 ** -12, 4.27014 * 10 ** -15,
             4.12342 * 10 ** -18]
        # polynomial coefficients
        hb_concentration = 0
        for i in range(0, len(a)):
            hb_concentration += a[i]*t**i
    else:
        hb_concentration = 3.6
    return hb_concentration
