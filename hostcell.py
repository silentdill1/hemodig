# TODO: improve polynomial fit


def get_hb_concentration(t):  # returns host cell hb concentration in mM
    if t <= 32:
        '''
        a = [7.03233, 0.326897, -0.495141, 0.278402, -0.0803537, 0.0137342, -0.00150074, 0.000109446,
             -5.44031 * 10 ** (-6), 1.84645 * 10 ** (-7), -4.18157 * 10 ** (-9), 5.87513 * 10 ** (-11),
             -3.95697 * 10 ** (-13), -1.00864 * 10 ** (-15), 3.55859 * 10 ** (-17), -1.74286 * 10 ** (-19)]
             funktioniert nicht ??
             '''
        a = [7.06231, -0.362935, 0.4128, -0.174946, 0.0377169, -0.00473409, 0.0003715, -0.0000189184,
             6.33731 * 10 ** -7, -1.38266 * 10 ** -8, 1.88175 * 10 ** -10, -1.4294 * 10 ** -12, 4.27014 * 10 ** -15,
             4.12342 * 10 ** -18]
        # polynomial coefficients
        hbconcentration = 0
        for i in range(0, len(a)):
            hbconcentration += a[i]*t**i
    else:
        hbconcentration = 3.6
    return hbconcentration
