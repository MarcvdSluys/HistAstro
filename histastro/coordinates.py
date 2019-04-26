def obliquity(JD):
    """Compute the obliquity of the ecliptic in radians from the JD(E)  (Seidelman 1992, Eq. 3.222-1)"""
    tJC = (JD - 2451545.0)/36525
    #eps = 23.4392911*d2r
    #eps += (-46.815*tJC - 0.00059*tJC**2 + 0.001813*tJC**3)*as2r
    eps = 0.409092804 - 2.269655e-4*tJC - 2.86e-9*tJC**2 + 8.78967e-9*tJC**3  # Checked
    return eps

