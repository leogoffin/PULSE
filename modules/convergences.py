def Pcontrol(ref, measure, last_errors=None, Kp=1.0, tol=1e-5):
    """
    Proportional controller for iterative convergence.

    Args:
        ref (float): Reference (target) value.
        measure (float): Measured value.
        last_errors (list, optional): Stored error history. Defaults to None.
        Kp (float, optional): Proportional gain. Defaults to 1.0.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (converged flag, control output, error history)
    """

    if last_errors is None:
        last_errors = []

    if ref == 0:
        error = measure
    else:
        error = (ref - measure) / ref

    last_errors.append(error)

    P = Kp * error
    output = P

    if abs(error) < tol:
        print("\n")
        print(f"P converged in {len(last_errors)} steps")
        print(f"Measured: {measure:.6f}\n")
        return True, 0.0, last_errors

    return False, output, last_errors

def secant(xn,xn_1,fxn,fxn_1,Min= None,Max= None):
    """Compute secant."""
    new_xn = xn - fxn * (xn-xn_1)/(fxn-fxn_1)
    if Min is not None :
        new_xn = max(Min,new_xn)
    if Max is not None : 
        new_xn = min(Max,new_xn)
    return new_xn,xn,0,fxn_1

def bisection(f, a, b, tol=1e-6, max_iter=100):
    """
    Bisection method for solving f(x) = 0.

    Parameters
    ----------
    f : function
        Function for which we seek a root.
    a, b : float
        Interval bounds such that f(a)*f(b) < 0.
    tol : float
        Convergence tolerance.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    root : float
        Estimated root.
    iterations : int
        Number of iterations performed.
    """

    fa = f(a)
    fb = f(b)

    if fa * fb > 0:
        raise ValueError("Root is not bracketed: f(a) and f(b) must have opposite signs.")

    for i in range(max_iter):

        c = 0.5 * (a + b)
        fc = f(c)

        if abs(fc) < tol or abs(b - a) < tol:
            return c, i + 1

        if fa * fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc

    raise RuntimeError("Maximum number of iterations reached.")