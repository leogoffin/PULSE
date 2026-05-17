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


def PIDcontrol(ref, measure, last_errors=None, Kp=1.0, Ki=0, Kd=0, tol=1e-5):
    """
    PID controller for iterative convergence.

    Args:
        ref (float): Reference (target) value.
        measure (float): Measured value.
        last_errors (list, optional): Stored error history. Defaults to None.
        Kp (float, optional): Proportional gain. Defaults to 1.0.
        Ki (float, optional): Integral gain. Defaults to 0.
        Kd (float, optional): Derivative gain. Defaults to 0.
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
    I = Ki * sum(last_errors)

    if len(last_errors) > 1:
        D = Kd * (last_errors[-1] - last_errors[-2])
    else:
        D = 0.0

    output = P + I + D

    if abs(error) < tol:
        print("\n")
        print(f"PID converged in {len(last_errors)} steps")
        print(f"Measured: {measure:.6f}\n")
        return True, 0.0, last_errors

    return False, output, last_errors


def secant(xn,xn_1,fxn,fxn_1,Min= None,Max= None):
    new_xn = xn - fxn * (xn-xn_1)/(fxn-fxn_1)
    if Min is not None :
        new_xn = max(Min,new_xn)
    if Max is not None : 
        new_xn = min(Max,new_xn)
    return new_xn,xn,0,fxn_1