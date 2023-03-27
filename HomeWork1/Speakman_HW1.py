from math import exp
from math import sqrt
from math import log
from math import cos as cosine
import cmath
import numpy as np
import sympy
import plotly as plt
import matplotlib.pyplot as plt



def quad_formula(a: float, b: float, c: float, tolerance: float = 10**10) -> tuple:
    '''
    Quadratic formula to avoid large round off errors. Uses cmath to handle imaginary numbers
    '''
    q1: float = -b
    q2: float = b * b
    q3: float = 4.*a*c
    q4: float = 2.*a

    # print(f"numerator plus = {q1 + sqrt(q2 - q3)}")
    # print(f"numerator minus = {q1 - sqrt(q2 - q3)}")

    x1, x2 = None, None 

    if (q2 - q3) >= 0.:
        # real solutions
        if (q1 + sqrt(q2 - q3)) <= tolerance:
            # values are too similar
            x1 = (2.*c) / (q1 - sqrt(q2 - q3))
        else:
            x1 = (q1 + sqrt(q1 - q3)) / q4

        if (q1 - sqrt(q2 - q3)) <= tolerance:
            x2 = (2.*c) / (q1 + sqrt(q2 - q3))
        else:
            x2 = (q1 - sqrt(q2 - q3)) / q4

    else:
        # imaginary solutions
        q23 = cmath.sqrt(q2 - q3) #"Factoring" out the imaginary component
        if (q1 + q23) <= tolerance:
            # values are too similar
            x1 = (2.*c) / (q1 - q23)
        else:
            x1 = (q1 + q23) / q4

        if (q1 - q23) <= tolerance:
            x2 = (2.*c) / (q1 + q23)
        else:
            x2 = (q1 - q23) / q4

    return x1, x2


def bisectional_method(f, lBoundary: int or float, rBoundary: int or float, tolerance: float = 10**-10, iterations = 10**9, 
                            *args, **kwargs) -> float:
    '''
    Takes a function and finds the roots. Args and Kwargs are used for additional inputs to the function
    '''

    fLeft = f(lBoundary, *args, **kwargs)
    fRight = f(rBoundary, *args, **kwargs)

    if (fLeft * fRight) > 0.:
        print("Unable to find roots: points must be on opposite site of y axis")
        print(f"f({lBoundary}) = {fLeft}, f({rBoundary}) = {fRight}")
        return None


    for i in range(iterations):
        midpoint = (lBoundary + rBoundary) / 2.

        fMidpoint = f(midpoint, *args, **kwargs)

        if (fMidpoint * fRight) < 0.:
            lBoundary, fLeft = midpoint, fMidpoint

        else:
            rBoundary, fRight = midpoint, fMidpoint

        if (rBoundary - lBoundary) < tolerance:
            # print(i)
            # print(rBoundary - lBoundary)
            return midpoint

    else:
        print(f"Failed to find root after {iterations} attempts")



def relaxation_method(f, guess: int or float, tolerance: float = 10**-10, max_iterations = 200, 
                      *args, **kwargs) -> float:
    '''
    Relaxation method for finding roots.

    Takes a function as an input (f), an inital guess (int or float), and will return a float for the root.

    Will, by default, only do 200 iterations, but can be adjusted. Tolerance for success if 10**-10 by default.
    '''
    xprev = xnext = guess

    for _ in range(max_iterations):
        xprev = xnext
        xnext = f(xprev, *args, **kwargs)

        if (abs(xnext - xprev) < tolerance):
            return xnext

    else:
        print(f"Failed to converge after {max_iterations} iterations")
        print(f"Values are:\n\tInitial guess = {guess}\n\tXprevious = {xprev}\n\tXnext = {xnext}\n\tError = {abs(xnext - xprev)}")

        return xnext




def Lagrangian_Basis(point: float, jndex: int, x_points: list or np.array) -> float:
    '''
    '''

    base = 1.

    for kndex in range(len(x_points)):
        if kndex != jndex:
            base *= (point - x_points[kndex]) / (x_points[jndex] - x_points[kndex])

    return base


def Lagrangian_Fit(data_points: list or np.array or dict or tuple, points = 100) -> dict:
    '''
    data points can be a numpy array, dictionary, list, or tuple. If anything other then a dictionary the x points are
    interpreted to be 0, 1, 2, 3, etc. while the input list is taken as the correspoding y points.
    '''
    if isinstance(data_points, dict):
        x = np.array([x for x in data_points.keys()])
        y = np.array([y for y in data_points.values()])
    else:
        x = np.array([x for x in range(len(data_points))])
        if isinstance(data_points, (list, tuple)):
            y = np.array(data_points)

    x_points = np.linspace(x[0], x[len(x) - 1], points)

    fit = {}

    for point in x_points:

        interpolation = 0.
        for j, weight in enumerate(y):
            interpolation += weight * Lagrangian_Basis(point, j, x)
            
        fit[point] = interpolation    

    return fit


def Lagrangian_Coefficents(data_points: list or np.array or dict or tuple) -> sympy.core.add:
    '''
    Uses sympy to find the 
    '''

    if isinstance(data_points, dict):
        x_points = np.array([x for x in data_points.keys()])
        y_points = np.array([y for y in data_points.values()])
    else:
        x_points = np.array([x for x in range(len(data_points))])
        if isinstance(data_points, (list, tuple)):
            y_points = np.array(data_points)

    x = sympy.Symbol("x")

    poly_fit = 0
    for j, weight in enumerate(y_points):

        p = 1.
        for k, x_point in enumerate(x_points):
            if j != k:
                p *= (x - x_points[k]) / (x_points[j] - x_points[k])

        poly_fit += sympy.simplify(weight*p)

    return poly_fit


def problem1():
    '''
    Find the roots of ax**2 + bx + c where 
    
    a = 1.6x10**-4, b = 1.5x10**4, c = 1.3x10**-4

    &

    a = 1.6x10**-4, b = -1.4x10**4, c = 1.3x10**-4
    '''
    # a, b, c = 1.6*10**-4, 1.5*10**4, 1.3*10**-4
    a, b, c = 10**-4, 10**4, 10**-4
    print(f"a = {a}, b = {b}, c = {c}")
    roots: tuple = quad_formula(a, b, c)
    print(f"x1 = {roots[0]}, x2 = {roots[1]}")

    # a, b, c = 1.6*10**-4, -1.4*10**4, 1.3*10**-4
    a, b, c = 10**-4, -10**4, 10**-4
    print(f"a = {a}, b = {b}, c = {c}")
    roots: tuple = quad_formula(a, b, c)
    print(f"x1 = {roots[0]}, x2 = {roots[1]}")


def problem2():
    '''
    Find the coefficents of p(x) = a0 + a1 x + a2 x**x + a3 x**3

    that interpolates cos(x)
    '''

    f = lambda x: cosine(x)

    cos = {0:  1.,
           1:  0.5403023058681398,
           2: -0.4161468365471424,
           3: -0.9899924966004454}

    # x = sympy.Symbol("x")

    poly_fit = Lagrangian_Coefficents(cos)  # finds an algebraic expression for the polynomial

    poly_fit = sympy.lambdify(sympy.Symbol("x"), poly_fit)

    x_expected = np.linspace(0, 5, 100)
    expected = {x: f(x) for x in x_expected}

    interportional = {x:poly_fit(x) for x in x_expected}

    # plt.title("Interpolation of Cos(x) from 0 to 3")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.plot(expected.keys(), expected.values(), label = "Expected")
    # plt.plot(interportional.keys(), interportional.values(), "--", label = "Interpolation")
    # plt.legend()
    # plt.show()

    plt.title("Interpolation of Cos(x) from 0 to 5")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(expected.keys(), expected.values(), label = "Expected")
    plt.plot(interportional.keys(), interportional.values(), "--", label = "Interpolation")
    plt.legend()
    plt.show()


def problem3():
    '''
    '''
    f = lambda x: -(x/2) + exp(1 - x**2)
    fln = lambda x: sqrt(-log(x) + log(2) + 1)

    root_bi = bisectional_method(f, 1, 2)
    print(root_bi)
    root_re = relaxation_method(fln, 1.0)
    print(root_re)




if __name__ in "__main__":
    problem1()
    problem2()
    problem3()