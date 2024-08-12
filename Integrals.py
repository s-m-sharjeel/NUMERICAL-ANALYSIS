import math
import numpy as np
import sympy as sp
import scipy.optimize as optimize

# returns the function value at x
def f(func, x):
    return eval(func, {"__builtins__": None},
                {"x": x, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e,
                 "log": math.log})

# returns the nth derivative
def fprime(func, n):
    x = sp.Symbol('x')
    function = sp.sympify(func)
    expression = str(function.diff(x, n))
    return expression


def integrate(func, a, b):
    x = sp.symbols("x")
    function = sp.sympify(func)
    return sp.integrate(function, (x, a, b))

def Trapezoidal(func, a, b, n):

    x = np.linspace(a, b, n + 1)
    y = []
    for element in x:
        y.append(f(func, element))

    integral = trapezoidal_general(x, y)

    h = (x[n] - x[0]) / n

    # error = 0
    second_derivative = fprime(func, 2)

    # direct formula:
    # error = -((b - a)/12)*(h ** 2)*maxPoint(second_derivative, a, b)[1]

    error = 0

    for i in range(0, n):
        error_term = -((h ** 3)/12) * maxPoint(second_derivative, x[i], x[i + 1])[1]
        error += error_term

    return integral, error

def trapezoidal_general(x, y):

    n = len(x) - 1
    h = (x[n] - x[0]) / n

    sum = (y[n] + y[0]) / 2

    for i in range(1, n):
        sum += y[i]

    return sum * h

def Simpson(func, a, b, n):

    x = np.linspace(a, b, n + 1)
    y = []
    for element in x:
        y.append(f(func, element))

    integral = simpson_general(x, y)
    
    h = (b - a) / n

    fourth_derivative = fprime(func, 4)

    # direct formula:
    # error = -((b - a) / 80) * (h ** 4) * maxPoint(fourth_derivative, a, b)[1]

    error = 0

    for i in range(0, n):
        if (i % 3 == 0):
            error_term = -(3*(h ** 5)/80) * maxPoint(fourth_derivative, x[i], x[i + 3])[1]
            error += error_term

    return integral, error

def simpson_general(x, y):

    n = len(x) - 1

    if n % 3 != 0:
        print("n must be a multiple of 3!")
        return

    h = (x[n] - x[0]) / n

    sum = y[n] + y[0]

    for i in range(1, n):
        if i % 3 == 0:
            sum += 2 * y[i]
        else:
            sum += 3 * y[i]

    return 3 * sum * h / 8

def maxPoint(func, a, b):

    interval = [a, b]
    interval = tuple(interval)
    negative_function = lambda x: -f(func, x)
    minimum = optimize.fmin(negative_function, interval[0], disp=False)

    maximum_point = minimum[0]
    maximum_value = f(func, maximum_point)
    
    return maximum_point, maximum_value

# 1:

func = "sin(x ** 2)"
n = 30
a = 0
b = 3

trapezoidal = Trapezoidal(func, a, b, n)
simpson = Simpson(func, a, b, n)

s = f"f(x) = {func}\na = {a}\nb = {b}\nn = {n}\n\n"
s += "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|"  + '\n'
s += "|" + f"{trapezoidal[0]}".center(25) + "|" + f"{simpson[0]}".center(25) + "|" + '\n'
print(s)

# 2:

s = f"f(x) = {func}\na = {a}\nb = {b}\n\n"
s += "|" + "n".center(10) + "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|" + '\n'

i = 1
while(i < 6):
    trapezoidal = Trapezoidal(func, a, b, i * 3)
    simpson = Simpson(func, a, b, i * 3)
    s += "|" + f"{i * 3}".center(10) + "|" + f"{trapezoidal[0]}".center(25) + "|" + f"{simpson[0]}".center(25) + "|" + '\n'
    i = i + 1

print(s)

# 3: done in pdf

# 4:

func = "(1 + x)**(1/2)"
a = 0
b = 0.1
n = 12

trapezoidal = Trapezoidal(func, a, b, n)
simpson = Simpson(func, a, b, n)

s = f"f(x) = {func}\na = {a}\nb = {b}\nn = {n}\n\n"
s += "|" + "Actual".center(25) + "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|"  + '\n'
s += "|" + f"{integrate(func, a, b)}".center(25) + "|" + f"{trapezoidal[0]}".center(25) + "|" + f"{simpson[0]}".center(25) + "|" + '\n'
print(s)

# 5:

s = "True Errors:\n\n" + "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|" + '\n'
s += "|" + f"{abs(trapezoidal[0] - integrate(func, a, b))}".center(25) + "|" + f"{abs(simpson[0] - integrate(func, a, b))}".center(25) + "|" + '\n'
print(s)

s = "Approximate Errors:\n\n" + "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|" + '\n'
s += "|" + f"{(trapezoidal[1])}".center(25) + "|" + f"{(simpson[1])}".center(25) + "|" + '\n'
print(s)

# 6:

# 15 data points
x = np.linspace(0, 84, 15)
y = [124, 134, 148, 156, 147, 133, 121, 109, 99, 85, 78, 89, 104, 116, 123]

# 14 intervals
n = 14
a = 0
b = 84

# for Simpson's, n must be a multiple of 3
print("Simpson's Method cannot be used")

trapezoidal = trapezoidal_general(x, y)
simpson = simpson_general(x, y)

s = f"x = {x}\ny = {y}\n"
s += f"a = {a}\nb = {b}\nn = {n}\n\n"
s += "|" + "Trapezoidal".center(25) + "|" + "Simpson's".center(25) + "|"  + '\n'
s += "|" + f"{trapezoidal}".center(25) + "|" + f"{simpson}".center(25) + "|" + '\n'
print(s)
