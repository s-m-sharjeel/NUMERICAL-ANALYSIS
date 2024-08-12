import math
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

# returns the function value at x
def f(func, x):
    return eval(func, {"__builtins__": None},
                {"x": x, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e,
                 "log": math.log})

# evaluates the first derivative at x
def fprime(func, x1):
    x = sp.Symbol('x')
    function = sp.sympify(func)
    expression = str(function.diff(x).subs(x, x1))
    return eval(expression, {"__builtins__": None}, {"x": x1, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e, "log": math.log})

def spline(x, y, clamped, boundary):

    n = len(x)

    S = []

    # initializing the cofficicents matrix
    A = np.zeros((n,n))

    # initializing the vetor h
    h = np.zeros((n,1))
    # initializing the constant vector b in Ax=b
    Cnstb=np.zeros((n,1))
    # initializing other coefficients b's and d's of the polynomials
    d = np.zeros((n,1))
    b = np.zeros((n,1))

    # defining height h
    for i in range(0,n-1):
        h[i] = x[i+1] - x[i]

    if clamped:

        # defining A and b for clamped boundary condition
        A[0,0] = 2 * h[0]
        A[0,1] = h[0]
        A[n-1,n-1] = 2 * h[n - 2]
        A[n-1,n-2] = h[n - 2]
        Cnstb[0] = 3 * (y[1] - y[0]) / h[0] - 3 * boundary[0]
        Cnstb[n-1] = 3 * boundary[1] - 3 * (y[n - 1] - y[n - 2]) / h[n - 2]

    else:

        # defining A and b for natural boundary condition 
        A[0,0] = 1
        A[n-1,n-1] = 1
        Cnstb[0] = 0
        Cnstb[n-1] = 0

    for i in range(1,n-1):
        A[i,i-1] = h[i-1] 
        A[i,i] = 2 * ((h[i-1] + h[i]))
        A[i,i+1] = h[i]
        Cnstb[i] = 3 * (y[i+1] - y[i]) / h[i] - 3 * (y[i] - y[i-1]) / h[i-1]  

    # solving A = xb
    c = np.linalg.solve(A,Cnstb)

    # calculate other two coefficients b's and d's of the polynomials
    for i in range(0,n-1):
        b[i] = (y[i+1] - y[i]) / h[i] - h[i] * (c[i+1] + 2 * c[i]) / 3
        d[i] = (c[i+1] - c[i]) / (3 * h[i])

    # append the partial spline polynomial to the array containing the collection 
    for i in range(0,n-1):
        s = f"{y[i]} + ({b[i]} * (x - {x[i]})) + ({c[i]} * (x - {x[i]}) ** 2) + ({d[i]} * (x - {x[i]}) ** 3)"
        S.append(s)

    return S

def equallySpacedVal(func, a, b, n):

    x = np.linspace(a, b, n)
    y = []
    for data in x:
        y.append(f(func, data))

    return x, y

func = "(e ** x) + (2 ** -x) + (2 * cos(x)) - 6"
a = -2
b = 2

# changing value of n here will get you a spline constructed from (n + 1) equally distanced data points
# and consequently n polynomials for n intervals (generic for any value of n)
n = 10

# since, n intervals require n + 1 data points
x, y = equallySpacedVal(func, a, b, n + 1)

boundary = [fprime(func, a), fprime(func, b)]

# can be made true for testing, which gives a better approximation (less error as can be seen in the error function E(x))
clamped = False

# constructs the spline
S = spline(x, y, clamped, boundary)

# plotting the data points
plt.plot(x, y, 'x', label="Data-Points")

# plotting the function for comparison using 100 points and a line which can be observed to be overlapped by the spline if n is sufficiently large
funcX, funcY = equallySpacedVal(func, a, b, 100)
plt.plot(funcX, funcY, label="f(x)")

errorValues = []

# plotting each polynomial si in its interval by taking n points (precision)
precision = 20
for i in range(0,n):
    var = np.linspace(x[i], x[i+1], precision)
    s = f(S[i], var)
    funcValues = []
    j = 0
    for element in var:
        funcVal = f(func, element)
        funcValues.append(funcVal)
        errorVal = s[j] - funcVal
        errorValues.append(errorVal)
        j = j + 1
        
    plt.plot(var,s, label=f"Polynomial {i + 1}")

# comparing the spline and the function itself by creating an error function
var1 = np.linspace(a,b,precision*n)
plt.plot(var1, errorValues, label="E(x)")
    
plt.title(f"Cubic Spline (Clamped = {clamped})")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()
plt.show()

# Article 3.5, question 32:

clamped = True

x1 = [1, 2, 5, 6, 7, 8, 10, 13, 17]
y1 = [3.0, 3.7, 3.9, 4.2, 5.7, 6.6, 7.1, 6.7, 4.5]
boundary1 = [1.0, -0.67]

S1 = spline(x1, y1, clamped, boundary1)

for i in range(0,len(x1) - 1):
    var=np.linspace(x1[i],x1[i+1],20)
    s = f(S1[i], var)
    plt.plot(var,s)
    
plt.plot(x1, y1, 'x', label="Data-Points for Spline 1")
plt.plot([x1[0], x1[len(x1) - 1]], [0, 0], label = "Curve 1 ^")

x2 = [17, 20, 23, 24, 25, 27, 27.7]
y2 = [4.5, 7.0, 6.1, 5.6, 5.8, 5.2, 4.1]
boundary2 = [3.0, -4.0]

S2 = spline(x2, y2, clamped, boundary2)

for i in range(0,len(x2) - 1):
    var=np.linspace(x2[i],x2[i+1],20)
    s = f(S2[i], var)
    plt.plot(var,s)
    
plt.plot(x2, y2, 'x', label="Data-Points for Spline 2")
plt.plot([x2[0], x2[len(x2) - 1]], [0, 0], label = "Curve 2 ^")

x3 = [27.7, 28, 29, 30]
y3 = [4.1, 4.3, 4.1, 3.0]
boundary3 = [0.33, -1.5]

S3 = spline(x3, y3, clamped, boundary3)

for i in range(0,len(x3) - 1):
    var=np.linspace(x3[i],x3[i+1],20)
    s = f(S3[i], var)
    plt.plot(var,s)
    
plt.plot(x3, y3, 'x', label="Data-Points for Spline 3")
plt.plot([x3[0], x3[len(x3) - 1]], [0, 0], label = "Curve 3 ^")

plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.legend()
plt.show()

print("The splines, and their polynomials, constructed are shown below:")

print("Spline 1:")
for i in range(0, len(x1) - 1):
    print(f"{x1[i]} - {x1[i + 1]} : {S1[i]}")

print("Spline 2:")
for i in range(0, len(x2) - 1):
    print(f"{x2[i]} - {x2[i + 1]} : {S2[i]}")

print("Spline 3:")
for i in range(0, len(x3) - 1):
    print(f"{x3[i]} - {x3[i + 1]} : {S3[i]}")