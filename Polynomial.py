import math
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

# returns the function value at x
def f(func, x):
    return eval(func, {"__builtins__": None},
                {"x": x, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e,
                 "log": math.log})

def Lagrange(x, y):

    n = len(x)

    # computing each term of the lagrange polynomial separately
    L = []
    for i in range(0, n):
        term = ""
        numerator = "1"
        denominator = 1
        for j in range(0, n):
            if j != i:
                numerator += f"*(x - ({x[j]}))"
                denominator *= (x[i] - x[j])
        term = f"(({numerator}) / ({denominator}))"
        L.append(term)

    # concatenating all the terms of the lagrange polynomial previously computed
    polynomial = "0"
    for i in range(0, n):
        polynomial += f" + (({L[i]}) * ({y[i]}))"

    return polynomial

def DDTable(x, y):
    
    # creating an array for the divided differences
    DD = []

    n = len(x)

    index = range(0, n)

    DD.append(index)
    DD.append(x)
    DD.append(y)

    # computing the divided differences
    for j in range (2, n + 1):
        nthDD = []
        for i in range(0, (n - j) + 1):
            val = (DD[j][i + 1] - DD[j][i]) / (DD[1][i + (j - 1)] - DD[1][i])
            nthDD.append(val)
            
        DD.append(nthDD)

    # arranging the data so that it can be printed in a triangular form:

    # taking transpose for tabular array and movine each element in the i th row to (i * 2)th row to form gaps between the data
    arr = [[None for _ in range(len(DD))] for _ in range(len(DD[0]) * 2)]
    for i in range(0, len(DD)):
        for j in range(0, len(DD[i])):
            arr[j * 2][i] = DD[i][j]

    tableArr = arr

    # moving all the elements in the 3rd+ columns (n - 2) places down to form a triangular arrangement
    arr = [[None for _ in range(len(tableArr[0]))] for _ in range(len(tableArr)*2)]
    for i in range(0, len(tableArr)):
        for j in range(0, len(tableArr[i])):
            if j > 2:
                arr[i + (j - 2)][j] = tableArr[i][j]
            else: arr[i][j] = tableArr[i][j]

    tableArr = arr

    # writing a string for the table
    table = ""

    # writing the header for the table
    for k in range(1, (len(tableArr[0]) * 26) + 2):
        table += "_"
    
    table += "\n|" + "i".center(25) + "|" + "xi".center(25) + "|"
    for i in range(0, n):
        val = "f[xi"
        for j in range(0, i):
            if (j < 1):
                val += f", xi+{j + 1}"
            else:
                val += f",... xi+{i}"
                break
        val += "]"
        table += val.center(25) + "|"

    table += "\n|"
    for k in range(1, (len(tableArr[0]) * 26) + 1):
        if k % 26 == 0:
            table += "|"
        else: table += "_"

    table += "\n"

    for i in range(0, (n * 4) - 2):
        table += "|"
        if (i % 2 != 0):
            for k in range(1, (len(tableArr[0]) * 26) + 1):
                if k % 26 == 0:
                    table += "|"
                else: table += "_"
        else: 
            for j in range(0, len(tableArr[(int)(i/2)])):
                if (tableArr[(int)(i/2)][j] != None):
                    table += f"{tableArr[(int)(i/2)][j]}".center(25)
                else: table += "".center(25)
                table += "|"
        table += "\n"

    # printing the table on console
    # print(table)

    # writing the table in a text file

    filepath = "DDTable.txt"

    with open(filepath, "w") as file:
        file.write(table)

    return DD


def DDPolynomial(x, y):

    DD = DDTable(x, y)

    DDTerms = []

    for i in range(2, len(DD)):
        term = f"({DD[i][0]})"
        for j in range(0, i - 2):
            term += f"*(x - {DD[1][j]})"
        DDTerms.append(term)

    polynomial = "0"

    for term in DDTerms:
        polynomial += f" + ({term})"

    filepath = "DDTable.txt"

    with open(filepath, "a") as file:
        file.write("\nDDPolynomial:\n" + polynomial)

    return polynomial

def equallySpacedVal(func, a, b, n):

    x = np.linspace(a, b, n)
    y = []
    for data in x:
        y.append(f(func, data))

    return x, y

def maxPoint(func, a, b):

    # making an array for all the critical points and the end-points
    points = [a, b]

    # getting all the critical points (real + complex roots for fprime)
    x = sp.Symbol('x')
    function = str(sp.sympify(func).diff(x))
    criticalPoints = sp.solve(function)

    # adding all the real roots to the points at which function value is to be checked for max point
    for element in criticalPoints:
        if (type(element) == sp.core.numbers.Float):
            points.append(element)

    # finding the point at which the function value is max
    maximum = points[0]
    maxFuncVal = f(func, points[0])
    for element in points:
        funcVal = f(func, element)
        if (funcVal > maxFuncVal):
            maxFuncVal = funcVal
            maximum = element

    return maximum

#_____________________________________________________________________________________________________________________________________#

# 1:

func = "(e ** x) + (2 ** -x) + (2 * cos(x)) - 6"
a = -2
b = 2

# value of n can be changed here (code is generic for any value of n)
n = 10

x, y = equallySpacedVal(func, a, b, n)

print(f"For n = {n}:")
print("DDTable: \n")
DDTable(x, y)
DDPoly = DDPolynomial(x, y)
print("DD Polynomial: \n" + DDPoly)
lagrangePoly = Lagrange(x, y)
print("Lagrange Polynomial: \n" + lagrangePoly)

#_____________________________________________________________________________________________________________________________________#

# 2 & 3:

func = "(e ** x) + (2 ** -x) + (2 * cos(x)) - 6"
a = -2
b = 2

# Note:
# the divided difference polynomial can become very unstable for large values of n. 
# this is because the divided differences of higher order can become very large, even for small changes in the data points. 
# this can lead to inaccurate approximations of the function, especially at the edges of the data range.

# we will get inaccurate approximations of the function (near the endpoints of the interval) for:
n = 100

# this, however, won't be the case for:
n = 50

x, y = equallySpacedVal(func, a, b, n)

print(f"For n = {n}:")
print("DDTable: \n(written on DDTable.txt)")
DDTable(x, y)
DDPoly = DDPolynomial(x, y)
print("DD Polynomial: \n" + DDPoly)
lagrangePoly = Lagrange(x, y)
print("Lagrange Polynomial: \n" + lagrangePoly)

# PLOTTING GRAPH:

precision = 500

# create an array of 500 x-values (precision) between a and b
xValues = np.linspace(a, b, precision)
lagrangeValues = []
DDValues = []
funcValues = []
lagrangeErrorValues = []
DDErrorValues = []

# creates an array for the function, lagrange and error values at those x-values
for element in xValues:
    lagrangeVal = f(lagrangePoly, element)
    lagrangeValues.append(lagrangeVal)
    DDVal = f(DDPoly, element)
    DDValues.append(DDVal)
    funcVal = f(func, element)
    funcValues.append(funcVal)
    lagrangeError = funcVal - lagrangeVal
    lagrangeErrorValues.append(lagrangeError)
    DDError = funcVal - DDVal
    DDErrorValues.append(DDError)

# plots and labels the Lagrange and DD polynomials and their error functions together on the same graph:
plt.plot(x, y, 'x', label="Data-Points")
plt.plot(xValues, funcValues, label="f(x)")
plt.plot(xValues, lagrangeValues, label="Pn(x) for Lagrange")
plt.plot(xValues, DDValues, label="Pn(x) for DD")
plt.plot(xValues, lagrangeErrorValues, label="E(x) for Lagrange")
plt.plot(xValues, DDErrorValues, label="E(x) for DD")

plt.legend()
plt.title(f"Lagrange and DD Polynomials (for n = {n})")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()

# plots the error functions for the Lagrange and DD Polynomials separately (has already been plotted with (2) above):
plt.plot(xValues, lagrangeErrorValues, label="E(x) for Lagrange")
plt.plot(xValues, DDErrorValues, label="E(x) for DD")

plt.legend()
plt.title(f"Error Functions (for n = {n})")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.grid(True)
plt.show()

#_____________________________________________________________________________________________________________________________________#

# 4: (Ex 3.1 Q 19)

x = [0, 6, 10, 13, 17, 20, 28]
y1 = [6.67, 17.33, 42.67, 37.33, 30.10, 29.31, 28.74]
y2 = [6.67, 16.11, 18.89, 15.00, 10.56, 9.44, 8.89]

# b) Finding an approximate maximum average weight for each sample by determining the maximum of the interpolating polynomial:

# For Sample # 1: (y1)

lagrangePoly1 = Lagrange(x, y1)
lagrangeMax1 = maxPoint(lagrangePoly1, 0, 28)

print(f"Maximum average weight for Sample 1: {f(lagrangePoly1, lagrangeMax1)} mg ({lagrangeMax1} days)")

# For Sample # 2: (y2)

lagrangePoly2 = Lagrange(x, y2)
lagrangeMax2 = maxPoint(lagrangePoly2, 0, 28)

print(f"Maximum average weight for Sample 2: {f(lagrangePoly2, lagrangeMax2)} mg ({lagrangeMax2} days)")

# a) Using Lagrange interpolation to approximate the average weight curve for each sample:

a = 0
b = 28
precision = 100

# create an array of 100 x-values (precision) between a and b
xValues = np.linspace(a, b, precision)
lagrangeValues1 = []
lagrangeValues2 = []

# creates an array for the function, lagrange and error values at those x-values
for element in xValues:
    lagrangeVal1 = f(lagrangePoly1, element)
    lagrangeValues1.append(lagrangeVal1)
    lagrangeVal2 = f(lagrangePoly2, element)
    lagrangeValues2.append(lagrangeVal2)

# Graphs plotted separately for each sample:

# for Sample 1:

# plots and labels the data-points and lagrange polynomial values 
plt.plot(x, y1, 'x', label="Data-Points")
plt.plot(xValues, lagrangeValues1, label="Lagrange")

# plotting the max point for the lagrange polynomials
plt.plot([a, b], [f(lagrangePoly1, lagrangeMax1), f(lagrangePoly1, lagrangeMax1)], '--', label="Max-Height")
plt.plot(lagrangeMax1, f(lagrangePoly1, lagrangeMax1), 'o', label="Max-Point")

plt.legend()
plt.title("Average weight curve for Sample 1")
plt.xlabel("Days")
plt.ylabel("Average Weight (mg)")
plt.grid(True)
plt.show()

# for Sample 2:

# plots and labels the data-points and lagrange polynomial values 
plt.plot(x, y2, 'x', label="Data-Points")
plt.plot(xValues, lagrangeValues2, label="Lagrange")

# plotting the max point for the lagrange polynomials
plt.plot([a, b], [f(lagrangePoly2, lagrangeMax2), f(lagrangePoly2, lagrangeMax2)], '--', label="Max-Height")
plt.plot(lagrangeMax2, f(lagrangePoly2, lagrangeMax2), 'o', label="Max-Point")

plt.legend()
plt.title("Average weight curve for Sample 2")
plt.xlabel("Days")
plt.ylabel("Average Weight (mg)")
plt.grid(True)
plt.show()

# plotted on the same graph:

# plots and labels the data-points and lagrange polynomial values 
plt.plot(x, y1, 'x', label="Data-Points for Sample 1")
plt.plot(x, y2, 'x', label="Data-Points for Sample 2")
plt.plot(xValues, lagrangeValues1, label="Lagrange for Sample 1")
plt.plot(xValues, lagrangeValues2, label="Lagrange for Sample 2")

# # plotting the max point for the lagrange polynomials
plt.plot([a, b], [f(lagrangePoly1, lagrangeMax1), f(lagrangePoly1, lagrangeMax1)], '--', label="Max-Height for Sample 1")
plt.plot([a, b], [f(lagrangePoly2, lagrangeMax2), f(lagrangePoly2, lagrangeMax2)], '--', label="Max-Height for Sample 2")
plt.plot(lagrangeMax1, f(lagrangePoly1, lagrangeMax1), 'o', label="Max-Point for Sample 1")
plt.plot(lagrangeMax2, f(lagrangePoly2, lagrangeMax2), 'o', label="Max-Point for Sample 2")

plt.legend()
plt.title("Average weight curve for each Sample")
plt.xlabel("Days")
plt.ylabel("Average Weight (mg)")
plt.grid(True)
plt.show()

#_____________________________________________________________________________________________________________________________________#