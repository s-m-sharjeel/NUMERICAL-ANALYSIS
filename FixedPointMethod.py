import sympy as sp
import math
import matplotlib.pyplot as plt
import numpy as np
import time

def FixedPointMethod(func, p0, tolerance):

    # returns the function value at x
    def f(x1):
        return eval(func, {"__builtins__": None}, {"x": x1, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e, "log": math.log})

    # evaluates the first derivative at x
    def fprime(x1):
        x = sp.Symbol('x')
        function = sp.sympify(func)
        expression = str(function.diff(x).subs(x, x1))
        return eval(expression, {"__builtins__": None}, {"x": x1, "sin": math.sin, "cos": math.cos, "tan": math.tan, "pi": math.pi, "e": math.e, "log": math.log})

    

    # initializing variables
    iterations = 0
    p = p0
    root_approx = []
    root_approx.append(p)

    start = time.process_time()

    # checking if root is at the initial guess
    if (abs(f(p0)) < tolerance):
        # root found at initial guess
        root = p0

    else:
        # iterate until convergence is achieved
        while True:

            iterations += 1

            # new root approximation
            p_new = f(p)
            root_approx.append(p_new)

            # bound on function value
            function_error = abs(f(p_new) - (p_new))

            # calculating relative error if pn != 0
            if (p != 0):
                relative_error = abs(p_new - p) / abs(p)
                if (relative_error < tolerance):
                    # relative error is less than tolerance
                    root = p_new
                    break

            if function_error < tolerance:
                # function bound less than tolerance
                root = p_new
                break

            # maximum number of iterations
            if iterations == 100:
                # exceeded maximum number of iterations i.e. 100
                root = p_new
                break

            # update root approximation
            p = p_new

    end = time.process_time()
    cpu_time = end - start

    # FILE WRITING:

    filepath = "FixedPointData.txt"

    with open(filepath, "w") as file:
        header = "FIXED POINT METHOD\n\niterations".center(10) + "root approximation".center(25) + "function value".center(25) + "|pn+1 - p|/|pn - p|".center(25) + "\n\n"
        file.write(header)
        for index, element in enumerate(root_approx):
            file.write(f"{index}".center(10) + f"{root_approx[index]}".center(25) + f"{f(root_approx[index])}".center(25))
            if index == 0:
                file.write("-".center(25) + "\n")
            else:
                expression = abs(root_approx[index] - root) / (abs(root_approx[index - 1] - root) ** 1)
                file.write(f"{expression}".center(25) + "\n")

    # PLOT GRAPH OF FUNCTION AND SUCCESSIVE ITERATIONS:

    # find the range for which the graph is to be plotted
    max = root
    min = root
    for element in root_approx:
        if element > max :
            max = element
        elif element < min :
            min = element

    # create an array of 100 x-values between the (min - 1) and (max + 1) value
    x2 = np.linspace(min - 1, max + 1, 100)
    y2 = []

    # creates an array for the function values at those x-values
    for element in x2:
        y2.append(f(element))

    # plots and labels the function, initial guess and the final root approximation
    plt.plot(x2, y2, label="Function")
    plt.plot([min - 1, max + 1], [min - 1, max + 1], label="f(x) = x")


    # plots all the root approximations and tangents
    for i in range(0, iterations):

        # plots a point at the nth root approximation
        plt.plot(root_approx[i], 0, 'o', label= f"p{i}")
        
        # plots a tangent to the function at the nth root approximation
        plt.plot([root_approx[i], root_approx[i], f(root_approx[i])], [0, f(root_approx[i]), f(root_approx[i])], '--', label="Iteration {}".format(i + 1))


    # plots the final root approx.
    plt.plot(root, f(root), 'x', label="Root")

    plt.legend()
    plt.title("Fixed Point method")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True)
    plt.show()

    return root, iterations, cpu_time


# Get input from user
# func1 = input("enter a function: ")
# p01 = float(input("Enter initial guess: "))
# tolerance1 = float(input("Enter tolerance: "))

func1 = "log(6 - (2**-x) - (2 * cos(x)))"
p01 = 1
tolerance1 = 10e-20


# Find root of function using Newton's method
root1, iterations1, cpu_time = FixedPointMethod(func1, p01, tolerance1)

# Print output
print(f"Root of function is: {root1}")
print(f"Number of iterations: {iterations1}")
print(f"The CPU time for convergence: {cpu_time} seconds")