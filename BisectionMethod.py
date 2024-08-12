import math


def bisection_method(func, a, b, tolerance):

    # returns the function value at x
    def f(x):
        return eval(func, {"__builtins__": None}, {"x": x, "sin": math.sin, "cos": math.cos, "tan": math.cos, "pi": math.pi, "e": math.e})

    if f(a) == 0:
        print(f"root already exists at {a}")
        return

    if f(b) == 0:
        print(f"root already exists at {b}")
        return

    while f(a) * f(b) > 0:
        print("none or multiple roots!\nplease try again:")
        a = float(input("lower boundary: "))
        b = float(input("upper boundary: "))

    filepath = "DATA.txt"

    with open(filepath, "w") as file:
        header = "iterations".center(10) + "interval".center(50) + "root approximation".center(25) + "function value".center(25) + "relative error".center(25) + "\n\n"
        file.write(header)

    i = 0   # no. of iterations
    pn = []
    init_a = a
    init_b = b

    # do-while
    while True:

        i = i + 1   # no. of iterations incremented

        # writing iteration number and interval in the file
        with open(filepath, "a") as file:
            file.write(f"{i}".center(10) + f"{a} to {b}".center(50))

        # previous term exists only after the first iteration
        if i > 1:
            prev_term = c

        c = a + (b - a) / 2   # mid-point of a and b

        # nth approximation
        pn.append(c)

        if f(c) * f(a) < 0:
            # root exists on the left-hand-side interval
            b = c

        elif f(c) * f(b) < 0:
            # root exists on the right-hand-side interval
            a = c

        # stopping criteria:

        if i == 1:
            error = abs(f(c))  # when there is no previous term to compute the relative error (first iteration)

        else:
            if c == 0:  # in order to prevent the divide-by-zero error
                error = abs(c - prev_term)
            else:
                error = abs((c - prev_term) / c)    # relative error

        # writing root approx. and its function value (+ the relative error) in the file
        with open(filepath, "a") as file:
            file.write(f"{c}".center(25) + f"{f(c)}".center(25) + f"{error}".center(25) + "\n")

        if abs(f(c)) < tolerance:
            # bound on the function value
            break

        if error < tolerance:
            # relative error
            break

        # maximum number of iterations set as 1000
        if i == 1000:
            print("maximum number of iterations reached!")
            break

    print(f"root: {c}")
    print(f"no. of iterations: {i}")

    with open(filepath, "a") as file:
        header = f"\nroot: {c} ; no. of iterations: {i}\n\nFor order of convergence and asymptotic convergence constant: \n\niterations".center(10) + "|pn - p|".center(25) + "((b - a) / 2^n)".center(25) + "\n\n"
        file.write(header)

    j = 1
    while j <= i:
        error_term = abs(pn.__getitem__(j - 1) - c)
        upper_bound = (init_b - init_a) / 2 ** j
        with open(filepath, "a") as file:
            file.write(f"{j}".center(10) + f"{error_term}".center(25) + f"{upper_bound}".center(25) + "\n")

        j = j + 1

    with open(filepath, "a") as file:
        file.write("\n(b - a)/2^n is clearly the upper bound for this error as illustrated in the table above\ntherefore the rate of convergence is O(1/(2^n))\nsince the interval is halved at each iteration, it linearly converges to the root with the asymptotic convergence constant 1/2")

    # (b - a)/2^n is clearly the upper bound for this error as illustrated in the table above
    # therefore the rate of convergence is O(1/(2^n))
    # since the interval is halved at each iteration it linearly converges to the root
    # and the asymptotic convergence constant is 1/2

    return c


bisection_method("x - cos(x)", 0, math.pi, 0.0001)

# assignment 1, question 6:
# bisection_method("(e ** x) - (x ** 2) + (3*x) - 2", 0, 1, 0.00001)

# in order to input a function, interval and tolerance from the user and apply the bisection method:

# func2 = input("enter a function: ")
# a2 = float(input("lower boundary: "))
# b2 = float(input("upper boundary: "))
# tolerance2 = float(input("tolerance: "))
#
# bisection_method(func2, a2, b2, tolerance2)
