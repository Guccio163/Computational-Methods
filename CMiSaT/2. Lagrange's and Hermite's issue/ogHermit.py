import numpy as np
import matplotlib.pyplot as plt
import math

wezly = 15

# "c" - podzial czebyszewa
# "r" - podzial na rownoodlegle punkty
podz = "c"


def f(x):
    return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)

def df(x):
    return x + 80 * math.sin(2 * x)

def factorial(x):
    res = 1
    for i in range(1, x + 1):
        res *= i
    return res


def a_s(x, y, dy):
    n = len(x)

    new_x = [x[i // 2] for i in range(2 * n)]

    matrix = [[None for _ in range(2 * n)] for __ in range(2 * n)]

    for i in range(n):
        matrix[2 * i][0] = y[i]
        matrix[2 * i + 1][0] = y[i]


    for i in range(1, 2 * n):
        for j in range(1, i + 1):
            if new_x[i] == new_x[i - j]:
                matrix[i][j] = dy[i // 2] / factorial(j)
            else:
                matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (new_x[i] - new_x[i - j])

    result = [matrix[i][i] for i in range(2 * n)]

    return result


def P(coef, x_values, x):
    product = 1
    result = coef[0]
    c_i = 1

    for i in range(0, len(x_values) - 1):
        product *= (x - x_values[i])
        result += product * coef[c_i]
        c_i += 1
        product *= (x - x_values[i])
        result += product * coef[c_i]
        c_i += 1

    product *= (x - x_values[len(x_values) - 1])
    result += product * coef[c_i]

    return result


def hermite(x, y, dy):
    coef = a_s(x, y, dy)

    x_result = np.linspace(min(x), max(x), 300)
    y_result = [P(coef, x, i) for i in x_result]

    error = 0
    max_her = 0
    for i in range(len(x_result)):
        error += (y_result[i] - f(x_result[i])) ** 2
        max_her = max(max_her, error)





    print(f" {wezly, max_her}")

    # print(f"błąd maksymalny: {max_her}")


    plt.plot(x, y, 'o', label = 'nodes')
    plt.plot(x_result, y_result, label = 'Hermit')
    plt.plot(x_result, [f(i) for i in x_result], label = 'Exact')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Interpolation with Hermit method for Chebyshev-spaced nodes')
    plt.show()


def equal_division(a, b, n):
    x_axis = np.linspace(a, b, n)
    return x_axis


def cheb_zeros(a, b, n):
    roots = np.zeros(n)
    for i in range(n):
        roots[i] = np.cos((2 * i + 1) / (2 * n) * np.pi)
    for i in range(n):
        roots[i] = ((b - a) / 2) * roots[i] + (a + b) / 2
    return roots


a, b = np.longdouble(-4*np.pi), np.longdouble(4*np.pi)

for wezly in range(10, 50):
    if podz == "c":
        x = cheb_zeros(a, b, wezly)
    else:
        x = equal_division(a, b, wezly)

    y = [f(i) for i in x]
    dy = [df(i) for i in x]

    hermite(x, y, dy)
