import math
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate



####################################################################################################
################################## FUNKCJE INTERPOLUJĄCE ###########################################
####################################################################################################


def a_s(xs, ys):
    n = len(ys)
    As = ys.copy()
    for j in range(1, n):
        for k in range(n - 1, j - 1, -1):  # ***
            As[k] = (As[k] - As[k - 1]) / (xs[k] - xs[k - j])
    return As


def my_newton(xs, ys, x):
    n = len(ys)
    a = a_s(xs, ys)
    iloczyns = [1 for _ in range(n)]
    for i in range(1, n):
        iloczyns[i] = (x - xs[i - 1]) * iloczyns[i - 1]

    results = [a[i] * iloczyns[i] for i in range(n)]
    res = 0
    for i in range(n):
        res += results[i]
    return res


# (***) robimy od końca do początku bo obliczając wartość na danym polu chcemy korzystać
# z tego co było na nim i na niższym w poprzednim wykonaniu pętli; jeśli będziemy
# iść od dołu to każde kolejne pole będzie czerpało z już zmienionego niższego,
# kiedy potrzebujemy wartości z poprzedniego przejścia


def L(k, x, xs):
    a, b = 1, 1
    for i in range(len(xs)):
        if i == k:
            continue
        a *= x - xs[i]
        b *= xs[k] - xs[i]
    return a / b


def my_lagrange(xs, ys, x):
    results = [ys[k] * L(k, x, xs) for k in range(n)]
    res = 1
    for i in range(n):
        res += results[i]

    return res


def L_polynomial(k, xs):
    a = np.polynomial.Polynomial(1)
    b = 1
    for i in range(len(xs)):
        if i==k:
            continue
        a *= np.polynomial.Polynomial([-xs[i], 1])
        b *= xs[k]-xs[i]
    return a, b


def lag_polynomial(xs, ys):
    n = len(xs)
    p = np.polynomial.Polynomial(0)
    for i in range(n):
        a, b = L_polynomial(i, xs)
        p += ys[i]/b*a
    return p


def chebyshev_nodes(a, b, n):
    k = np.arange(1, n + 1)
    x = np.cos((2 * k - 1) * np.pi / (2 * n))  # węzły Czebyszewa na przedziale [-1,1]
    x = (x + 1) * (b - a) / 2 + a  # przekształcenie liniowe na [a, b]
    return x



def my_function(x):
    return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)


####################################################################################################
############################### FUNKCJE RYSUJĄCE WYKRESY ###########################################
####################################################################################################

def plot_points(xx):
    yy = [1 for _ in range(len(xx))]
    for i in range(len(xx)):
        yy[i] = my_function(xx[i])
    plt.plot(xx, yy, label='function')


def plot_nodes(nodes, nodes_y):
    plt.plot(nodes, nodes_y, '.', label='function nodes')


def plot_mynewton(x, y, xx):
    y_new = [0 for _ in range(len(xx))]
    for i in range(len(xx)):
        y_new[i] = my_newton(x, y, xx[i])

    plt.plot(xx, y_new, label='newton')


def plot_lagrangePOL(x, y, xx):
    y_lagpol = [0 for _ in range(len(xx))]
    for i in range(len(xx)):
        y_lagpol[i] = lag_polynomial(x, y)(xx[i])

    plt.plot(xx, y_lagpol, label='lagrange')


def print_evr(x, opt):


    max_node = 20
    errors = [[0,0,0] for _ in range(max_node - 4)]
    for j in range(5, max_node):
        if opt == 1:
            nodes = np.linspace(a, b, j)
        else:
            nodes = chebyshev_nodes(a, b, j)
        nodes_y = [0 for _ in range(j)]
        for i in range(j):
            nodes_y[i] = my_function(nodes[i])
        max_lag = 0
        max_new = 0
        lag_pol = lag_polynomial(nodes, nodes_y)
        for i in range(200):
            lagrange_value = lag_pol(x[i])
            newton_value = my_newton(nodes, nodes_y, x[i])
            exact_value = my_function(x[i])
            lag_err = abs(lagrange_value-exact_value)
            new_err = abs(newton_value-exact_value)
            max_lag = max(max_lag, lag_err)
            max_new = max(max_new, new_err)
        errors[j-5] = [j, max_lag, max_new]
        print(errors[j-5])
        if(j == max_node-1):
            polyn = [[0,0] for _ in range(max_node)]
            for i in range(max_node-1):
                polyn[i] = [i, lag_pol.coef[i]]
            head = ['x\'s power', 'coefficient' ]
            print(tabulate(polyn, headers = head))

    # headers = ['nodecount', 'Lag', 'New']
    #
    # with open('obliczenia.txt', 'w') as f:
    #     print(tabulate(errors, headers=headers), file = f)

        # print('max error for lagrange: ', max_lag, " max error for newton: ", max_new, "\n", file=f)

    plot_points(x)
    plot_nodes(nodes, nodes_y)
    #plot_lagrangePOL(nodes, nodes_y, x)
    plot_mynewton(nodes, nodes_y, x)
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Interpolation with Newton method for 20 Chebyshev-spaced nodes')
    plt.show()

def print_pol():
    pol_arr = []
    sides = []
    for i in range(5,40):
        x = chebyshev_nodes(-4 * np.pi, 4 * np.pi, i)
        y = [1 for _ in range(i)]
        for k in range(i):
            y[k] = my_function(x[k])
        pol_i = lag_polynomial(x, y)
        pol_arr.append(pol_i)
        side = [0, i]
        for j in range(len(x)):
            side[0] += abs(pol_i(x[j])-my_function(x[j]))
            sides.append(side)
    print(min(sides)[0], "\n", pol_arr[min(sides)[1]])


####################################################################################################
###################################### WYKONANIE FUNKCJI ###########################################
####################################################################################################

a = -4 * np.pi
b = 4 * np.pi
n = 20
# nodes_lin = np.linspace(a, b, n)
# nodes_ch = chebyshev_nodes(a, b, n)
x = np.linspace(a, b, 200)
nodes_y = [0 for i in range(0, n)]      # placeholder na wyniki zależnie od rozłożenia nodów

opt = int(input("opcja 1: linspace, opcja 2: chebyshev"))

print_evr(x, opt)
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #
# #

import numpy as np
import matplotlib.pyplot as plt
import math

#
# # "c" - podzial czebyszewa
# # "r" - podzial na rownoodlegle punkty
# podz = "c"
#
#
# def f(x):
#     return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)
#
#
# def df(x):
#     return x + math.sin(2 * x)
#
#
# def factorial(x):
#     res = 1
#     for i in range(1, x + 1):
#         res *= i
#     return res
#
#
# ########################################################################################################
#
#
# def a_s2(xs, ys):
#     n = len(ys)
#     As = ys.copy()
#
#     for i in range(1, n):
#         xd = xs[i]
#         As[i] = dy(xd)
#
#     for j in range(2, n):
#         for k in range(n - 1, j - 1, -1):  # ***
#             As[k] = (As[k] - As[k - 1]) / (xs[k] - xs[k - j])
#     return As
#
#
# def my_newton2(xs, ys, x):
#     n = len(ys)
#     a = a_s2(xs, ys)
#     iloczyns = [1 for _ in range(n)]
#     for i in range(1, n):
#         iloczyns[i] = (x - xs[i - 1]) * iloczyns[i - 1]
#
#     results = [a[i] * iloczyns[i] for i in range(n)]
#     res = 0
#     for i in range(n):
#         res += results[i]
#     return res
#
#
#
# ###################################################################################################
#
#
# def get_matrix(x, y, dy):
#     n = len(x)
#
#     new_x = [x[i // 2] for i in range(2 * n)]       # new_x = [0,0,1,1,2,2,3,3,4,4,5,5 ... n-2,n-2,n-1,n-1]
#     matrix = [[None for _ in range(2 * n)] for __ in range(2 * n)]
#
#     for i in range(2*n):
#         matrix[i][0] = y[i//2]        # pierwsza kolumna macierzy to [0,0,1,1,2,2,3,3,4,4 ... n-2,n-2,n-1,n-1]
#
#     for i in range(1, 2 * n):           # od drugiej kolumny bo pierwszą już mamy
#         for j in range(1, i + 1):               # robimy tylko dolny trójkąt do przekątnej, tyle potrzebujemy
#             if new_x[i] == new_x[i - j]:
#                 matrix[i][j] = dy[i // 2] / factorial(j)    # pochodna i//2  /  potęga i
#             else:
#                 matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (new_x[i] - new_x[i - j])
#                                                             # różnica nad nami i po skosie od nas/ różnice x[i] i x[i-j]
#
#     result = [matrix[i][i] for i in range(2 * n)]       # bierzemy przekątną macierzy
#
#     return result
#
#
#
# def P(coef, x_values, x):
#     product = 1
#     result = coef[0]
#     c_i = 1
#
#     for i in range(0, len(x_values) - 1):
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#
#     product *= (x - x_values[len(x_values) - 1])
#     result += product * coef[c_i]
#
#     return result
#
#
# def hermite(x, y, dy):
#     coef = get_matrix(x, y, dy)         # współczynniki
#
#     x_result = np.linspace(min(x), max(x), 300)
#     y_result = [P(coef, x, i) for i in x_result]
#
#
#     error = 0
#     for i in range(len(x_result)):
#         error += (y_result[i] - f(x_result[i])) ** 2
#
#     max_dev = -99999999
#
#     for i in range(len(x_result)):
#         max_dev = max(max_dev, abs(y_result[i] - f(x_result[i])))
#
#     error /= 300
#     error = error ** 0.5
#
#     # print(f"wezly: {wezly}")
#     # print(f"blad: {error}")
#     # print(f"odch: {max_dev}")
#     # print("=============================")
#
#     plt.plot(x_result, y_result, label = 'inter')
#     plt.scatter(x, y, label='nodes')
#     plt.plot(x_result, [f(i) for i in x_result], label='kox')
#     plt.legend()
#     plt.show()
#
#
#
# def chebyshev_nodes(a, b, n):
#     k = np.arange(1, n + 1)
#     x = np.cos((2 * k - 1) * np.pi / (2 * n))  # węzły Czebyszewa na przedziale [-1,1]
#     x = (x + 1) * (b - a) / 2 + a  # przekształcenie liniowe na [a, b]
#     return x
#
#
#
# a = -4 * np.pi
# b = 4 * np.pi
# n = 17
#
# nodes = chebyshev_nodes(a, b, n)
# # nodes = np.linspace(a, b, n)
# y = np.array([f(nodes[i]) for i in range(len(nodes))])
# dy = np.array([df(nodes[i]) for i in range(len(nodes))])
#
#
# hermite(nodes, y, dy)
#
# # plt.scatter(nodes, y, color='red', label='Węzły')
# # plt.scatter(nodes, dy, color='green', label='pochodna')
# # plt.plot(x, y_exact, color = 'lightgreen')
# # plt.legend()
# # plt.show()







#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################


import numpy as np
import matplotlib.pyplot as plt
import math
from tabulate import tabulate



# def factorial(x):
#     res = 1
#     for i in range(1, x + 1):
#         res *= i
#     return res
#
# def my_function(x):
#     return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)
#
# def d_my_function(x):
#     return x + 80* math.sin(2 * x)
#
# def a_s(x, y, dy):
#     n = len(x)
#
#     new_x = [x[i // 2] for i in range(2 * n)]
#
#     matrix = [[0 for _ in range(2 * n)] for __ in range(2 * n)]
#
#     for i in range(2*n):
#         matrix[i][0] = y[i//2]
#
#     for i in range(1, 2 * n):
#         for j in range(1, i + 1):
#             if new_x[i] == new_x[i - j]:
#                 matrix[i][j] = dy[i // 2] / factorial(j)
#             else:
#                 matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (new_x[i] - new_x[i - j])
#
#     a = [matrix[i][i] for i in range(2 * n)]
#
#     return a
#
#
# def Hermit(a, x_values, x):
#     product = 1
#     result = a[0]
#     c_i = 1
#
#     for i in range(0, len(x_values) - 1):
#         product *= (x - x_values[i])
#         result += product * a[c_i]
#         c_i += 1
#         product *= (x - x_values[i])
#         result += product * a[c_i]
#         c_i += 1
#
#     product *= (x - x_values[len(x_values) - 1])
#     result += product * a[c_i]
#
#     return result
#
#
#
# def chebyshev_nodes(a, b, n):
#     k = np.arange(1, n + 1)
#     x = np.cos((2 * k - 1) * np.pi / (2 * n))  # węzły Czebyszewa na przedziale [-1,1]
#     x = (x + 1) * (b - a) / 2 + a  # przekształcenie liniowe na [a, b]
#     return x
#
#
# #  0 - podzial na rownoodlegle punkty
# #  1 - podzial czebyszewa
#
# podz = 0
#
#
# a = -4 * np.pi
# b = 4 * np.pi
# max_nodes = 17
#
# # errors = [[0,0] for _ in range(max_nodes-4)]
# # max_her = 0
#
# for i in range(5, max_nodes):
#     if podz == 1:
#         xs = chebyshev_nodes(a, b, i)
#     else:
#         xs = np.linspace(a, b, i)
#
#     ys = [0 for x in xs]
#     for k in range(len(xs)):
#         ys[k] = my_function(xs[k])
#     dys = [0 for x in xs]
#     for k in range(len(xs)):
#         dys[k] = d_my_function(xs[k])
#     x_result = np.linspace(min(xs), max(xs), 300)
#
#     a = a_s(xs, ys, dys)
#     y_result = [Hermit(a, xs, x) for x in x_result]
#
#     # error = 0
#     # for j in range(len(x_result)):
#     #     error = max(abs(y_result[j] - my_function(x_result[j])), error)
#     # errors[i-5] = [i, error]
#     # max_her = max(max_her, error)
#
#     # print(f"wezly: {wezly}")
#     # print(f"blad: {error}")
#     # print(f"odch: {max_dev}")
#     # print("=============================")
#
#
#
#     plt.plot(xs, ys, 'o', label='nodes')
#     plt.plot(x_result, y_result, label='Hermit')
#     plt.plot(x_result, [my_function(i) for i in x_result], label='exact')
#     plt.legend()
#     plt.show()
#
# # headers = ['nodecount', 'Lag', 'New']
# #
# # with open('obliczenia2.txt', 'w') as f:
# #     print(tabulate(errors, headers=headers), file = f)
# #
# #     print('max error for Hermit: ', max_her, "\n", file=f)


# import numpy as np
# import matplotlib.pyplot as plt
# import math
#
#
#
#
# def my_function(x):
#     return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)
#
# def d_my_function(x):
#     return x + 80* math.sin(2 * x)
#
#
# def factorial(x):
#     res = 1
#     for i in range(1, x + 1):
#         res *= i
#     return res
#
#
# def a_s(x, y, dy):
#     n = len(x)
#
#     new_x = [x[i // 2] for i in range(2 * n)]
#
#     matrix = [[None for _ in range(2 * n)] for __ in range(2 * n)]
#
#     for i in range(n):
#         matrix[2 * i][0] = y[i]
#         matrix[2 * i + 1][0] = y[i]
#
#     for i in range(1, 2 * n):
#         for j in range(1, i + 1):
#             if new_x[i] == new_x[i - j]:
#                 matrix[i][j] = dy[i // 2] / factorial(j)
#             else:
#                 matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (new_x[i] - new_x[i - j])
#
#     result = [matrix[i][i] for i in range(2 * n)]
#
#     return result
#
#
# def Hermit(coef, x_values, x):
#     product = 1
#     result = coef[0]
#     c_i = 1
#
#     for i in range(0, len(x_values) - 1):
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#
#     product *= (x - x_values[len(x_values) - 1])
#     result += product * coef[c_i]
#
#     return result
#
#
#
#
# def cheb_zeros(a, b, n):
#     roots = np.zeros(n)
#     for i in range(n):
#         roots[i] = np.cos((2 * i + 1) / (2 * n) * np.pi)
#     for i in range(n):
#         roots[i] = ((b - a) / 2) * roots[i] + (a + b) / 2
#     return roots
#
#
# #  0 - podzial na rownoodlegle punkty
# #  1 - podzial czebyszewa
#
# podz = 1
#
# a = -4 * np.pi
# b = 4 * np.pi
# max_nodes = 17
#
#
# for wezly in range(5, max_nodes):
#     if podz == 1:
#         x = cheb_zeros(a, b, wezly)
#     else:
#         x = np.linspace(a, b, wezly)
#
#     y = [my_function(i) for i in x]
#     dy = [d_my_function(i) for i in x]
#
#     a = a_s(x, y, dy)
#
#     x_result = np.linspace(min(x), max(x), 300)
#     y_result = [Hermit(a, x, i) for i in x_result]
#
#     error = 0
#     for i in range(len(x_result)):
#         error += (y_result[i] - my_function(x_result[i])) ** 2
#
#     max_dev = -99999999
#
#     for i in range(len(x_result)):
#         max_dev = max(max_dev, abs(y_result[i] - my_function(x_result[i])))
#
#     error /= 300
#     error = error ** 0.5
#
#     print(f"wezly: {wezly}")
#     print(f"blad: {error}")
#     print(f"odch: {max_dev}")
#     print("=============================")
#
#     plt.plot(x, y, 'ro', x_result, y_result)
#     plt.plot(x_result, [my_function(i) for i in x_result])
#     plt.show()


# import numpy as np
# import matplotlib.pyplot as plt
# import math
#
# wezly = 15
#
# # "c" - podzial czebyszewa
# # "r" - podzial na rownoodlegle punkty
# podz = "r"
#
#
# def my_function(x):
#     return 10 * 4 + (x ** 2) / 2 - 10 * 4 * math.cos(2 * x)
#
# def d_my_function(x):
#     return x + 80* math.sin(2 * x)
#
#
# def factorial(x):
#     res = 1
#     for i in range(1, x + 1):
#         res *= i
#     return i
#
#
# def get_matrix(x, y, dy):
#     n = len(x)
#
#     new_x = [x[i // 2] for i in range(2 * n)]
#
#     matrix = [[None for _ in range(2 * n)] for __ in range(2 * n)]
#
#     for i in range(n):
#         matrix[2 * i][0] = y[i]
#         matrix[2 * i + 1][0] = y[i]
#
#     for i in range(1, 2 * n):
#         for j in range(1, i + 1):
#             if new_x[i] == new_x[i - j]:
#                 matrix[i][j] = dy[i // 2] / factorial(j)
#             else:
#                 matrix[i][j] = (matrix[i][j - 1] - matrix[i - 1][j - 1]) / (new_x[i] - new_x[i - j])
#
#     result = [matrix[i][i] for i in range(2 * n)]
#
#     return result
#
#
# def P(coef, x_values, x):
#     product = 1
#     result = coef[0]
#     c_i = 1
#
#     for i in range(0, len(x_values) - 1):
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#         product *= (x - x_values[i])
#         result += product * coef[c_i]
#         c_i += 1
#
#     product *= (x - x_values[len(x_values) - 1])
#     result += product * coef[c_i]
#
#     return result
#
#
# def hermite(x, y, dy):
#     coef = get_matrix(x, y, dy)
#
#     x_result = np.linspace(min(x), max(x), 300)
#     y_result = [P(coef, x, i) for i in x_result]
#
#     error = 0
#     for i in range(len(x_result)):
#         error += (y_result[i] - f(x_result[i])) ** 2
#
#     max_dev = -99999999
#
#     for i in range(len(x_result)):
#         max_dev = max(max_dev, abs(y_result[i] - f(x_result[i])))
#
#     error /= 300
#     error = error ** 0.5
#
#     print(f"wezly: {wezly}")
#     print(f"blad: {error}")
#     print(f"odch: {max_dev}")
#     print("=============================")
#
#     plt.plot(x, y, 'ro', x_result, y_result)
#     plt.plot(x_result, [f(i) for i in x_result])
#     plt.show()
#
#
# def equal_division(a, b, n):
#     x_axis = np.linspace(a, b, n)
#     return x_axis
#
#
# def cheb_zeros(a, b, n):
#     roots = np.zeros(n)
#     for i in range(n):
#         roots[i] = np.cos((2 * i + 1) / (2 * n) * np.pi)
#     for i in range(n):
#         roots[i] = ((b - a) / 2) * roots[i] + (a + b) / 2
#     return roots
#
#
# a, b = np.longdouble(-7), np.longdouble(7)
#
# for wezly in range(30, 50):
#     if podz == "c":
#         x = cheb_zeros(a, b, wezly)
#     else:
#         x = equal_division(a, b, wezly)
#
#     y = [f(i) for i in x]
#     dy = [df(i) for i in x]
#
#     hermite(x, y, dy)
