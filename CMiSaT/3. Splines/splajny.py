import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate


def f(x):
    return 10 * 4 + (x ** 2) / 2 - 10 * 4 * np.cos(2 * x)


def df(x):
    return x + 80 * np.sin(2 * x)


def f_star(xi, yi, xj, yj, x):
    psi1 = (xj - x) / (xj - xi)
    psi2 = (x - xi) / (xj - xi)
    return psi1 * yi + psi2 * yj


def h(i):
    return x_nodes[i + 1] - x_nodes[i]


def delta(i):
    return (y_nodes[i + 1] - y_nodes[i]) / h(i)


def sigmatable(end_var):
    b = [0 for _ in range(n)]
    for i in range(1, n - 1):
        b[i] = delta(i) - delta(i - 1)

    if (end_var == 2):
        b[0] = (df(x_nodes[1]) - df(x_nodes[0])) / 6 * (x_nodes[1] - x_nodes[0])
        b[n - 1] = (df(x_nodes[n - 1]) - df(x_nodes[n - 2])) / 6 * (x_nodes[n - 1] - x_nodes[n - 2])

    A = [[0 for _ in range(n)] for __ in range(n)]
    A[0][0] = 1
    for i in range(1, n - 1):
        A[i][i - 1] = h(i - 1)
        A[i][i] = 2 * (h(i - 1) + h(i))
        A[i][i + 1] = h(i)
    A[n - 1][n - 1] = 1
    return np.linalg.solve(A, b)


# def sigma(table):
#     for i in range(1, len(table) - 1):
#         table[i] = (delta(i) - delta(i - 1) - h(i - 1) * table[i - 1] - h(i) * table[i + 1]) / (2 * (h(i - 1) + h(i)))


def s_c(x, i, varr):
    if varr == 1:
        return y_nodes[i] + b_c_1[i] * (x - x_nodes[i]) + c_c_1[i] * (x - x_nodes[i]) ** 2 + d_c_1[i] * (
                    x - x_nodes[i]) ** 3
    else:
        return y_nodes[i] + b_c_2[i] * (x - x_nodes[i]) + c_c_2[i] * (x - x_nodes[i]) ** 2 + d_c_2[i] * (
                    x - x_nodes[i]) ** 3


def s_q(x, i, varr):
    if varr == 1:
        return y_nodes[i] + b_q_1[i] * (x - x_nodes[i]) + c_q_1[i] * (x - x_nodes[i]) ** 2
    else:
        return y_nodes[i] + b_q_2[i] * (x - x_nodes[i]) + c_q_2[i] * (x - x_nodes[i]) ** 2


a = -4 * np.pi
b = 4 * np.pi
nodecount = 20
# tabele na max errory dla każdej ilości węzłów od 5 do zadanej
errors_c_nat = [0 for _ in range(4, nodecount)]
errors_c_cla = [0 for _ in range(4, nodecount)]
errors_q_nat = [0 for _ in range(4, nodecount)]
errors_q_cla = [0 for _ in range(4, nodecount)]

opt = int(input("1 for cubic, 2 for quadratic, 3 for both"))
end_var = int(input("1 for natural 2 for clamped, 3 for both"))

for n in range(5, nodecount+1):
    x_nodes = np.linspace(a, b, n)
    y_nodes = f(x_nodes)
    xs = np.linspace(a, b, 300)
    ys_c_1 = [0 for _ in range(len(xs))]
    ys_c_2 = [0 for _ in range(len(xs))]
    ys_q_1 = [0 for _ in range(len(xs))]
    ys_q_2 = [0 for _ in range(len(xs))]
    y_true = f(xs)


    ######################## CUBIC


    if end_var == 1 or end_var == 3:

        sigma_table_n = sigmatable(1)

        # b, c i d są określone dla przedziałów więc te tablice mają długości n-1
        b_c_1 = [0 for _ in range(n - 1)]
        c_c_1 = [0 for _ in range(n - 1)]
        d_c_1 = [0 for _ in range(n - 1)]

        for i in range(n - 1):
            b1 = (y_nodes[i + 1] - y_nodes[i]) / h(i)
            b2 = h(i) * (sigma_table_n[i + 1] + 2 * sigma_table_n[i])
            b_c_1[i] = b1 - b2
            c_c_1[i] = 3 * sigma_table_n[i]
            d_c_1[i] = (sigma_table_n[i + 1] - sigma_table_n[i]) / h(i)


    if end_var == 2 or end_var == 3:

        sigma_table_c = sigmatable(2)

        b_c_2 = [0 for _ in range(n - 1)]
        c_c_2 = [0 for _ in range(n - 1)]
        d_c_2 = [0 for _ in range(n - 1)]

        for i in range(n - 1):
            b1 = (y_nodes[i + 1] - y_nodes[i]) / h(i)
            b2 = h(i) * (sigma_table_c[i + 1] + 2 * sigma_table_c[i])
            b_c_2[i] = b1 - b2
            c_c_2[i] = 3 * sigma_table_c[i]
            d_c_2[i] = (sigma_table_c[i + 1] - sigma_table_c[i]) / h(i)

    ######################################## CUBIC END
    ######################################## QUADRATIC


    b_q_1 = [0 for _ in range(n)]
    b_q_2 = [0 for _ in range(n)]
    c_q_1 = [0 for _ in range(n - 1)]
    c_q_2 = [0 for _ in range(n - 1)]

    b_q_1[n - 1] = 0  # natural spline (left) and clamped spline (down)
    b_q_2[n - 1] = (y_nodes[n - 1] - y_nodes[n - 2]) / 6 * (x_nodes[n - 1] - x_nodes[n - 2])

    for p in range(n - 2, -1, -1):
        b_q_1[p] = 2 * (y_nodes[p + 1] - y_nodes[p]) / (x_nodes[p + 1] - x_nodes[p]) - b_q_1[p + 1]
        b_q_2[p] = 2 * (y_nodes[p + 1] - y_nodes[p]) / (x_nodes[p + 1] - x_nodes[p]) - b_q_2[p + 1]

    for p in range(n - 1):
        c_q_1[p] = (b_q_1[p + 1] - b_q_1[p]) / (2 * (x_nodes[p + 1] - x_nodes[p]))
        c_q_2[p] = (b_q_2[p + 1] - b_q_2[p]) / (2 * (x_nodes[p + 1] - x_nodes[p]))

    ######################################## QUADRATIC END

    j = 1


    max_c_nat = 0
    max_c_cla = 0
    max_q_nat = 0
    max_q_cla = 0




    if opt == 1:
        if end_var == 1:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_1[k] = s_c(xs[k], j - 1, 1)
            plt.plot(xs, ys_c_1, label="3-go stopnia natural")
        elif end_var == 2:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_2[k] = s_c(xs[k], j - 1, 2)
            plt.plot(xs, ys_c_2, label="3-go stopnia clamped")
        else:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_1[k] = s_c(xs[k], j - 1, 1)
                ys_c_2[k] = s_c(xs[k], j - 1, 2)
            plt.plot(xs, ys_c_1, label="3-go stopnia natural")
            plt.plot(xs, ys_c_2, label="3-go stopnia clamped")

    elif opt == 2:
        if end_var == 1:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_q_1[k] = s_q(xs[k], j - 1, 1)
            plt.plot(xs, ys_q_1, label="2-go stopnia natural")
        elif end_var == 2:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_q_2[k] = s_q(xs[k], j - 1, 2)
            plt.plot(xs, ys_q_2, label="2-go stopnia clamped")
        else:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_q_1[k] = s_q(xs[k], j - 1, 1)
                ys_q_2[k] = s_q(xs[k], j - 1, 2)
            plt.plot(xs, ys_q_1, label="2-go stopnia natural")
            plt.plot(xs, ys_q_2, label="2-go stopnia clamped")



    else:
        if end_var == 1:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_1[k] = s_c(xs[k], j - 1, 1)
                ys_q_1[k] = s_q(xs[k], j - 1, 1)
            plt.plot(xs, ys_c_1, label="3-go stopnia natural")
            plt.plot(xs, ys_q_1, label="2-go stopnia natural")
        elif end_var == 2:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_2[k] = s_c(xs[k], j - 1, 2)
                ys_q_2[k] = s_q(xs[k], j - 1, 2)
            plt.plot(xs, ys_c_2, label="3-go stopnia clamped")
            plt.plot(xs, ys_q_2, label="2-go stopnia clamped")
        else:
            for k in range(len(ys_c_1)):
                if xs[k] > x_nodes[j]:
                    j += 1
                ys_c_1[k] = s_c(xs[k], j - 1, 1)
                ys_c_2[k] = s_c(xs[k], j - 1, 2)
                ys_q_1[k] = s_q(xs[k], j - 1, 1)
                ys_q_2[k] = s_q(xs[k], j - 1, 2)
            plt.plot(xs, ys_c_1, label="3-go stopnia natural")
            plt.plot(xs, ys_c_2, label="3-go stopnia clamped")
            plt.plot(xs, ys_q_1, label="2-go stopnia natural")
            plt.plot(xs, ys_q_2, label="2-go stopnia clamped")

    plt.plot(xs, y_true, label="funkcja interpolowana")
    plt.plot(x_nodes, y_nodes, ".", label="węzły")
    plt.legend()
    plt.title("Interpolacja funkcją sklejaną dla "+str(n)+" węzłów")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()


    for i in range(len(xs)):
        if opt==1:
            if end_var==1 or end_var == 3:
                max_c_nat = max(abs(ys_c_1[i] - y_true[i]), max_c_nat)
            if end_var==2 or end_var == 3:
                max_c_cla = max(abs(ys_c_2[i] - y_true[i]), max_c_cla)

        elif opt==2:
            if end_var==1 or end_var==3:
                max_q_nat = max(abs(ys_q_1[i] - y_true[i]), max_q_nat)
            if end_var==2 or end_var==3:
                max_q_cla = max(abs(ys_q_2[i] - y_true[i]), max_q_cla)

        else:
            if end_var==1 or end_var == 3:
                max_c_nat = max(abs(ys_c_1[i] - y_true[i]), max_c_nat)
                max_q_nat = max(abs(ys_q_1[i] - y_true[i]), max_q_nat)
            if end_var==2 or end_var == 3:
                max_c_cla = max(abs(ys_c_2[i] - y_true[i]), max_c_cla)
                max_q_cla = max(abs(ys_q_2[i] - y_true[i]), max_q_cla)

    if opt == 1:
        if end_var == 1 or end_var == 3:
            errors_c_nat[n-5] = max_c_nat
        elif end_var == 2 or end_var == 3:
            errors_c_cla[n-5] = max_c_cla
    elif opt == 2:
        if end_var == 1 or end_var == 3:
            errors_q_nat[n-5] = max_q_nat
        if end_var == 2 or end_var == 3:
            errors_q_cla[n-5] = max_q_cla
    else:
        if end_var == 1 or end_var == 3:
            errors_c_nat[n-5] = max_c_nat
            errors_q_nat[n-5] = max_q_nat
        if end_var == 2 or end_var == 3:
            errors_c_cla[n-5] = max_c_cla
            errors_q_cla[n-5] = max_q_cla





intigers = [i for i in range(5, nodecount+1)]
if opt == 1:
    if end_var == 1:
        table_combined = [list(x) for x in zip(intigers, errors_c_nat)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural c"]))
    elif end_var == 2:
        table_combined = [list(x) for x in zip(intigers, errors_c_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max clamped c"]))
    else:
        table_combined = [list(x) for x in zip(intigers, errors_c_nat, errors_c_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural c", "Błąd max clamped c"]))
elif opt == 2:
    if end_var == 1:
        table_combined = [list(x) for x in zip(intigers, errors_q_nat)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural q"]))
    elif end_var == 2:
        table_combined = [list(x) for x in zip(intigers, errors_q_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max clamped q"]))
    else:
        table_combined = [list(x) for x in zip(intigers, errors_q_nat, errors_q_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural q", "Błąd max clamped q"]))
else:
    if end_var == 1:
        table_combined = [list(x) for x in zip(intigers, errors_c_nat, errors_q_nat)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural c", "Bład max natural q"]))
    elif end_var == 2:
        table_combined = [list(x) for x in zip(intigers, errors_c_cla, errors_q_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max clamped c", "Błąd max clamped q"]))
    else:
        table_combined = [list(x) for x in zip(intigers, errors_c_nat, errors_c_cla, errors_q_nat, errors_q_cla)]
        print(tabulate(table_combined, headers=["Liczba węzłów", "Błąd max natural c", "Błąd max clamped c", "Błąd max natural q", "Błąd max clamped q"]))

