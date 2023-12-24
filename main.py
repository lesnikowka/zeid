import numpy as np

nMax = 10000
s = 0
eps = 0.00000001
eps_max = 0.
eps_cur = 0.
n = 9
x = [0.]*n
x_old = 0.
x_new = 0.

a = [
    [-832., 16., 0., 400., 0., 0., 0., 0., 0.],
    [16., -832., 16., 0., 400., 0., 0., 0., 0.],
    [0., 16., -832., 0, 0., 400., 0., 0., 0.],
    [400., 0., 0., -832., 16., 0., 400., 0., 0.],
    [0., 400., 0., 16., -832., 16., 0., 400., 0.],
    [0., 0., 400., 0., 16., -832., 0., 0., 400.],
    [0., 0., 0., 400., 0., 0., -832., 16., 0.],
    [0., 0., 0., 0., 400., 0., 16., -832., 16.],
    [0., 0., 0., 0., 0., 400., 0., 16., -832.]
]

old_a = a.copy()

b = [
    -836.452,
    -846.7,
    -1011.952,
    -29.916,
    3.6,
    -42.916,
    -839.104,
    -849.3,
    -1014.604
]

gauss_sol = np.linalg.solve(a, b)

print("Решение методом Гаусса")
print(gauss_sol)
print()

def trueSol(x_, y_):
    return x_**3 + y_**3 + 2

#trueSolArr = [trueSol(x_i, y_i) for x_i in [0.25, 0.5, 0.75] for y_i in [0.05, 0.1, 0.15]]
trueSolArr2 = [trueSol(x_i, y_i)  for y_i in [0.05, 0.1, 0.15] for x_i in [0.25, 0.5, 0.75]]

print("Истинное решение: ")
#print(trueSolArr)
print(trueSolArr2)
print()

endCauseIsEps = False
endCauseIsS = False

end = False

while not end:
    eps_max = 0.
    for i in range(n):
        x_old = x[i]
        x_new = b[i]
        for j in range(n):
            if j != i:
                x_new -= a[i][j] * x[j]
        x_new /= a[i][i]
        eps_cur = abs(x_old - x_new)
        eps_max = max(eps_cur, eps_max)
        x[i] = x_new
    s += 1
    if eps_max <= eps:
        endCauseIsEps = True
        end = True
    if s >= nMax:
        endCauseIsS = True
        end = True


def fillInfo(text, data):
    for val in data:
        text = text.replace("%", str(val), 1)
    return text

ss = ""

if endCauseIsEps:
    ss += "Выход по точности\n"
if endCauseIsS:
    ss += "Выход по числу итераций\n"

residual = (np.array(old_a) @ np.array(x)) - np.array(b)

residual_abs = [abs(i) for i in residual]

inf_norm = max(residual_abs)

ss += "Норма невязки: " + str(inf_norm) + "\n"

ss += """При решении СЛАУ Ax=b методом Зейделя с критериями остановки Nmax=% и eps = %
 за S=% итераций достигнута точность eps_max=% и получено численное решение:\n\n"""

ss = fillInfo(ss, [nMax, eps, s, eps_max])

for i in range(n):
    ss += "x[" + str(i) + "] = " + str(x[i]) + "\n"







print(ss)