nMax = 10000
s = 0
eps = 0.0000001
eps_max = 0
eps_cur = 0
n = 4
x = [0]*n
x_old = 0
x_new = 0

a = [[1]*n]*n
b = [2]*n

while True:
    eps_max = 0
    for i in range(n):
        x_old = x[i]
        x_new = b[i]
        for j in range(n):
            if j != i:
                x_new -= a[i][j] * x[j]
            x_new /= a[i][i]
            eps_cur = abs(x_old - x_new)
            if eps_cur > eps_max:
                eps_max = eps_cur
            x[i] = x_new
    s += 1
    if eps_max <= eps or s >= nMax:
        break


def fillInfo(text, data):
    for val in data:
        text = text.replace("%", str(val), 1)
    return text

ss = """При решении СЛАУ Ax=b методом Зейделя с критериями остановки Nmax=% и eps = %
 за S=% итераций достигнута точность eps_max=% и получено численное решение:\n"""

ss = fillInfo(ss, [nMax, eps, s, eps_max])

for i in range(n):
    ss += "x[" + str(i) + "] = " + str(x[i]) + "\n"

print(ss)