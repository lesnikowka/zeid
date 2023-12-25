import numpy as np

nMax = 100
s = 0
eps = 1e-8
endCauseIsEps = False
endCauseIsS = False
end = False

A = 0
B = 1
C = 0
D = 0.2
M = 4
N = 4
h = (B-A)/N
k = (D-C)/M

def trueSol(x_, y_):
    return x_**3 + y_**3 + 2
def u(x_, y_):
    return trueSol(x_, y_)
def u0y(y_):
    return y_**3 + 2
def u1y(y_):
    return y_**3 + 3
def ux0(x_):
    return x_**3 + 2
def ux0_2(x_):
    return x_**3 + 251./125.
def f(x_, y_):
    return -6 * (x_ + y_)

#Создаем нулевую матрицу V
v = [ ([0.]*(M+1)).copy() for i in range(N + 1)]

#Заполняем граничные значения
for i in range(M + 1):
    v[0][i] = u0y(i * k)
    v[len(v)-1][i] = u1y(i * k)

for i in range(N + 1):
    v[i][0] = ux0(i * h)
    v[i][len(v[0])-1] = ux0_2(i * h)

AA = -2 * ( 1 / (h ** 2) + 1 / (k ** 2))

#метод Зейделя
while not end:
    eps_max = 0.
    for j_ in range(1, M):
        for i_ in range(1, N):
            v_old = v[i_][j_]
            v[i_][j_] = -f(i_*h, j_*k)
            v[i_][j_] -= v[i_ - 1][j_] * (1 / (h**2))
            v[i_][j_] -= v[i_ + 1][j_] * (1 / (h**2))
            v[i_][j_] -= v[i_ ][j_ - 1] * (1 / (k**2))
            v[i_][j_] -= v[i_ ][j_ + 1] * (1 / (k**2))
            v[i_][j_] /= AA
            eps_cur = abs(v_old - v[i_][j_])
            eps_max = max(eps_cur, eps_max)
    s += 1

    if eps_max <= eps:
        endCauseIsEps = True
        end = True

    if s >= nMax:
        endCauseIsS = True
        end = True

# ИНФОРМАЦИЯ ДЛЯ ОТЧЕТА

def fillInfo(text, data):
    for val in data:
        text = text.replace("%", str(val), 1)
    return text
def getInfo():
    ss = ""

    if endCauseIsEps:
        ss += "Выход по точности\n"
    if endCauseIsS:
        ss += "Выход по числу итераций\n"

    residual_ = []


    #считаем невязку
    for j_ in range(1, M):
        for i_ in range(1, N):
            res_i = f(i_ * h, j_ * k)
            res_i += v[i_ - 1][j_] * (1 / (h ** 2))
            res_i += v[i_ + 1][j_] * (1 / (h ** 2))
            res_i +=  v[i_][j_ - 1] * (1 / (k ** 2))
            res_i +=  v[i_][j_ + 1] * (1 / (k ** 2))
            res_i += v[i_][j_] * AA
            residual_.append(res_i)

    inf_norm_ = 0

    for ri in residual_:
        inf_norm_ += ri * ri

    inf_norm_ = inf_norm_ ** 0.5

    ss += "Норма невязки: " + str(inf_norm_) + "\n"

    ss += """При решении СЛАУ Ax=b методом Зейделя с критериями остановки Nmax=% и eps = %
     за S=% итераций достигнута точность eps_max=% и получено численное решение:\n\n"""

    ss = fillInfo(ss, [nMax, eps, s, eps_max])

    for i in range(n):
        ss += "x[" + str(i) + "] = " + str(x[i]) + "\n"

    return ss


n = (M-1)*(N-1)
x = [0.]*n

for j in range(M-1):
    for i in range(N-1):
        x[j*(N-1)+i] = v[i+1][j+1]

trueSolArr = [trueSol(x_i * h, y_i * k)  for y_i in range(1, M) for x_i in range(1, N)]
print("Истинное решение: ")
print(trueSolArr)
print()

info = getInfo()
print(info)

resTrue = np.array(x) - np.array(trueSolArr)
print("Норма разности истинного решения с V(N)")
print(max(abs(i)  for i in resTrue))