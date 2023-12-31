import numpy as np

nMax = 10000
s = 0
eps = 1e-12
endCauseIsEps = False
endCauseIsS = False
end = False

A = 0
B = 1
C = 0
D = 0.2
#M = int(input("Введите число разбиений по y: "))
#N = int(input("Введите число разбиений по x: "))

M = 0
N = 0

def trueSol(x_, y_):
    return x_**3 + y_**3 + 2
def u(x_, y_):
    return trueSol(x_, y_)
def f(x_, y_):
    return -6 * (x_ + y_)

def print_v():
    temp_v = np.array([[0.]*(N+1)]*(M+1))
    temp_true_sol_arr = np.array([[0.]*(N+1)]*(M+1))
    for j in range(M,-1,-1):
        for i in range(0, N + 1):
            temp_v[M - j, i] = v[i][j]
            temp_true_sol_arr[M - j, i] = trueSol(i * h, j * k)
    print("V(N) в виде сетки:")
    print(temp_v)
    print("Истинное решение в виде сетки")
    print(temp_true_sol_arr)

#Заполняем граничные значения

file = open("grid.txt")

v = []

#for j in range(M + 1):
#    v.append(file.readline().split())

while True:
    line = file.readline()
    if not line:
        break
    splitted = line.split()
    N = len(splitted)
    M += 1
    v.append(splitted)

N -= 1
M -= 1

h = (B-A)/N
k = (D-C)/M


new_v = [(['a']*(M+1)).copy() for i in range(N + 1)]

for i in range(M + 1):
    for j in range(N + 1):
          new_v[j][M - i] = v[i][j]

v = new_v

print(np.array(v))

for i in range(N + 1):
    for j in range(M + 1):
        if v[i][j] == '0':
            v[i][j] = u(i * h, j * k)
        else:
            v[i][j] = 0.

print(np.array(v))

AA = -2 * ( 1 / (h ** 2) + 1 / (k ** 2))
kh = (1 / (h**2))
kk = (1 / (k**2))

#метод Зейделя
while not end:
    eps_max = 0.
    for j_ in range(1, M):
        for i_ in range(1, N):
            v_old = v[i_][j_]
            v[i_][j_] = -f(i_*h, j_*k)
            v[i_][j_] -= v[i_ - 1][j_] * kh
            v[i_][j_] -= v[i_ + 1][j_] * kh
            v[i_][j_] -= v[i_ ][j_ - 1] * kk
            v[i_][j_] -= v[i_ ][j_ + 1] * kk
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




print_v()


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
            res_i += v[i_][j_ - 1] * (1 / (k ** 2))
            res_i += v[i_][j_ + 1] * (1 / (k ** 2))
            res_i += v[i_][j_] * AA
            residual_.append(res_i)

    inf_norm_ = 0

    for ri in residual_:
        inf_norm_ += ri * ri

    inf_norm_ = inf_norm_ ** 0.5

    ss += "Норма невязки на выходе: " + str(inf_norm_) + "\n"

    residual_ = [abs(r) for r in residual_]

    ss += "Max|U-V|: " + str(max(residual_)) + "\n"

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
info = getInfo()
print(info)

resTrue = np.array(x) - np.array(trueSolArr)
print("Норма общей погрешности с V(N)")
print(max([abs(i)  for i in resTrue]))