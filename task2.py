nMax = 100 # максимальное количество итераций
s = 0 # количество итераций
eps = 1e-8 # точность

endCauseIsEps = False # параметры для определения причины завершения работы алгоритма
endCauseIsS = False
end = False

A = 0 # левая граница по x
B = 1 # правая граница по x
C = 0 # левая граница по y
D = 0.2 # правая граница по y
M = 8 # число разбиений по y
N = 8 # число разбиений по x
h = (B-A)/N # шаг по x
k = (D-C)/M # шаг по y
minimum_M = 8 # минимальный размер сетки по y для сохранения общего вида сетки
minimum_N = 8 # минимальный размер сетки по x для сохранения общего вида сетки
M_ratio = M / minimum_M # отношение реального M к минимальному
N_ratio = N / minimum_N # отношение реального N к минимальному
eps_max = 0. # максимальная точность, в конце работы алгоритма равна максимальной точности на последней итерации

# определяем принадлежность Vij сетке и границе
def isMainV(i, j):
    if i > 0 and i < 4 * N_ratio and j > 0 and j < 6 * M_ratio:
        return True
    if i >= 4 * N_ratio and i < 8 * N_ratio and j > 0 and j < 2 * M_ratio:
        return True
    if i > 6 * N_ratio and i < 8 * N_ratio and j >= 2 * M_ratio and  j <= 6 * M_ratio:
        return True
    if i > 2 * N_ratio and i < 8 * N_ratio and j > 6 * M_ratio and j < 8 * M_ratio:
        return True
    if i > 2 * N_ratio and i < 4 * N_ratio and j > 5 * M_ratio and j < 7 * M_ratio:
        return True
    return False

# задаем функцию f
def f(x_, y_):
    return ...

#Создаем нулевую матрицу V
v = [ ([0.]*(M+1)).copy() for i in range(N + 1)]

#Заполняем граничные значения в матрице v
#.........

AA = -2 * ( 1 / (h ** 2) + 1 / (k ** 2))

#метод Зейделя
while not end:
    eps_max = 0.
    for j_ in range(1, M):
        for i_ in range(1, N):
            if isMainV(i_, j_):
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

    # выход по точности
    if eps_max <= eps:
        endCauseIsEps = True
        end = True

    # выход по числу итераций
    if s >= nMax:
        endCauseIsS = True
        end = True
