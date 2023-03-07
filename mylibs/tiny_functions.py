# Standard libraries
import numpy as np


# Умные функции
def kd_from_kp(k):
    return 2 * np.sqrt(k)


def forse_from_beam(a, diam, n, tau, b, f0, f1, f2):    # Возвращает в ССК
    rate = 1e-5
    if (f0 > -1) and (f0 < 1):
        a1 = a - f0 * n / 2
    else:
        a1 = a - np.sign(f0) * n / 2
        # tmp = rate * / (np.linalg.norm(a1) - np.linalg.norm(n)/2)**5  # **5  # np.sqrt
    tmp = rate / (np.linalg.norm(a1) - diam) ** 5  # np.sqrt
    return a1 / np.linalg.norm(a1) * tmp


def get_C_hkw(r, v, w):
    """Возвращает константы C[0]..C[5] движения Хилла-Клохесси-Уилтштира"""
    return np.double([2 * r[2] + v[0] / w,
                      v[2] / w,
                      -3 * r[2] - 2 * v[0] / w,
                      r[0] - 2 * v[2] / w,
                      v[1] / w,
                      r[1]])


def r_HKW(C, mu, w, t):
    """Возвращает вектор координат в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    r = np.zeros(3)
    r[0] = -3 * C[0] * w * t + 2 * C[1] * np.cos(w * t) - 2 * C[2] * np.sin(w * t) + C[3]
    r[1] = C[5] * np.cos(w * t) + C[4] * np.sin(w * t)
    r[2] = 2 * C[0] + C[2] * np.cos(w * t) + C[1] * np.sin(w * t)
    return np.double(r)


def v_HKW(C, mu, w, t):
    """Возвращает вектор скоростей в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    v = np.zeros(3)
    v[0] = -3 * C[0] * w - 2 * w * C[1] * np.sin(w * t) - 2 * w * C[2] * np.cos(w * t)
    v[1] = w * C[4] * np.cos(w * t) - w * C[5] * np.sin(w * t)
    v[2] = -w * C[2] * np.sin(w * t) + w * C[1] * np.cos(w * t)
    return np.double(v)


def quart2dcm(L):
    """Функция ищет матрицу поворота из кватерниона поворота; \n
    Кватернион L передаётся вектором длины 4; \n
    Возвращает матрицу 3х3."""
    w, x, y, z = L
    A = np.double(np.eye(3))
    A[0][0] = np.double(1 - 2 * y ** 2 - 2 * z ** 2)
    A[0][1] = np.double(2 * x * y + 2 * z * w)
    A[0][2] = np.double(2 * x * z - 2 * y * w)
    A[1][0] = np.double(2 * x * y - 2 * z * w)
    A[1][1] = np.double(1 - 2 * x ** 2 - 2 * z ** 2)
    A[1][2] = np.double(2 * y * z + 2 * x * w)
    A[2][0] = np.double(2 * x * z + 2 * y * w)
    A[2][1] = np.double(2 * y * z - 2 * x * w)
    A[2][2] = np.double(1 - 2 * x ** 2 - 2 * y ** 2)
    return A


def q_dot(L1, L2):
    """Функция является кватернионным умножением; \n
    Кватернион L1,L2 передаются векторами длины 4; \n
    Возвращает кватернион L[0]..L[3]."""
    return np.double(np.array([L1[0]*L2[0] - L1[1]*L2[1] - L1[2]*L2[2] - L1[3]*L2[3],
                               L1[0]*L2[1] + L1[1]*L2[0] + L1[2]*L2[3] - L1[3]*L2[2],
                               L1[0]*L2[2] + L1[2]*L2[0] + L1[3]*L2[1] - L1[1]*L2[3],
                               L1[0]*L2[3] + L1[3]*L2[0] + L1[1]*L2[2] - L1[2]*L2[1]]))


# Не такие умные функции
def clip(a, bot, top):
    if a < bot:
        return bot
    if a > top:
        return top
    return a


def my_cross(a, b):
    """Функция векторного произведения"""
    return np.double(np.array([a[1]*b[2] - a[2]*b[1],
                               a[2]*b[0] - a[0]*b[2],
                               a[0]*b[1] - a[1]*b[0]]))


def flatten(lst):
    """Функция берёт 2D массив, делает 1D"""
    return [item for sublist in lst for item in sublist]


def get_v(v, phi, theta):
    return np.array([v*np.cos(phi)*np.cos(theta),
                     v*np.sin(phi)*np.cos(theta),
                     v*np.sin(theta)])


def kronoker(a, b, tolerance=1e-6):
    """Функция является функцией кронокера;"""
    tmp = abs(np.linalg.norm(a - np.array(b)))
    return 1 if tmp < tolerance else 0

def print_time(t0, simple=False):
    if simple:
        t = t0
    else:
        t = t0.seconds
    s = t % 60
    t = int(t/60)
    m = t % 60
    h = int(t/60)
    if h > 0:
        return f"{h} час{okonchanye(h)}, {m} минут, {s} сенунд"
    else:
        if m > 0:
            return f"{m} минут, {s} сенунд"
        else:
            return f"{s} сенунд"
