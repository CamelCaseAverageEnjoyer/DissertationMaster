# Standard libraries
from datetime import datetime
from p_tqdm import p_map
import scipy

# Local libraries
from mylibs.tiny_functions import *
from mylibs.plot_functions import *
from mylibs.im_sample import *


def call_crash_internal_func(r, r1, r2, diam, return_force=False):
    """ Дополнительная функция для функции call_crash; \n
    Проверяет наличие точки r в цилиндре с концами r1,r2, диаметром diam; \n
    Возвращает {True,False} соответственно при отсутствии и наличии. """
    r1 = np.array(r1.copy())
    r2 = np.array(r2.copy())
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    n = np.array([x2 - x1, y2 - y1, z2 - z1])  # Вектор вдоль стержня
    length = np.linalg.norm(n)
    tau = my_cross(np.array([0, 0, 1]), n)
    if np.linalg.norm(tau) < 1e-6:
        tau = my_cross(np.array([1, 0, 0]), n)
        if np.linalg.norm(tau) < 1e-6:
            tau = my_cross(np.array([0, 1, 0]), n)
    b = my_cross(tau, n)
    a = r - (r1 + r2) / 2
    f0 = np.dot(a, n) / (np.linalg.norm(n)**2 / 2)
    f1 = np.dot(a, tau) / (np.linalg.norm(tau) * diam)
    f2 = np.dot(a, b) / (np.linalg.norm(b) * diam)

    if return_force:
        return forse_from_beam(a, diam, n, tau, b, f0, f1, f2)
    if not ((f0 > -1) and (f0 < 1) and (f1**2 + f2**2 < 1)):
        return False
    return True


def call_crash(o, r_sat, R, S, taken_beams=np.array([])):
    """ Функция проверяет наличие тела внутри балки (назовём это соударением);      \n
    Input:                                                                          \n
    -> o - AllObjects класс;                                                        \n
    -> r_sat - радиус-вектор центра аппарата; (отдельно для эпизодов)               \n
    -> R - вектор до центра масс конструкции; (отдельно для эпизодов)               \n
    -> S - матрица поворота ОСК -> ССК;       (отдельно для эпизодов)               \n
    -> taken_beams - вектор id стержней, которых учитывать не надо                  \n
    Output:                                                                         \n
    -> True в случае соударения, False иначе. """
    r = o.o_b(r_sat, S=S, R=R)
    for i in range(o.N_beams):
        if not(np.any(taken_beams == i)):
            if np.sum(o.X.flag[i]) > 0:
                r1 = o.X.r1[i]
                r2 = o.X.r2[i]
            else:
                r1 = [o.X.r_st[i][0], o.X.r_st[i][1], o.X.r_st[i][2]]
                r2 = [o.X.r_st[i][0]-o.X.length[i], o.X.r_st[i][1], o.X.r_st[i][2]]
            if call_crash_internal_func(r, r1, r2, o.d_crash):
                return True

    for i in range(o.N_cont_beams):  # Пробег по каждому элементу грузового отсека
        r1 = o.X_cont.r1[i]
        r2 = o.X_cont.r2[i]
        if call_crash_internal_func(r, r1, r2, o.X_cont.diam[i]*1.05):
            return True
    return False


def call_possible_transport(o):
    """ Функция осуществляет последовательность сборки;                                             \n
    Input:                                                                                          \n
    -> o - AllObjects класс;                                                                        \n
    Output:                                                                                         \n
    -> Отсортированный по необходимости массив id стержней, которые можно брать с грузового отсека. \n
    В большинстве случаев нужен нулевой элемент возвращаемого массива"""
    beams_to_take = np.array([])

    mask_half_fixed = np.zeros(o.N_beams)
    mask_full_fixed = np.zeros(o.N_beams)
    mask_non_fixed = np.zeros(o.N_beams)
    for i in range(o.N_beams):
        if np.sum(o.X.flag[i]) == 2:
            mask_full_fixed[i] += 1
        if np.sum(o.X.flag[i]) == 1:
            mask_half_fixed[i] += 1
        if np.sum(o.X.flag[i]) == 0:
            mask_non_fixed[i] += 1

    needed_number_nodes = np.zeros(o.N_nodes)
    current_number_nodes = np.zeros(o.N_nodes)
    mask_open_nodes = np.zeros(o.N_nodes)
    needed_number_nodes[o.X.id_node[0][0]] = 1      # Костыль на точку стыковки коснтрукции
    current_number_nodes[o.X.id_node[0][0]] = 1     # с грузовым контейнером

    for i in range(o.N_beams):
        needed_number_nodes[o.X.id_node[i][0]] += 1  # Сколько стержней приходят в узел
        needed_number_nodes[o.X.id_node[i][1]] += 1
        current_number_nodes[o.X.id_node[i][0]] += o.X.flag[i][0]  # Сколько стержней в узле находятся
        current_number_nodes[o.X.id_node[i][1]] += o.X.flag[i][1]

    for i in range(o.N_nodes):  # В каких узлах неполное кол-во стержней, но есть хоть один
        if (needed_number_nodes[i] - current_number_nodes[i] > 0) and (current_number_nodes[i] > 0):
            mask_open_nodes[i] = 1  # Основная маска

    for i in range(o.N_beams):
        if mask_non_fixed[i] > 0:  # берём нетронутые со склада
            if mask_open_nodes[o.X.id_node[i][0]] + mask_open_nodes[o.X.id_node[i][1]] > 0:  # проверка надобности балки
                beams_to_take = np.append(beams_to_take, o.X.id[i])

    i = 0
    while i < len(beams_to_take):  # Удалить те, которые уже взяты
        if beams_to_take[i] in o.taken_beams:
            beams_to_take = np.delete(beams_to_take, i)
            i -= 1
        i += 1

    return [int(i) for i in beams_to_take]


def call_inertia(o, id_s=np.array([]), app_y=None, app_n=None):
    """Функция считает тензор инерции в собственной ск конструкции и центр масс;    \n
    Input:                                                                          \n
    -> o - AllObjects класс;                                                        \n
    -> id_s - вектор id стержней, которых учитывать не надо                         \n
    -> app_y - номера аппаратов, которых стоит учесть при расчёте инерции;          \n
    -> app_n - номера аппаратов, которых не стоит учитывать при расчёте инерции;    \n
    Output:                                                                         \n
    -> Тензор инерции, вектор центра масс. """
    J = np.double(np.zeros((3, 3)))
    N_points_m = 10                                                 # количество материальных точек на стержень
    N_m = N_points_m * (o.N_beams + o.N_cont_beams)                 # материальных точек всего
    m = np.double([0. for i in range(N_m + o.N_app)])
    r = np.double([[0., 0., 0.] for i in range(N_m + o.N_app)])

    # Выделение координат стержня - на месте или в грузовом отсеке
    for n in range(o.N_beams):
        if not np.any(id_s == n):
            if np.sum(o.X.flag[n]) > 0:
                r_1 = np.array(o.X.r1[n])
                r_2 = np.array(o.X.r2[n])
            else:
                r_1 = o.X.r_st[n]
                r_2 = o.X.r_st[n] - np.array([o.X.length[n], 0, 0])

            for i in range(N_points_m):
                m[n * N_points_m + i] = o.X.mass[n] / N_points_m
                r[n * N_points_m + i] = r_1 + (r_2 - r_1) * i / N_points_m

    # Учёт грузового отсека на тензор инерции
    for n in range(o.N_cont_beams):
        r_1 = np.array(o.X_cont.r1[n])
        r_2 = np.array(o.X_cont.r2[n])
        for i in range(N_points_m):
            m[(o.N_beams + n) * N_points_m + i] = o.X_cont.mass[n] / N_points_m
            if n > 0:
                r[(o.N_beams + n) * N_points_m + i] = r_1 + (r_2 - r_1) * i / N_points_m
            else:
                r_0 = (r_1 + r_2) / 2
                b = o.X_cont.diam[n] * my_cross((r_2 - r_1), [0, 1, 0]) / np.linalg.norm(r_2 - r_1)
                k = o.X_cont.diam[n] * my_cross((r_2 - r_1), b) / np.linalg.norm(r_2 - r_1) / np.linalg.norm(b)
                r[(o.N_beams + n) * N_points_m + i] = r_0 + b * np.sin(2*np.pi * i / N_points_m) + \
                                                      k * np.cos(2*np.pi * i / N_points_m)

    # Учёт аппаратов
    tmp_i = 0
    for a in range(o.N_app):
        if (not o.X_app.flag_fly[a] and app_n != a) or app_y == a:
            m[N_m + tmp_i] = o.X_app.mass[a]
        else:
            m[N_m + tmp_i] = 0
            r[N_m + tmp_i] = o.X_app.target[a]
        tmp_i += 1

    # Вычисление центра масс
    tmp_numerator = np.double(np.zeros(3))
    tmp_denominator = np.double(np.zeros(3))
    for k in range(N_m):
        tmp_numerator += m[k] * np.array(r[k])
        tmp_denominator += m[k]
    r_mass_center = tmp_numerator / tmp_denominator
    for k in range(N_m):
        r[k] -= r_mass_center  # С заботой о ваших нервах J считается относительно центра масс

    for k in range(N_m):
        for i in range(3):
            for j in range(3):
                J[i][j] += m[k] * (
                        kronoker(i, j) * (r[k][0] ** 2 + r[k][1] ** 2 + r[k][2] ** 2) - r[k][i] * r[k][j])
    return J, np.array(r_mass_center)


def call_middle_point(X_cont, r_right):
    """ Функция возвращает id удобной ручки грузового отсека как промежуточная цель;    \n
    Input:                                                                              \n
    -> информационная матрица X_cont, содержашая столбцы:                               \n
    id | mass | length | diam | r1 | r2 | flag_grab                                     \n
    -> r_right - радиус-вектор таргетной точки на конструкции ССК                       \n
    Output:                                                                             \n
    -> id-номер удобной ручки крепления с конструкцией """
    X_middles = X_cont[X_cont.flag_grab]  # Стержни - узлы захвата
    difference = [np.array(X_middles.r1[i]) / 2 + np.array(X_middles.r2[i]) / 2 for i in X_middles.id]
    N_middles = len(difference)
    difs = np.zeros(N_middles)
    for i in range(N_middles):
        difs[i] = np.linalg.norm(np.array(difference[i]) - r_right)
    i_needed = np.argmin(difs)

    return np.array(X_middles.id)[i_needed]


def calculation_motion(o, u0, t, T_max, id_app, interaction=True, line_return=False):
    """Функция расчитывает угловое и поступательное движение одного аппарата и конструкции; \n
    Других аппаратов и их решения она не учитывает;                                         \n
    Input:                                                                                  \n
    -> o - AllObjects класс;                                                                \n
    -> u0 - начальная скорость аппарата (ССК при interaction=True, ОСК иначе)               \n
    -> t - время отталкивания/импульса (для расчёта положения/скорости ХКУ)                 \n
    -> T_max - время эпизода                                                                \n
    -> id_app - номер аппарата                                                              \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции        \n
    -> line_return - возвращать ли линию полёта аппарата ССК (для красивых картинок)        \n
    Output:                                                                                     \n
    -> F_min - минимальная невязка-вектор от нужной точки до аппарата                           \n
    -> F_min - средняя невязка-вектор от нужной точки до аппарата                               \n
    -> w_modeling/w_after - максимальный модуль угловой скорости станции по ходу/после эпизода  \n
    -> V_max - максимальный модуль поступательной скорости станции по ходу/после эпизода        \n
    -> R_max - максимальный модуль положения ОСК станции по ходу эпизода                        \n
    -> j_max - максимальный угол отклонения от целевого положения станции по ходу эпизода       \n
    -> N_crashes - количество соударений (не должно превосходить 1)                             \n
    -> line - линия аппарата ССК (при line_return=True)"""
    La = o.La.copy()
    u = u0.copy()
    V = np.array(v_HKW(o.C_R, o.mu, o.w_hkw, t - o.t_start[o.N_app]))
    R = np.array(r_HKW(o.C_R, o.mu, o.w_hkw, t - o.t_start[o.N_app]))
    w = np.array(o.w.copy())

    r = np.array(o.X_app.r[id_app])
    r0 = o.o_b(r)
    r_right = np.array(o.X_app.target[id_app])
    id_beam = o.X_app.flag_beam[id_app]
    m_app = o.X_app.mass[id_app]

    # Initialization
    m_beam = 0. if (id_beam is None) else o.X.mass[id_beam]
    M = np.sum(o.X.mass.to_numpy()) + np.sum(o.X_cont.mass.to_numpy())
    M_without = M - m_beam
    m_extra = m_app + m_beam
    w_hkw = np.linalg.norm(o.W_hkw)

    J_p, r_center_p = call_inertia(o, o.taken_beams_p, app_y=id_app)
    J, r_center = call_inertia(o, o.taken_beams, app_n=id_app)
    J_1 = np.linalg.inv(J)
    U, S, A, R_e = o.call_rotation_matrix(La, t)

    if interaction:
        R_p = R.copy()
        V_p = V.copy()
        R = R_p + S.T @ (r_center - r_center_p)                 # ORF
        u_rot = my_cross(w, S.T @ r0)                         # ORF
        V_rot = my_cross(w, S.T @ (r_center - r_center_p))    # ORF
        V0 = - u * m_extra / M_without                         # BRF
        u = V_p + S.T @ u + u_rot                              # ORF
        V = V_p + S.T @ V0 + V_rot                              # ORF
        w = o.b_o(J_1, S) @ (                                   # ORF
            o.b_o(J_p, S) @ w - m_extra * my_cross(r, u) +
            (M_without + m_extra) * my_cross(R_p, V_p) - M_without * my_cross(R, V)) 

    Om = U.T @ w + o.W_hkw
    C_r = get_C_hkw(r, u, w_hkw)
    C_R = get_C_hkw(R, V, w_hkw)

    F_min = np.array([1e10, 0, 0]) 
    w_modeling = 0
    N_crashes = 0
    F_mean = 0. 
    line = r0                       
    V_max = 0.
    R_max = 0.
    j_max = 0.
    dt = o.dt
    
    # Iterations
    for i in range(round(T_max / dt)):
        t_i = (i + 1) * dt
        r = np.array(r_HKW(C_r, o.mu, o.w_hkw, t_i))
        R = np.array(r_HKW(C_R, o.mu, o.w_hkw, t_i))
        V = np.array(v_HKW(C_R, o.mu, o.w_hkw, t_i))

        La, Om = o.rk4_w(La, Om, J, t + t_i)
        U, S, A, R_e = o.call_rotation_matrix(La, t + t_i)
        w = U @ (Om - o.W_hkw)

        # Writing parameters
        if interaction:
            F = o.o_b(r, S=S, R=R, r_center=r_center) - r_right  # BRF
        else:
            F = r - o.b_o(r_right, S=S, R=R, r_center=r_center)  # ORF
        if np.linalg.norm(F_min) > np.linalg.norm(F):
            F_min = F.copy()
        w_modeling = max(w_modeling, np.linalg.norm(w))
        V_max = max(V_max, np.linalg.norm(V))
        R_max = max(R_max, np.linalg.norm(R))
        r_max = max(R_max, np.linalg.norm(r))
        j_max = max(j_max, 180/np.pi*np.arccos((np.trace(S)-1)/2))
        F_mean += np.linalg.norm(F)
        iter = i
        line = np.append(line, S @ (r - R) + r_center)

        # Stop Criteria
        if o.d_to_grab is not None:
            if np.linalg.norm(F) <= o.d_to_grab * 0.95:
                break

        if o.d_crash is not None:
            if t_i > 20:
                if call_crash(o, r, R, S, o.taken_beams):
                    N_crashes = N_crashes + 1
                    break

    if i > 0:
        F_mean /= i

    if interaction is False:
        V_p = np.array(v_HKW(C_R, o.mu, o.w_hkw, t_i))
        u_p = np.array(v_HKW(C_r, o.mu, o.w_hkw, t_i))
        J_p, r_center_p = call_inertia(o, o.taken_beams_p, app_n=id_app)
        J, r_center = call_inertia(o, o.taken_beams, app_y=id_app)
        R_p = R.copy()
        R += S.T @ (r_center - r_center_p)
        V = (V_p * M_without + u_p * m_extra) / (M_without + m_extra)
        w_after = np.linalg.inv(o.b_o(J, S)) @ (o.b_o(J_p, S) @ w - (M_without + m_extra) * my_cross(R, V) +
                                   m_extra * my_cross(r, u_p) + M_without * my_cross(R_p, V_p))
        w_modeling = np.linalg.norm(w_after)
        Om = U.T @ w_after + o.W_hkw
        C_R = get_C_hkw(R, V, o.w_hkw)
        V_max = np.linalg.norm(V)
        R_max = np.linalg.norm(R)

        for i in range(round(iter / 20)):
            t_i += dt
            R = np.array(r_HKW(C_R, o.mu, o.w_hkw, t_i))
            V = np.array(v_HKW(C_R, o.mu, o.w_hkw, t_i))

            La, Om = o.rk4_w(La, Om, J, t + t_i)
            U, S, A, R_e = o.call_rotation_matrix(La, t + t_i)
            w = U @ (Om - o.W_hkw)

            w_modeling = max(w_modeling, np.linalg.norm(w))
            V_max = max(V_max, np.linalg.norm(V))
            R_max = max(R_max, np.linalg.norm(R))
            r_max = max(R_max, np.linalg.norm(r))
            j_max = max(j_max, 180/np.pi*np.arccos((np.trace(S)-1)/2))

    w = w_modeling if interaction else  np.linalg.norm(w_after)
    if line_return:
        return F_min, F_mean, w, V_max, R_max, j_max, N_crashes, line
    else:
        return F_min, F_mean, w, V_max, R_max, j_max, N_crashes


def diff_evolve_sample(o, t, T_max, id_app, interaction, convex, tmp1x, tmp1y, tmp1z, tmp2x, tmp2y, tmp2z,
                       tmp3x, tmp3y, tmp3z, ajx, ajy, ajz, a_recjx, a_recjy, a_recjz, iter):
    """- Зачем я существую?         \n
    - Ты - условная единица парпрога"""
    tmp1 = np.array([tmp1x, tmp1y, tmp1z])
    tmp2 = np.array([tmp2x, tmp2y, tmp2z])
    tmp3 = np.array([tmp3x, tmp3y, tmp3z])
    aj = np.array([ajx, ajy, ajz])
    a_recj = np.array([a_recjx, a_recjy, a_recjz])
    u_p = o.X_app.v[id_app]
    mu_IPM = o.mu_IPM / 2**iter

    tmp = tmp1 + o.diff_evolve_F * (tmp2 - tmp3)  # mutant vector
    if random.uniform(0, 1) < o.diff_evolve_chance:  # crossing
        tmp[0] = aj[0]
    if random.uniform(0, 1) < o.diff_evolve_chance:
        tmp[1] = aj[1]
    if random.uniform(0, 1) < o.diff_evolve_chance:
        tmp[2] = aj[2]
    u_pre = np.array([aj[0] * np.cos(aj[1]) * np.cos(aj[2]),
                      aj[0] * np.sin(aj[1]) * np.cos(aj[2]),
                      aj[0] * np.sin(aj[2])])
    u_new = np.array([tmp[0] * np.cos(tmp[1]) * np.cos(tmp[2]),
                      tmp[0] * np.sin(tmp[1]) * np.cos(tmp[2]),
                      tmp[0] * np.sin(tmp[2])])
    u_pre = o.cases['diff_vel_control'](u_pre, np.linalg.norm(u_p) == 0.)
    u_new = o.cases['diff_vel_control'](u_new, np.linalg.norm(u_p) == 0.)

    if a_recj[0] is None:
        if convex:
            dr_p, _, _ = dr_gradient(o=o, u0=u_pre, t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=0.1)
        else:
            dr_p = dr_evolve(u_pre, o, t, T_max, id_app, interaction)
    else:
        dr_p = a_recj

    if np.linalg.norm(u_p) == 0.:
        du = u_pre - np.array(o.X_app.v[id_app])
        du_m = np.linalg.norm(du)
        dr_p = dr_p - mu_IPM * np.log([o.du_impulse_max - du_m, o.du_impulse_max - du_m, o.du_impulse_max - du_m]) \
            if o.du_impulse_max > du_m else np.array([1e10, 1e10, 1e10])

    if convex:
        dr_n, _, real_dr = dr_gradient(o=o, u0=u_new, t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=o.mu_e)
    else:
        dr_n = dr_evolve(u_new, o, t, T_max, id_app, interaction)

    if np.linalg.norm(u_p) == 0.:
        du = u_new - np.array(o.X_app.v[id_app])
        du_m = np.linalg.norm(du)
        dr_n = dr_n - mu_IPM * np.log([o.du_impulse_max - du_m, o.du_impulse_max - du_m, o.du_impulse_max - du_m]) \
            if o.du_impulse_max > du_m else np.array([1e10, 1e10, 1e10])*clip(10*(du_m - o.du_impulse_max), 0, 1)
    a_recj = dr_n if np.linalg.norm(dr_n) < np.linalg.norm(dr_p) else dr_p
    aj = tmp if np.linalg.norm(dr_n) < np.linalg.norm(dr_p) else aj
    return [aj[0], aj[1], aj[2], a_recj[0], a_recj[1], a_recj[2]]


def diff_evolve(o, t, T_max, id_app, interaction=True, convex=False):
    """Функция дифференциальной эволюции - неградиентный поиск скорости отталкивания/импульса   \n
    Input:                                                                                      \n
    -> o - AllObjects класс;                                                                    \n
    -> t - время отталкивания/импульса (для расчёта положения/скорости ХКУ)                     \n
    -> T_max - время эпизода                                                                    \n
    -> id_app - номер аппарата                                                                  \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции            \n
    Output:                                                                                     \n
    -> u - вектор отталкивания/импульса (ССК при interaction=True, ОСК иначе)"""
    if o.if_talk:
        talk()
    n = o.diff_evolve_vectors
    lst_errors = []

    a = [[np.exp(random.uniform(np.log(o.u_min), np.log(o.u_max))), 
          random.uniform(-np.pi, np.pi), 
          random.uniform(-np.pi / 2, np.pi / 2)]
          for _ in range(n)]
    a_rec = [[None, None, None] for _ in range(n)]
    tmp = [[0. for i in range(6)] for _ in range(n)]
    tmp1_ind = [-1 for i in range(n)]
    tmp2_ind = [-1 for i in range(n)]
    tmp3_ind = [-1 for i in range(n)]
    for i in range(o.diff_evolve_times):
        if o.if_testing_mode:
            print(Style.RESET_ALL + f"Шаг {i + 1}/{o.diff_evolve_times} дифференциальной эволюции")
        for j in range(n):
            to_choice = np.delete(range(n), j)
            tmp1_ind[j] = random.choice(to_choice)
            tmp2_ind[j] = random.choice(to_choice)
            tmp3_ind[j] = random.choice(to_choice)
        if o.if_multiprocessing:
            tmp = p_map(diff_evolve_sample,
                        [o for ii in range(n)],
                        [t for ii in range(n)],
                        [T_max for ii in range(n)],
                        [id_app for ii in range(n)],
                        [interaction for ii in range(n)],
                        [convex for ii in range(n)],
                        [a[tmp1_ind[ii]][0] for ii in range(n)],
                        [a[tmp1_ind[ii]][1] for ii in range(n)],
                        [a[tmp1_ind[ii]][2] for ii in range(n)],
                        [a[tmp2_ind[ii]][0] for ii in range(n)],
                        [a[tmp2_ind[ii]][1] for ii in range(n)],
                        [a[tmp2_ind[ii]][2] for ii in range(n)],
                        [a[tmp3_ind[ii]][0] for ii in range(n)],
                        [a[tmp3_ind[ii]][1] for ii in range(n)],
                        [a[tmp3_ind[ii]][2] for ii in range(n)],
                        [a[ii][0] for ii in range(n)],
                        [a[ii][1] for ii in range(n)],
                        [a[ii][2] for ii in range(n)],
                        [a_rec[ii][0] for ii in range(n)],
                        [a_rec[ii][1] for ii in range(n)],
                        [a_rec[ii][2] for ii in range(n)],
                        [n for ii in range(n)])
        else:
            for ii in range(n):
                tmp[ii] = diff_evolve_sample(o,t,T_max,id_app,interaction,
                        a[tmp1_ind[ii]][0],a[tmp1_ind[ii]][1],a[tmp1_ind[ii]][2],
                        a[tmp2_ind[ii]][0],a[tmp2_ind[ii]][1],a[tmp2_ind[ii]][2],
                        a[tmp3_ind[ii]][0],a[tmp3_ind[ii]][1],a[tmp3_ind[ii]][2],
                        a[ii][0],a[ii][1],a[ii][2],a_rec[ii][0],a_rec[ii][1],a_rec[ii][2])
        
        a = [[tmp[ii][0],tmp[ii][1],tmp[ii][2]] for ii in range(n)]
        a_rec = [[tmp[ii][3],tmp[ii][4],tmp[ii][5]] for ii in range(n)]
        tmp2 = [np.linalg.norm(a_rec[ii]) for ii in range(n)]
        if o.if_testing_mode:
            print(Style.RESET_ALL + f"ошибки проб {tmp2}")
        u_ = a[np.argmin(tmp2)]
        lst_errors.append(np.min(tmp2))
    u = np.array([u_[0] * np.cos(u_[1]) * np.cos(u_[2]),
                  u_[0] * np.sin(u_[1]) * np.cos(u_[2]),
                  u_[0] * np.sin(u_[2])])
    u = o.cases['diff_vel_control'](u, True)
    if o.if_testing_mode:
        print(Fore.GREEN + f"Сделали скорость:{u}, прогресс целевой функции: {lst_errors}")
    return np.array(u)


def dr_gradient(o, u0, t, T_max, id_app, interaction, mu_IPM, mu_e):
    dr_, _, w_, V_, R_, j_, _ = calculation_motion(o=o, u0=u0, t=t, T_max=T_max, id_app=id_app, interaction=interaction)
    reserve_rate = 1.5
    e_w = 0. if (o.w_max/reserve_rate - w_) > 0 else abs(w_ - o.w_max/reserve_rate)
    e_V = 0. if (o.V_max/reserve_rate - V_) > 0 else abs(V_ - o.V_max/reserve_rate)
    e_R = 0. if (o.R_max/reserve_rate - R_) > 0 else abs(R_ - o.R_max/reserve_rate)
    e_j = 0. if (o.j_max/reserve_rate - j_) > 0 else abs(j_ - o.j_max/reserve_rate)
    if o.if_testing_mode:
        print(Style.RESET_ALL + f"e: {e_w/o.w_max} | {e_V/o.V_max} | {e_R/o.R_max} | {e_j/o.j_max}")
    cnd = (e_w + e_V + e_R + e_j > 0)
    return dr_ - \
           mu_IPM * (np.log(o.w_max - w_ + e_w) + np.log(o.V_max - V_ + e_V) + np.log(o.R_max - R_ + e_R) + np.log(o.j_max - j_ + e_j)) + \
           mu_e * (e_w/o.w_max + e_V/o.V_max + e_R/o.R_max + e_j/o.j_max), \
           cnd, np.linalg.norm(dr_)


def dr_evolve(u, o, t, T_max, id_app, interaction):
    mu_IPM = 0.01
    dr_p, dr_m_p, w_p, V_p, R_p, j_p, N_crashes_p = calculation_motion(o, u, t, T_max, id_app, interaction=interaction)

    # Нет проблем с ограничениями
    if (o.w_max - w_p > 0) and (o.V_max - V_p > 0) and (o.R_max - R_p > 0) and (o.j_max - j_p > 0):  
        return dr_p - mu_IPM * (np.log(o.w_max - w_p) + np.log(o.V_max - V_p) + np.log(o.R_max - R_p) + np.log(o.j_max - j_p)) + \
            np.array([1e50, 1e50, 1e50]) * N_crashes_p
    # Нарушение ограничений
    else:  
        return np.array([1e5, 1e5, 1e5]) * \
            (clip(10*(w_p - o.w_max), 0, 1) + clip(10*(V_p - o.V_max), 0, 1) + clip(1e-1*(R_p - o.R_max), 0, 1) + clip(1e-2*(j_p - o.j_max), 0, 1)) + \
            np.array([1e50, 1e50, 1e50]) * N_crashes_p


def calc_shooting(o, t, id_app, r_right, interaction=True, shooting_amount=None, T_max=None, convex_domain=True, mu_e=None):  
    """ Функция выполняет пристрелочный/спектральный поиск оптимальной скорости отталкивания/импульса аппарата; \n
    Input:                                                                                  \n
    -> o - AllObjects класс;                                                                \n
    -> t - время отталкивания/импульса (для расчёта положения/скорости ХКУ)                 \n
    -> id_app - номер аппарата                                                              \n
    -> r_right - радиус-вектор ССК положения цели                                           \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции        \n
    -> shooting_amount - количество шагов пристрелки/спектрального пучка                    \n
    -> T_max - время эпизода                                                                \n
    -> convex_domain - считать ли пространство целевой функции выпуклым                     \n
    Output:                                                                                 \n
    -> u - оптимальный вектор скорости отталкивания/импульса (ССК при interaction=True, ОСК иначе)"""
    # shooting_amount = (o.shooting_amount_repulsion if interaction else o.shooting_amount_impulse) if (shooting_amount is None) else shooting_amount
    shooting_amount = o.shooting_amount_repulsion if (shooting_amount is None) else shooting_amount
    if interaction:
        tmp = r_right - np.array(o.o_b(o.X_app.r[id_app]))
    else:
        tmp = o.b_o(r_right) - np.array(o.X_app.r[id_app])
    mu_e = o.mu_e if (mu_e is None) else mu_e
    T_max = 2*(np.linalg.norm(tmp)/o.u_min if o.X_app.flag_fly[id_app] else o.T_max)  # if (T_max is None) else T_max
    T_max_hard_limit = T_max*16 if o.X_app.flag_fly[id_app] else o.T_max_hard_limit
    u = o.u_min * tmp / np.linalg.norm(tmp)

    # Differential evolution
    start_time = datetime.now()
    # if interaction:  # это кому-то было лень делать нормально
    # u = diff_evolve(o=o, t=t, T_max=T_max, id_app=id_app, interaction=interaction, convex=not interaction)
    '''cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
            {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
            {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})
    u = sp.optimize.minimize(dr_evolve, u, o, t, T_max, id_app, interaction, method='COBYLA', constraints=cons, 
                             callback='callable')'''

    if o.if_testing_mode:
        print(Style.RESET_ALL + 'Время, потраченное на дифференциальную эволюцию:', datetime.now() - start_time)

    if convex_domain:
        # Метод пристрелки
        mu_IPM = o.mu_IPM
        # for i_iteration in range(shooting_amount):
        i_iteration = 0
        i_iteration_for_stable_T = 0
        while i_iteration < shooting_amount:
            du = 0.00001
            u_i = np.eye(3) * du
            mu_IPM /= 5
            mu_e /= 2
            i_iteration += 1
            i_iteration_for_stable_T += 1

            dr, c, tol = dr_gradient(o=o, u0=u, t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=mu_e)
            dr_x, cx, _ = dr_gradient(o=o, u0=u + u_i[:, 0], t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=mu_e)
            dr_y, cy, _ = dr_gradient(o=o, u0=u + u_i[:, 1], t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=mu_e)
            dr_z, cz, _ = dr_gradient(o=o, u0=u + u_i[:, 2], t=t, T_max=T_max, id_app=id_app, interaction=interaction, mu_IPM=mu_IPM, mu_e=mu_e)
            if o.if_any_print:
                print(Style.RESET_ALL + f"Точность пристрелки: {tol} метров, разрешённое время {T_max} секунд")
            if tol < o.d_to_grab*0.98:  # and not (c or cx or cy or cz): 
                break

            # Коррекция скорости отталкивания
            Jacobian = np.array([[dr_x[0] - dr[0], dr_y[0] - dr[0], dr_z[0] - dr[0]],
                                 [dr_x[1] - dr[1], dr_y[1] - dr[1], dr_z[1] - dr[1]],
                                 [dr_x[2] - dr[2], dr_y[2] - dr[2], dr_z[2] - dr[2]]]) / du
            Jacobian_1 = np.linalg.inv(Jacobian)
            correction = Jacobian_1 @ dr
            correction = correction/np.linalg.norm(correction)*clip(np.linalg.norm(correction), 0, o.u_max/2)
            u = u - correction

            # Ограничение по модулю скорости
            if np.linalg.norm(u) > o.u_max:
                if o.if_any_print:
                    print(Style.RESET_ALL + f'attention: shooting speed reduction: {np.linalg.norm(u)/o.u_max*100} %')
                u = u / np.linalg.norm(u) * o.u_max * 0.9

            if c or cx or cy or cz:
                mu_IPM *= 10
                i_iteration -= 1
            '''if i_iteration_for_stable_T >= 5 and ((c and cx and cy and cz) > 0 or np.linalg.norm(dr) > o.d_to_grab):  # (tolerance > o.d_to_grab) or (
                # T_max *= 2 if T_max < T_max_hard_limit else 1
                # u = o.u_max / 5 * tmp / np.linalg.norm(tmp)
                i_iteration_for_stable_T = 0
                i_iteration -= 5'''
        target_is_reached = True  # tol <= o.d_to_grab and ((c or cx or cy or cz) is False)
    else: 
        # Метод спектрального пучка
        """N_bundle = 3  # обязательно нечётный
        du = np.linalg.norm(u) / 5
        print('Метод спектрального пучка, du =', du)
        for i_iteration in range(shooting_amount):
            v1 = my_cross(u, [0, 0, 1]) / np.linalg.norm(my_cross(u, [0, 0, 1]))
            v2 = my_cross(u, v1) / np.linalg.norm(my_cross(u, v1))
            u_map = np.array([[(u + v1 * (i - (N_bundle - 1) / 2) / N_bundle * du + v2 * (
                        j - (N_bundle - 1) / 2) / N_bundle * du) for i in range(N_bundle)] for j in range(N_bundle)])
            tolerance = np.zeros((N_bundle, N_bundle))
            for i in range(N_bundle):
                for j in range(N_bundle):
                    tmp = u_map[i, j]
                    if np.linalg.norm(tmp) > u_max:
                        tmp = tmp / np.linalg.norm(tmp) * u_max
                    calc_args.u = np.array(tmp)
                    dr, w_modeling_0, N_crashes = calculation_motion(calc_args)
                    tolerance[i, j] = np.linalg.norm(dr)
            i_best = int(np.argmin(tolerance) / N_bundle)
            j_best = np.argmin(tolerance) % N_bundle
            print('Точность', tolerance[i_best, j_best])
            u = np.array(u_map[i_best, j_best])
            du = du / N_bundle"""
    return u, target_is_reached


def repulsion(o, t, id_app, u_a_priori=None):
    """Input:                                                                   \n
    -> o - AllObjects класс                                                     \n
    -> t - время принятия решения (для расчёта положения и скорости ХКУ)        \n
    -> id_app - номер аппарата                                                  \n
    -> u_a_priori - заданный вектор скорости (на случай вынужденного решения,   \n
    используется для быстрой отладки алгоритмов управления"""

    # Параметры до отталкивания
    J_p, r_center_p = call_inertia(o, o.taken_beams, app_y=id_app)
    o.taken_beams_p = o.taken_beams.copy()

    # Алгоритм выбора цели
    if o.X_app.flag_start[id_app]:  # IF ON START
        id_beam = int(np.round(call_possible_transport(o)[0]))
        o.taken_beams = np.append(o.taken_beams, id_beam)
        o.my_print(f"Аппарат {id_app} забрал стержень {id_beam}", mode="blue")
    else:
        if o.X_app.flag_beam[id_app] is None:   # GOING HOME
            id_beam = None
            o.my_print(f"Аппарат {id_app} возвращается на базу", mode="blue")
        else:                                   # GOING TO POINT
            id_beam = o.X_app.flag_beam[id_app]
            o.my_print(f"Аппарат {id_app} летит установливать стержень {id_beam}", mode="blue")

    if id_beam is not None:
        r_right = np.array(o.X.r1[id_beam])
    else:
        tmp_id = call_possible_transport(o)[0]
        r_right = np.array(o.X.r_st[tmp_id] - np.array([o.X.length[tmp_id], 0.0, 0.0]))

    if o.X_app.flag_start[id_app] or (not o.X_app.flag_to_mid[id_app]):
        m = call_middle_point(o.X_cont, r_right)
        if np.linalg.norm(o.X_app.target[id_app] - np.array(o.X_cont.r1[m]) / 2 - np.array(o.X_cont.r2[m]) / 2) > 1e-2:
            r_right = np.array(o.X_cont.r1[m]) / 2 + np.array(o.X_cont.r2[m]) / 2
        o.X_app.loc[id_app, 'flag_to_mid'] = True
        o.my_print(f"Аппарат {id_app} целится в промежуточный узел захвата", mode="blue")
    else:
        o.X_app.loc[id_app, 'flag_to_mid'] = False

    o.X_app.loc[id_app, 'target_p'][0], o.X_app.loc[id_app, 'target_p'][1], o.X_app.loc[id_app, 'target_p'][2] = \
        o.X_app.target[id_app]
    o.X_app.loc[id_app, 'target'][0], o.X_app.loc[id_app, 'target'][1], o.X_app.loc[id_app, 'target'][2] = r_right
    o.X_app.loc[id_app, 'flag_beam'] = id_beam
    o.X_app.loc[id_app, 'flag_start'] = False

    # Нахождение вектора скорости отталкивания
    if o.control:
        u0 = u_a_priori if (u_a_priori is not None) else \
            diff_evolve(o=o, t=t, T_max=o.T_max, id_app=id_app, interaction=True)
    else:
        u0, _ = (u_a_priori, 1) if (u_a_priori is not None) else \
            calc_shooting(o=o, t=t, id_app=id_app, r_right=r_right, interaction=True,
                          shooting_amount=o.shooting_amount_repulsion, T_max=o.T_max)

    R_p = o.R.copy()
    V_p = o.V.copy()
    o.J, o.r_center = call_inertia(o, o.taken_beams, app_n=id_app)
    o.J_1 = np.linalg.inv(o.J)
    m_beam = 0. if (id_beam is None) else o.X.mass[id_beam]
    M_without = np.sum(o.X.mass.to_numpy()) + np.sum(o.X_cont.mass.to_numpy()) - m_beam
    m_extra = o.X_app.mass[id_app] + m_beam
    R = R_p + o.S.T @ (o.r_center - r_center_p)
    r = np.array(o.X_app.r[id_app])
    r0 = o.o_b(r, S=o.S, R=R_p, r_center=r_center_p) 

    u_rot = my_cross(o.w, o.S.T @ r0)                         # ORF
    V_rot = my_cross(o.w, o.S.T @ (o.r_center - r_center_p))  # ORF
    V0 = - u0 * m_extra / M_without                           # BRF
    u = V_p + o.S.T @ u0 + u_rot                              # ORF
    V = V_p + o.S.T @ V0 + V_rot                              # ORF
    o.w = o.b_o(o.J_1) @ (
            o.b_o(J_p) @ o.w - m_extra * my_cross(r, u) +
            (M_without + m_extra) * my_cross(R_p, V_p) - M_without * my_cross(R, V))  # ORF
    o.Om_update()
    o.C_r[id_app] = get_C_hkw(r, u, o.w_hkw)
    o.C_R = get_C_hkw(R, V, o.w_hkw)
    o.flag_vision[id_app] = False
    o.t_start[0] = t
    o.t_start[o.N_app] = t

    o.X_app.loc[id_app, 'r'][0], o.X_app.loc[id_app, 'r'][1], o.X_app.loc[id_app, 'r'][2] = r
    o.X_app.loc[id_app, 'v'][0], o.X_app.loc[id_app, 'v'][1], o.X_app.loc[id_app, 'v'][2] = u
    o.X_app.loc[id_app, 'r_0'] = o.get_discrepancy(id_app)
    o.X_app.loc[id_app, 'flag_fly'] = True
    o.X_app.loc[id_app, 'flag_hkw'] = False if (o.if_PID_control or o.if_LQR_control) else True
    o.my_print(f"App id:{id_app} pushed off with u={np.linalg.norm(u)}, w={np.linalg.norm(o.w)}", mode="blue", test=True)
    o.my_print(f"Taken beams list: {o.taken_beams}", mode="green", test=True)
    talk_decision(o.if_talk)

    return u0


def capturing(o, id_app, t):
    if o.if_any_print:
        print(Fore.BLUE + f'Аппарат id:{id_app} захватился')
    o.a_self[id_app] = np.array(np.zeros(3))
    talk_success(o.if_talk)

    # Параметры до захвата
    id_beam = o.X_app.flag_beam[id_app]
    m_beam = 0 if id_beam is None else o.X.mass[id_beam]
    M = np.sum(o.X.mass.to_numpy()) + np.sum(o.X_cont.mass.to_numpy())
    M_without = M - m_beam
    m_extra = o.X_app.mass[id_app] + m_beam
    r = np.array(o.X_app.r[id_app])
    u_p = np.array(o.X_app.v[id_app])
    V_p = o.V.copy()
    R_p = o.R.copy()
    J_p, r_center_p = call_inertia(o, o.taken_beams, app_n=id_app)

    # Процесс захвата
    if id_beam is not None:
        if np.linalg.norm(np.array(o.X.r1[id_beam]) - np.array(o.X_app.target[0])) < 1e-2:
            o.my_print(f"Стержень id:{id_beam} устанавливается", mode="blue")
            o.X.loc[id_beam, 'flag'][0] = 1
            o.X.loc[id_beam, 'flag'][1] = 1
            o.X_app.loc[id_app, 'flag_beam'] = None
            o.taken_beams = np.delete(o.taken_beams, np.argmax(o.taken_beams == id_beam))
    else:
        if o.X_app.target[id_app][0] < -0.6:  # Если "слева" нет промежуточных точек, то окей
            o.my_print(f'Аппарат id:{id_app} в грузовом отсеке')
            o.X_app.loc[id_app, 'flag_start'] = True
    
    # Пересчёт параметров 
    o.J, o.r_center = call_inertia(o, o.taken_beams, app_y=id_app)
    o.J_1 = np.linalg.inv(o.J)
    o.R += o.S.T @ (o.r_center - r_center_p)
    o.V = (V_p * M_without + u_p * m_extra) / (M_without + m_extra)
    o.w = o.b_o(o.J_1, o.S) @ (o.b_o(J_p, o.S) @ o.w - (M_without + m_extra) * my_cross(o.R, o.V) +
                               m_extra * my_cross(r, u_p) + M_without * my_cross(R_p, V_p))
    o.Om_update()
    o.C_R = get_C_hkw(o.R, o.V, o.w_hkw)

    o.X_app.loc[id_app, 'busy_time'] = o.time_to_be_busy
    o.X_app.loc[id_app, 'flag_hkw'] = True
    o.X_app.loc[id_app, 'flag_fly'] = False
    o.t_reaction_counter = o.t_reaction
    o.t_start[o.N_app] = t
    o.flag_impulse = True
    o.taken_beams_p = o.taken_beams.copy()
