from mylibs.im_sample import *
from mylibs.numerical_methods import *
from mylibs.control_function import control_condition
import scipy


def call_crash_internal_func(r, r1, r2, diam, return_force=False, k_av=None, level: int = 5):
    """ Дополнительная функция для функции call_crash; \n
    Проверяет наличие точки r в цилиндре с концами r1,r2, диаметром diam; \n
    Возвращает {True,False} соответственно при отсутствии и наличии. """
    r1 = np.array(r1.copy())
    r2 = np.array(r2.copy())
    x1, y1, z1 = r1
    x2, y2, z2 = r2
    n = np.array([x2 - x1, y2 - y1, z2 - z1])  # Вектор вдоль стержня
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
        return forсe_from_beam(a, diam, n, tau, b, f0, f1, f2, k_av, level)
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
            if np.sum(o.s.flag[i]) > 0:
                r1 = o.s.r1[i]
                r2 = o.s.r2[i]
            else:
                r1 = [o.s.r_st[i][0], o.s.r_st[i][1], o.s.r_st[i][2]]
                r2 = [o.s.r_st[i][0] - o.s.length[i], o.s.r_st[i][1], o.s.r_st[i][2]]
            if call_crash_internal_func(r, r1, r2, o.d_crash):
                return True
    for i in range(o.N_cont_beams):
        r1 = o.c.r1[i]
        r2 = o.c.r2[i]
        if call_crash_internal_func(r, r1, r2, o.c.diam[i] * 1.05):
            return True
    return False


def call_inertia(o, id_s=np.array([]), app_y=None, app_n=None):
    """Функция считает тензор инерции в собственной ск конструкции и центр масс;    \n
    Input:                                                                          \n
    -> o - AllObjects класс;                                                        \n
    -> id_s - вектор id стержней, которых учитывать не надо                         \n
    -> app_y - номера аппаратов, которых стоит учесть при расчёте инерции;          \n
    -> app_n - номера аппаратов, которых не стоит учитывать при расчёте инерции;    \n
    Output:                                                                         \n
    -> Тензор инерции, вектор центра масс. """
    J = np.zeros((3, 3))
    N_points_m = 10                                                 # количество материальных точек на стержень
    N_m = N_points_m * (o.N_beams + o.N_cont_beams)                 # материальных точек всего
    m = [0. for _ in range(N_m + o.N_app)]
    r = [np.zeros(3) for _ in range(N_m + o.N_app)]

    # Выделение координат стержня - на месте или в грузовом отсеке
    for n in range(o.N_beams):
        if not np.any(id_s == n):
            if np.sum(o.s.flag[n]) > 0:
                r_1 = np.array(o.s.r1[n])
                r_2 = np.array(o.s.r2[n])
            else:
                r_1 = o.s.r_st[n]
                r_2 = o.s.r_st[n] - np.array([o.s.length[n], 0, 0])

            for i in range(N_points_m):
                m[n * N_points_m + i] = o.s.mass[n] / N_points_m
                r[n * N_points_m + i] = r_1 + (r_2 - r_1) * i / N_points_m

    # Учёт грузового отсека на тензор инерции
    for n in range(o.N_cont_beams):
        r_1 = np.array(o.c.r1[n])
        r_2 = np.array(o.c.r2[n])
        for i in range(N_points_m):
            m[(o.N_beams + n) * N_points_m + i] = o.c.mass[n] / N_points_m
            if n > 0:
                r[(o.N_beams + n) * N_points_m + i] = r_1 + (r_2 - r_1) * i / N_points_m
            else:
                r_0 = (r_1 + r_2) / 2
                b = o.c.diam[n] * my_cross((r_2 - r_1), [0, 1, 0]) / np.linalg.norm(r_2 - r_1)
                k = o.c.diam[n] * my_cross((r_2 - r_1), b) / np.linalg.norm(r_2 - r_1) / np.linalg.norm(b)
                r[(o.N_beams + n) * N_points_m + i] = r_0 + b * np.sin(2*np.pi * i / N_points_m) + \
                                                            k * np.cos(2*np.pi * i / N_points_m)

    # Учёт аппаратов
    for a in range(o.N_app):
        if (not o.a.flag_fly[a] and app_n != a) or app_y == a:
            m[N_m + a] = o.a.mass[a]
        else:
            m[N_m + a] = 0
            r[N_m + a] = o.a.target[a]

    # Вычисление центра масс
    tmp_numerator = np.zeros(3)
    tmp_denominator = np.zeros(3)
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
    X_middles = []
    for i in range(len(X_cont.id)):
        if X_cont.flag_grab[i]:
            X_middles.append(i)
    difference = [X_cont.r1[X_middles[i]] / 2 + X_cont.r2[X_middles[i]] / 2 for i in range(len(X_middles))]
    N_middles = len(difference)
    difs = np.zeros(N_middles)
    for i in range(N_middles):
        difs[i] = np.linalg.norm(difference[i] - r_right)
    i_needed = np.argmin(difs)

    return np.array(X_middles)[i_needed]


def calculation_motion(o, u, T_max, id_app, interaction=True, line_return=False, check_visible=False):
    """Функция расчитывает угловое и поступательное движение одного аппарата и конструкции; \n
    Других аппаратов и их решения она не учитывает;                                         \n
    Input:                                                                                  \n
    -> o - AllObjects класс;                                                                \n
    -> u - начальная скорость аппарата (ССК при interaction=True, ОСК иначе)                \n
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
    o_lcl = o.copy()

    e_max, V_max, R_max, j_max, f_mean = np.zeros(5)
    o_lcl.flag_vision[id_app] = False
    n_crashes = 0
    line = []

    if interaction:
        o_lcl.repulsion_change_params(id_app=id_app, u0=u)

    f_min = o_lcl.S @ o_lcl.get_discrepancy(id_app=id_app, vector=True) if interaction else \
        o_lcl.get_discrepancy(id_app=id_app, vector=True)

    # Iterations
    for i in range(int(T_max // o_lcl.dt)):
        # Writing parameters
        f = o_lcl.S @ o_lcl.get_discrepancy(id_app=id_app, vector=True) if interaction else \
            o_lcl.get_discrepancy(id_app=id_app, vector=True)
        if np.linalg.norm(f_min) > np.linalg.norm(f):
            f_min = f.copy()
        e_max = max(e_max, o_lcl.get_e_deviation())
        V_max = max(V_max, np.linalg.norm(o_lcl.V))
        R_max = max(R_max, np.linalg.norm(o_lcl.R))
        j_max = max(j_max, 180/np.pi*np.arccos((np.trace(o_lcl.S)-1)/2))
        f_mean += np.linalg.norm(f)
        line = np.append(line, o_lcl.o_b(o_lcl.a.r[id_app]))
        if check_visible:
            tmp_a = control_condition(o_lcl, id_app)
            # print(f"----{tmp_a} {o_lcl.flag_vision[id_app]}")

        # Stop Criteria
        if o_lcl.d_to_grab is not None:
            if np.linalg.norm(f) <= o_lcl.d_to_grab * 0.95:
                break
        if o_lcl.d_crash is not None:
            if (o_lcl.t - o.t) > 10:
                if call_crash(o, o_lcl.a.r[id_app], o_lcl.R, o_lcl.S, o_lcl.taken_beams):
                    n_crashes += 1
                    break
        if not o_lcl.control and \
                np.linalg.norm(o_lcl.a.r[id_app] - o_lcl.b_o(o_lcl.a.target_p[id_app])) > 7 * o_lcl.a.r_0[id_app]:
            break
        if i % 5 == 0:
            o_lcl.file_save(f'график {id_app} {np.linalg.norm(f)} {o_lcl.get_e_deviation()} '
                            f'{np.linalg.norm(180 / np.pi * np.arccos(clip((np.trace(o_lcl.S) - 1) / 2, -1, 1)))} '
                            f'{np.linalg.norm(o_lcl.V)} {np.linalg.norm(o_lcl.R)} '
                            f'{np.linalg.norm(o_lcl.a_self[id_app])}')

        # Iterating
        o_lcl.time_step()

    if (o_lcl.iter - o.iter) > 0:
        f_mean /= (o_lcl.iter - o.iter)

    if interaction is False:   # Пристыковка
        o_lcl.capturing_change_params(id_app)
        V_max = np.linalg.norm(o_lcl.V)
        R_max = np.linalg.norm(o_lcl.R)
        e_max = o_lcl.get_e_deviation()

        for i in range(int((o_lcl.iter - o.iter) // 20)):
            o_lcl.time_step()

            e_max = max(e_max, o_lcl.get_e_deviation())
            V_max = max(V_max, np.linalg.norm(o_lcl.V))
            R_max = max(R_max, np.linalg.norm(o_lcl.R))
            j_max = max(j_max, 180/np.pi*np.arccos((np.trace(o_lcl.S)-1)/2))

    if line_return:
        return f_min, f_mean, e_max, V_max, R_max, j_max, n_crashes, line
    return f_min, f_mean, e_max, V_max, R_max, j_max, n_crashes, o_lcl.flag_vision[id_app]


def find_repulsion_velocity(o, id_app: int, target=None, interaction: bool = True, method: str = 'trust-constr'):
    if interaction and not o.control or not interaction:
        func = f_to_capturing
        tol = o.d_to_grab
        o.my_print(f"Попадание: tol={tol}", test=True)
    else:
        func = f_to_detour
        tol = o.a.r_0[0]
        o.my_print(f"Огибание: tol={tol}, cont={o.s.container_length}", test=True)

    if interaction:
        target = o.a.target[id_app] if target is None else target
        tmp = target - np.array(o.o_b(o.a.r[id_app]))
    else:
        target = o.b_o(o.a.target[id_app]) if target is None else target
        tmp = o.b_o(target) - np.array(o.a.r[id_app])
    u = o.u_min * tmp / np.linalg.norm(tmp)
    mtd = method  # 'SLSQP' 'TNC' 'trust-constr'
    opt = {'verbose': 3}  # , 'xtol': o.u_max, 'gtol': tol}  # ,
    T_max = 2 * (np.linalg.norm(tmp) / o.u_min)

    nonlinear_constraint = scipy.optimize.NonlinearConstraint(fun=lambda x: np.linalg.norm(x), lb=o.u_min, ub=o.u_max)
    res = scipy.optimize.minimize(func, u, args=(o, T_max, id_app, interaction, False),
                                  tol=tol, method=mtd, options=opt,
                                  bounds=((0, o.u_max), (0, o.u_max), (0, o.u_max)),
                                  constraints={'type': 'ineq',
                                               'fun': lambda x: (np.linalg.norm(x) - o.u_min) / (o.u_max - o.u_min)})
    # constraints=nonlinear_constraint
    u = o.cases['repulse_vel_control'](res.x)
    # u = res.x
    return u


def repulsion(o, id_app, u_a_priori=None):
    """Input:                                                                   \n
    -> o - AllObjects класс                                                     \n
    -> id_app - номер аппарата                                                  \n
    -> u_a_priori - заданный вектор скорости (для отладки управления)"""
    # Параметры до отталкивания
    o.taken_beams_p = o.taken_beams.copy()
    o.repulse_app_config(id_app=id_app)
    r_1 = o.a.target[id_app]

    target_is_reached = False
    if u_a_priori is not None:
        u0 = u_a_priori
    else:
        method_comps = o.method.split('+')
        if 'diffevolve' in method_comps:
            u0, _ = diff_evolve(f_to_detour, [o.u_min, o.u_max], True, o, o.T_max, id_app, True, False, True,
                                n_vec=o.diff_evolve_vectors, chance=0.5, f=0.8, len_vec=3, n_times=o.diff_evolve_times,
                                multiprocessing=True, print_process=True)
        else:
            tmp = r_1 - np.array(o.o_b(o.a.r[id_app]))
            u0 = o.u_min * tmp / np.linalg.norm(tmp)
        if 'shooting' in method_comps:
            u1, target_is_reached = calc_shooting(o=o, id_app=id_app, r_right=r_1, interaction=True, u0=u0)
            if target_is_reached:
                u0 = u1
        elif 'const-propulsion' in method_comps:
            v = calc_shooting(o=o, id_app=id_app, r_right=r_1, interaction=True, func=f_controlled_const, n=6,
                              u0=np.append(u0, np.zeros(3)))
            u0 = v[0:3]
            o.a_self[id_app] = o.cases['acceleration_control'](v[3:6])
            print(f"TEPER {o.a_self[0]}")
        else:
            u0 = find_repulsion_velocity(o=o, id_app=id_app, target=r_1, interaction=True)

    if not target_is_reached and 'pd' in method_comps:
        o.my_print('Я не справляюсь, включаю ПД', mode="test")
        o.control = True
        o.if_PID_control = True
    if not target_is_reached and 'imp' in method_comps:
        o.my_print('Я не справляюсь, включаю импульсы', mode="test")
        o.control = True
        o.if_impulse_control = True

    o.repulsion_change_params(id_app=id_app, u0=u0)

    o.my_print(f"App id:{id_app} pushed off with u={np.linalg.norm(u0)}, w={np.linalg.norm(o.w)}", mode="b", test=True)
    o.my_print(f"Taken beams list: {o.taken_beams}", mode="g", test=True)
    talk_decision(o.if_talk)

    return u0


def capturing(o, id_app):
    if o.method == 'shooting+pd':
        o.control = False
        o.if_PID_control = False
    if o.method == 'shooting+imp':
        o.control = False
        o.if_impulse_control = False

    if o.if_any_print:
        print(Fore.BLUE + f'Аппарат id:{id_app} захватился')
    o.a_self[id_app] = np.array(np.zeros(3))
    # talk_success(o.if_talk)

    id_beam = o.a.flag_beam[id_app]
    if id_beam is not None:
        if np.linalg.norm(np.array(o.s.r1[id_beam]) - np.array(o.a.target[0])) < 1e-2:
            o.my_print(f"Стержень id:{id_beam} устанавливается", mode="b")
            o.s.flag[id_beam] = np.array([1., 1.])
            o.a.flag_beam[id_app] = None
            o.taken_beams = np.delete(o.taken_beams, np.argmax(o.taken_beams == id_beam))
    else:
        if o.a.target[id_app][0] < -0.6:  # Если "слева" нет промежуточных точек, то окей
            o.my_print(f'Аппарат id:{id_app} в грузовом отсеке')
            o.a.flag_start[id_app] = True

    o.capturing_change_params(id_app)
