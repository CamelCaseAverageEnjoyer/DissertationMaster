import random
from colorama import Fore, Style
from p_tqdm import p_map
import copy
from mylibs.tiny_functions import *


def capturing_penalty(o, dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper, mu_ipm):
    reserve_rate = 0.0001
    '''e_V = reserve_rate if (o.V_max - V) > 0 else abs((o.V_max - V)/o.V_max + reserve_rate*np.exp((o.V_max - V)/o.V_max))
    e_R = reserve_rate if (o.R_max - R) > 0 else abs((o.R_max - R)/o.R_max + reserve_rate*np.exp((o.R_max - R)/o.R_max))
    e_w = reserve_rate if (o.w_max - w) > 0 else abs((o.w_max - w)/o.w_max + reserve_rate*np.exp((o.w_max - w)/o.w_max))
    e_j = reserve_rate if (o.j_max - j) > 0 else abs((o.j_max - j)/o.j_max + reserve_rate*np.exp((o.j_max - j)/o.j_max))'''
    e_w = 1 - abs(w)/o.w_max + reserve_rate if (1 - abs(w)/o.w_max) > 0 else reserve_rate * np.exp(1 - abs(w)/o.w_max)
    e_j = 1 - abs(j)/o.j_max + reserve_rate if (1 - abs(j)/o.j_max) > 0 else reserve_rate * np.exp(1 - abs(j)/o.j_max)
    e_V = 1 - abs(V)/o.V_max + reserve_rate if (1 - abs(V)/o.V_max) > 0 else reserve_rate * np.exp(1 - abs(V)/o.V_max)
    e_R = 1 - abs(R)/o.R_max + reserve_rate if (1 - abs(R)/o.R_max) > 0 else reserve_rate * np.exp(1 - abs(R)/o.R_max)
    # print(f"eV={e_V}, eR={e_R}, ew={e_w}, ej={e_j} {(o.w_max - w)/o.w_max}")
    id_app = 0
    a = o.a.target[id_app] - o.o_b(o.r_ub)
    tau = my_cross(dr, a + dr)
    tau /= np.linalg.norm(tau)
    b = my_cross(tau, dr)
    b /= np.linalg.norm(b)
    tmp = abs(np.dot(dr / np.linalg.norm(dr), a / np.linalg.norm(a))) * 1e2 
    # anw = np.linalg.norm(dr) * (clip(1 - crhper, 0, 1) + clip(crhper, 0, 1) * tmp)
    anw = dr * (clip(1 - crhper, 0, 1) + clip(crhper, 0, 1) * tmp)
    anw += (-mu_ipm * (np.log(e_w) + np.log(e_V) + np.log(e_R) + np.log(e_j))) * dr / np.linalg.norm(dr)
    '''min_local = -1e10
    anw += (mu_ipm * (clip(w/o.w_max - 1, min_local, 1e10)**2 + clip(j/o.j_max - 1, min_local, 1e10)**2 +
                      clip(R/o.R_max - 1, min_local, 1e10)**2 + clip(V/o.V_max - 1, min_local, 1e10)**2)) * \
           dr / np.linalg.norm(dr)'''
    return anw

def detour_penalty(o, dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper, mu_ipm):
    '''# Constraint's fulfillment
    if False:  # (o.e_max - e > 0) and (o.V_max - V > 0) and (o.R_max - R > 0) and (o.j_max - j > 0):
        anw = dr + dr_average - o.mu_ipm * (np.log(o.e_max - e) + np.log(o.V_max - V) +
                               np.log(o.R_max - R) + np.log(o.j_max - j)) + \
              np.array([1e3, 1e3, 1e3]) * n_crashes'''
    # np.array([1e2, 1e2, 1e2]) * (clip(10*(e - o.e_max), 0, 1) + clip(10*(V - o.V_max), 0, 1) + clip(1e-1 * (R - o.R_max), 0, 1) + clip(1e-2 * (j - o.j_max), 0, 1)) + \

    anw = np.linalg.norm(dr) + dr_average * 0.1 + 1 * n_crashes + 1e2 * (not visible)
    params = [[o.e_max, e], [o.j_max, j], [o.V_max, V], [o.R_max, R]]
    for i in range(2):
        if params[i][0] - params[i][1] > 0:  # Constraint's fulfillment
            anw -= o.mu_ipm * np.log(params[i][0] - params[i][1])
        else:
            anw += 1e2
    # print(f"{n_crashes} {1e2 * visible} |  {np.linalg.norm(dr)}:{dr_average}  | {np.linalg.norm(anw)}")
    return anw

def f_scipy(u, *args):
    from mylibs.calculation_functions import calculation_motion
    o, T_max, id_app, interaction, check_visible, mu_ipm = args
    f_min, e_max, V_max, R_max, w_max, j_max, g = calculation_motion(o=o, u=u, T_max=T_max, id_app=id_app,
                                                                     interaction=interaction, to_scipy=True)

def f_dr(u, *args):
    from mylibs.calculation_functions import calculation_motion
    o, T_max, id_app, interaction, check_visible, mu_ipm = args
    dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper = calculation_motion(o=o, u=u, T_max=T_max, id_app=id_app,
                                                                                   interaction=interaction,
                                                                                   check_visible=False)
    return np.linalg.norm(dr)**2

def f_to_capturing(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, check_visible, mu_ipm = args
    u = o.cases['repulse_vel_control'](u)
    dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper = calculation_motion(o=o, u=u, T_max=T_max, id_app=id_app,
                                                                                   interaction=interaction,
                                                                                   check_visible=False)
    return capturing_penalty(o, dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper, mu_ipm), np.linalg.norm(dr)

def f_to_detour(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, check_visible, mu_ipm = args
    u = o.cases['repulse_vel_control'](u)
    dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper = calculation_motion(o, u, T_max, id_app,
                                                                                   interaction=interaction,
                                                                                   check_visible=check_visible)
    anw = detour_penalty(o, dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper, mu_ipm)
    return anw

def f_controlled_const(v, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, check_visible, mu_ipm = args
    u = o.cases['repulse_vel_control'](v[0:3])
    if o.if_T_in_shooting:
        tmp = len(v) - 1
    else:
        tmp = len(v) 
    dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper = calculation_motion(o, u, T_max, id_app,
                                                                                   interaction=interaction,
                                                                                   check_visible=check_visible,
                                                                                   control=v[3:tmp])
    return capturing_penalty(o, dr, dr_average, e, V, R, w, j, n_crashes, visible, crhper, mu_ipm)

def my_calc_shooting(o, id_app, r_1, interaction: bool = True, u: any = None, func: any = f_dr):
    """u = [T, u_x, u_y, u_z]"""
    shooting_amount = o.shooting_amount_repulsion if interaction else o.shooting_amount_impulse
    mu_e = o.mu_e
    i_iteration = 0

    def local_correction(t, w, v_x, v_y, v_z, x0, y0, z0, x1, y1, z1):
        from numpy import sin, cos
        return np.array([2*(-(3*t*w - 4*sin(t*w))*(2*v_z*cos(t*w) - 2*v_z - w*(3*t*(v_x + 2*w*z0) - x0 + x1) +
                                                   (4*v_x + 6*w*z0)*sin(t*w)) -
                            2*(cos(t*w) - 1)*(2*v_x + v_z*sin(t*w) + w*(4*z0 - z1) +
                                              (-2*v_x - 3*w*z0)*cos(t*w)))/w**2,
                         2*(-(3*t*w - 4*sin(t*w))*(2*v_z*cos(t*w) - 2*v_z - w*(3*t*(v_x + 2*w*z0) - x0 + x1) +
                                                   (4*v_x + 6*w*z0)*sin(t*w)) - 2*(cos(t*w) - 1) *
                            (2*v_x + v_z*sin(t*w) + w*(4*z0 - z1) + (-2*v_x - 3*w*z0)*cos(t*w)))/w**2,
                         2*(v_y*sin(t*w) + w*(y0*cos(t*w) - y1))*sin(t*w)/w**2,
                         2*(2*(cos(t*w) - 1)*(2*v_z*cos(t*w) - 2*v_z - w*(3*t*(v_x + 2*w*z0) - x0 + x1) +
                                              (4*v_x + 6*w*z0)*sin(t*w)) + (2*v_x + v_z*sin(t*w) + w*(4*z0 - z1) +
                                                                            (-2*v_x - 3*w*z0)*cos(t*w))*sin(t*w))/w**2])

    while i_iteration < shooting_amount:
        i_iteration += 1
        f = func(u[1:4], o, u[0], id_app, interaction, False, o.mu_ipm)
        print(f"Итерация {i_iteration}|dr={round(np.sqrt(f),5)}, T={u[0]}, u={u[1:4]}")
        if f > o.d_to_grab**2:
            v = o.b_o(u[1:4])
            r0 = o.b_o(o.a.target_p[id_app])
            r1 = o.b_o(r_1)
            u -= local_correction(u[0], o.w_hkw, v[0], v[1], v[2], r0[0], r0[1], r0[2], r1[0], r1[1], r1[2])
        pass

def calc_shooting_sample(u, o, T_max, id_app, interaction, mu_ipm, func):
    T_max = u[len(u) - 1] if o.if_T_in_shooting else T_max
    f, tol = func(u, o, T_max, id_app, interaction, o.method in ['linear-propulsion', 'const-propulsion'], mu_ipm)
    return [i for i in f] + [tol]

def calc_shooting(o, id_app, r_1, interaction: bool = True, u0: any = None, n: int = 3, func: any = f_to_capturing, T_max=None):
    """ Функция выполняет пристрелочный/спектральный поиск оптимальной скорости отталкивания/импульса аппарата; \n
    Input:                                                                                  \n
    -> o - AllObjects класс;                                                                \n
    -> id_app - номер аппарата                                                              \n
    -> r_right - радиус-вектор ССК положения цели                                           \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции        \n
    -> u0 - начальное приближение                                                           \n
    -> n - длина рабочего вектора                                                           \n
    -> func - функция, выдающая вектор длины 3, минимизируемая оптимальным вход-вектором    \n
    Output:                                                                                 \n
    -> u - оптимальный вектор скорости отталкивания/импульса (ССК при interaction=True, ОСК иначе)"""
    shooting_amount = o.shooting_amount_repulsion if interaction else o.shooting_amount_impulse
    tmp = r_1 - o.o_b(o.a.r[id_app]) if interaction else o.b_o(r_1) - np.array(o.a.r[id_app])
    mu_e = o.mu_e
    T_max = o.T_max if T_max is None else T_max
    u = o.u_min * tmp / np.linalg.norm(tmp) if u0 is None else u0
    if n > 3:
        u = np.zeros(n)
        u[0:3] = o.u_min * tmp / np.linalg.norm(tmp) if u0 is None or np.linalg.norm(u0) < 1e-5 else u0[0:3]
    u_anw = u.copy()
    tol_anw = 1e5
    tol_list = []

    # Метод пристрелки
    mu_ipm = o.mu_ipm
    i_iteration = 0
    n_full = n
    if o.if_T_in_shooting:
        u = np.append(u, T_max)
        n_full += 1
    while i_iteration < shooting_amount:
        du = 1e-7
        if o.if_T_in_shooting:
            u_i = np.diag([du * o.u_max] * 3 + [du * o.a_pid_max] * (n - 3) + [10])
        else:
            u_i = np.diag([du * o.u_max] * 3 + [du * o.a_pid_max] * (n - 3))
        mu_ipm /= 2
        mu_e /= 2
        i_iteration += 1

        anw = p_map(calc_shooting_sample,
                    [u] + [u + u_i[:, i] for i in range(n_full)],
                    [o for _ in range(1 + n_full)],
                    [T_max for _ in range(1 + n_full)],
                    [id_app for _ in range(1 + n_full)],
                    [interaction for _ in range(1 + n_full)],
                    [mu_ipm for _ in range(1 + n_full)],
                    [func for _ in range(1 + n_full)])
        dr = np.array(anw[0][0:3])
        dr_ = [anw[1 + i] for i in range(n_full)]
        tol = anw[0][3]
        tol_list += [tol]
        if o.if_T_in_shooting:
            T_max = u[len(u) - 1]
        o.my_print(f"Точность {round(tol,5)}, целевая функция {round(np.linalg.norm(dr),5)}, T_max={T_max}",
                   mode=None)
        if np.linalg.norm(dr) < tol_anw:
            tol_anw = np.linalg.norm(dr)
            u_anw = u.copy()
   
        if o.d_to_grab is not None and tol < o.d_to_grab*0.98:  # and not (c or cx or cy or cz):
            tol_anw = tol
            u_anw = u.copy()
            break
        else:
            Jacobian = np.array([[(dr_[j][i] - dr[i]) / u_i[j][j] for j in range(n_full)] for i in range(3)])
            if n_full == 3 and np.linalg.det(Jacobian) > 1e-7:
                Jacobian_1 = np.linalg.inv(Jacobian)
            else:
                Jacobian_1 = np.linalg.pinv(Jacobian)
            correction = Jacobian_1 @ dr
            if np.linalg.norm(correction) > 1e-12:
                correction = correction / np.linalg.norm(correction) * clip(np.linalg.norm(correction), 0, o.u_max / 4)
            else:
                print(f"Пристрелка закончена, градиента нет")
                break
            u -= correction

            # Ограничение по модулю скорости
            if np.linalg.norm(u[0:3]) > o.u_max:
                o.my_print(f'attention: shooting speed reduction: {np.linalg.norm(u[0:3])/o.u_max*100} %')
                u[0:3] = u[0:3] / np.linalg.norm(u[0:3]) * o.u_max * 0.9
            if o.method != 'linear-angle':
                if n > 3 and np.linalg.norm(u[3:6]) > o.a_pid_max:
                    o.my_print(f'attention: acceleration reduction: {np.linalg.norm(u[3:6])/o.a_pid_max*100} %')
                    u[3:6] = u[3:6] / np.linalg.norm(u[3:6]) * o.a_pid_max * 0.9
                if n > 6 and np.linalg.norm(u[6:9]) > o.a_pid_max:
                    o.my_print(f'attention: acceleration reduction: {np.linalg.norm(u[6:9])/o.a_pid_max*100} %')
                    u[6:9] = u[6:9] / np.linalg.norm(u[6:9]) * o.a_pid_max * 0.9

            '''if np.linalg.norm(dr) > tol * 1e2 and i_iteration == shooting_amount:
                mu_ipm *= 10
                i_iteration -= 1
                T_max *= 1.1 if T_max < o.T_max_hard_limit else 1'''
    method_comps = o.method.split('+')
    file_local = open('storage/iteration_docking.txt', 'a')
    if True:  # tol < o.d_to_grab*0.999:
        for i in range(len(tol_list)):
            file_local.write(f"{i} {tol_list[i]}\n")
    file_local.close()
    return u_anw, tol < o.d_to_grab*0.999  # and abs(tol - np.linalg.norm(dr)) < 1e-10

def diff_evolve_sample(j: int, func: any, v, target_p, comp_index: list,
                       chance: float = 0.5, f: float = 1., *args):
    args = args[0]
    mutant = v[comp_index[0]].copy() + f * (v[comp_index[1]] - v[comp_index[2]])
    for i in range(len(mutant)):
        if random.uniform(0, 1) < chance:
            mutant[i] = v[j][i].copy()
    target_p = func(v[j], args) if target_p is None else target_p
    target = func(mutant, args)
    v[j] = mutant.copy() if target < target_p else v[j]
    target_p = target if target < target_p else target_p
    return np.append(v[j], target_p)

def diff_evolve(func: any, search_domain: list, vector_3d: bool = False, *args, **kwargs):
    """Функция дифференциальной эволюции.
    :param func: целевая функция
    :param search_domain: 2-мерный список разброса вектора аргументов: [[v[0]_min, v[0]_max], [v[1]_min, v[1]_max],...]
    :param vector_3d: bla bla bla
    :return: v_best: len_vec-мерный список"""
    chance = 0.5 if 'chance' not in kwargs.keys() else kwargs['chance']
    f = 1. if 'f' not in kwargs.keys() else kwargs['f']
    n_vec = 10 if 'n_vec' not in kwargs.keys() else kwargs['n_vec']
    len_vec = 3 if 'len_vec' not in kwargs.keys() else kwargs['len_vec']
    n_times = 5 if 'n_times' not in kwargs.keys() else kwargs['n_times']
    multiprocessing = True if 'multiprocessing' not in kwargs.keys() else kwargs['multiprocessing']
    print_process = False if 'print_process' not in kwargs.keys() else kwargs['print_process']
    lst_errors = []
    # попробовать tuple(search_domain[i])
    if vector_3d:
        if len_vec == 3:
            v = np.array([polar2dec(np.exp(random.uniform(np.log(search_domain[0]), np.log(search_domain[1]))),
                                    random.uniform(0, 2 * np.pi),
                                    random.uniform(- np.pi / 2, np.pi / 2)) for _ in range(n_vec)])
        else:
            v = np.array([np.append(polar2dec(np.exp(random.uniform(np.log(search_domain[0]), np.log(search_domain[1]))),
                                              random.uniform(0, 2 * np.pi), random.uniform(- np.pi / 2, np.pi / 2)),
                                    [random.uniform(search_domain[2], search_domain[3])
                                     for _ in range(len_vec - 3)]) for _ in range(n_vec)])
    else:
        v = np.array([np.array([random.uniform(search_domain[i][0], search_domain[i][1]) for i in range(len_vec)])
                      for _ in range(n_vec)])
    v_record = [copy.deepcopy(v)]
    target_prev = [None for _ in range(n_vec)]
    v_best = None
    for i in range(n_times):
        if print_process:
            print(Fore.CYAN + f"Шаг {i + 1}/{n_times} дифференциальной эволюции" + Style.RESET_ALL)
        comp_index = [[] for _ in range(n_vec)]
        for j in range(n_vec):
            complement = list(range(n_vec))
            complement.remove(j)
            for _ in range(3):
                comp_index[j].append(random.choice(complement))
                complement.remove(comp_index[j][len(comp_index[j]) - 1])
        anw = p_map(diff_evolve_sample,
                    [j for j in range(n_vec)],
                    [func for _ in range(n_vec)],
                    [v for _ in range(n_vec)],
                    [target_prev[j] for j in range(n_vec)],
                    [comp_index[j] for j in range(n_vec)],
                    [chance for _ in range(n_vec)],
                    [f for _ in range(n_vec)],
                    [args for _ in range(n_vec)]) if multiprocessing else \
            [diff_evolve_sample(j, func, v,
                                target_prev[j],
                                comp_index[j],
                                chance, f, args) for j in range(n_vec)]
        v = np.array([np.array(anw[j][0:len_vec]) for j in range(n_vec)])
        v_record += [copy.deepcopy(v)]
        target_prev = [anw[j][len_vec] for j in range(n_vec)]
        lst_errors.append(np.min(target_prev))
        v_best = v[np.argmin(target_prev)]
    print(Fore.MAGENTA + f"Ошибка: {lst_errors}" + Style.RESET_ALL)
    print(Fore.MAGENTA + f"Ответ: {v_best}" + Style.RESET_ALL)
    # plt.plot(range(len(lst_errors)), lst_errors, c='indigo')
    # plt.show()
    return v_best, v_record
