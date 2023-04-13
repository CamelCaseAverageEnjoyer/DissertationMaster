from colorama import Fore, Style
import matplotlib.pyplot as plt
from p_tqdm import p_map
import numpy as np
import random
import copy
from mylibs.tiny_functions import *


def capturing_penalty(o, dr, dr_average, e, V, R, j, n_crashes, visible):
    reserve_rate = 1.5
    e_e = 0. if (o.e_max / reserve_rate - e) > 0 else abs(e - o.e_max / reserve_rate)
    e_V = 0. if (o.V_max / reserve_rate - V) > 0 else abs(V - o.V_max / reserve_rate)
    e_R = 0. if (o.R_max / reserve_rate - R) > 0 else abs(R - o.R_max / reserve_rate)
    e_j = 0. if (o.j_max / reserve_rate - j) > 0 else abs(j - o.j_max / reserve_rate)
    return dr + dr_average * 0.3 - o.mu_ipm * (np.log(o.e_max - e + e_e) + np.log(o.V_max - V + e_V) +
                            np.log(o.R_max - R + e_R) + np.log(o.j_max - j + e_j)) + \
           o.mu_e * (e_e / o.w_max + e_V / o.V_max + e_R / o.R_max + e_j / o.j_max)

def detour_penalty(o, dr, dr_average, e, V, R, j, n_crashes, visible):
    if False:  # (o.e_max - e > 0) and (o.V_max - V > 0) and (o.R_max - R > 0) and (o.j_max - j > 0):  # Constraint's fulfillment
        anw = dr + dr_average - o.mu_ipm * (np.log(o.e_max - e) + np.log(o.V_max - V) +
                               np.log(o.R_max - R) + np.log(o.j_max - j)) + \
              np.array([1e3, 1e3, 1e3]) * n_crashes
    else:
        anw = dr + dr_average + 1e2 * n_crashes - abs(1e2 * visible) # np.array([1e2, 1e2, 1e2]) * (clip(10*(e - o.e_max), 0, 1) + clip(10*(V - o.V_max), 0, 1) + clip(1e-1 * (R - o.R_max), 0, 1) + clip(1e-2 * (j - o.j_max), 0, 1)) + \
        print(f"{n_crashes} {1e2 * visible} |  {np.linalg.norm(dr + dr_average * 0.1)}  | {np.linalg.norm(anw)}")
    return anw

def f_to_capturing(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, return_to_shooting_method, check_visible = args
    u = o.cases['repulse_vel_control'](u)
    dr, dr_average, e, V, R, j, n_crashes, visible = calculation_motion(o=o, u=u, T_max=T_max, id_app=id_app,
                                                                        interaction=interaction,
                                                                        check_visible=check_visible)
    anw = capturing_penalty(o, dr, dr_average, e, V, R, j, n_crashes, visible)
    if return_to_shooting_method:
        return anw, np.linalg.norm(dr)
    return np.linalg.norm(anw)

def f_to_detour(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, return_to_shooting_method, check_visible = args
    u = o.cases['repulse_vel_control'](u)
    dr, dr_average, e, V, R, j, n_crashes, visible = calculation_motion(o, u, T_max, id_app, interaction=interaction,
                                                                        check_visible=check_visible)
    anw = detour_penalty(o, dr, dr_average, e, V, R, j, n_crashes, visible)
    if return_to_shooting_method:
        return anw, np.linalg.norm(dr)
    return np.linalg.norm(anw)

def f_controlled_const(v, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, return_to_shooting_method, check_visible = args
    u = o.cases['repulse_vel_control'](v[0:3])
    o.a_self = o.cases['acceleration_control'](v[3:6])
    dr, dr_average, e, V, R, j, n_crashes, visible = calculation_motion(o, u, T_max, id_app, interaction=True,
                                                                        check_visible=check_visible)
    return capturing_penalty(o, dr, dr_average, e, V, R, j, n_crashes, visible), np.linalg.norm(dr)

def calc_shooting_sample(u, o, T_max, id_app, interaction, func):
    dr, tol = func(u, o, T_max, id_app, interaction, True, False)
    return [i for i in dr] + [tol]

def calc_shooting(o, id_app, r_right, interaction: bool = True, u0: any = None, n: int = 3, func: any = f_to_capturing):
    """ Функция выполняет пристрелочный/спектральный поиск оптимальной скорости отталкивания/импульса аппарата; \n
    Input:                                                                                  \n
    -> o - AllObjects класс;                                                                \n
    -> id_app - номер аппарата                                                              \n
    -> r_right - радиус-вектор ССК положения цели                                           \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции      test  \n
    Output:                                                                                 \n
    -> u - оптимальный вектор скорости отталкивания/импульса (ССК при interaction=True, ОСК иначе)"""
    shooting_amount = o.shooting_amount_repulsion if interaction else o.shooting_amount_impulse
    if interaction:
        tmp = r_right - np.array(o.o_b(o.a.r[id_app]))
    else:
        tmp = o.b_o(r_right) - np.array(o.a.r[id_app])
    mu_e = o.mu_e
    T_max = o.T_max  # 2*(np.linalg.norm(tmp)/o.u_min)  # if o.a.flag_fly[id_app] else o.T_max
    u = np.array(u0)
    if n > 3:
        u[0:3] = o.u_min * tmp / np.linalg.norm(tmp) if u0 is None or np.linalg.norm(u0) < 1e-5 else u0[0:3]
    else:
        u = o.u_min * tmp / np.linalg.norm(tmp) if u0 is None else u0
    print(f"u: {u}")

    # Метод пристрелки
    mu_ipm = o.mu_ipm
    i_iteration = 0
    while i_iteration < shooting_amount:
        du = 0.001
        u_i = np.eye(n) * du
        mu_ipm /= 5
        mu_e /= 2
        i_iteration += 1

        anw = p_map(calc_shooting_sample,
                    [u] + [u + u_i[:, i] for i in range(n)],
                    [o for _ in range(1 + n)],
                    [T_max for _ in range(1 + n)],
                    [id_app for _ in range(1 + n)],
                    [interaction for _ in range(1 + n)],
                    [func for _ in range(1 + n)])
        dr = np.array(anw[0][0:3])
        dr_ = [anw[1 + i] for i in range(n)]
        '''dr_x = np.array(anw[1][0:3])
        dr_y = np.array(anw[2][0:3])
        dr_z = np.array(anw[3][0:3])'''
        tol = anw[0][3]
        # print(f"dr {dr}")
        # print(f"dr_ {dr_}")
        # print(f"[u] + [u + u_i[:, i] for i in range(n)] {[u] + [u + u_i[:, i] for i in range(n)]}")
        o.my_print(f"Точность пристрелки {tol}, разрешённое время {T_max} секунд, целевая функция {np.linalg.norm(dr)}",
                   mode=None)
        if o.d_to_grab is not None and tol < o.d_to_grab*0.98:  # and not (c or cx or cy or cz):
            break

        Jacobian = np.array([[dr_[j][i] - dr[i] for j in range(n)] for i in range(3)]) / du
        # print(f"J {Jacobian}")
        if True:  # np.linalg.det(Jacobian) > 1e-7:
            # Jacobian_1 = np.linalg.inv(Jacobian)
            Jacobian_1 = np.linalg.pinv(Jacobian)
            correction = Jacobian_1 @ dr
            correction = correction/np.linalg.norm(correction)*clip(np.linalg.norm(correction), 0, o.u_max/2)
        else:
            print(f"ВCЁ ПЛОХО ТАК ТО ГОВОРЯ")
            correction = u * np.array([random.uniform(-1, 1)] * n) * np.linalg.norm(u) * o.k_u
        u = u - correction

        # Ограничение по модулю скорости
        if np.linalg.norm(u[0:3]) > o.u_max:
            o.my_print(f'attention: shooting speed reduction: {np.linalg.norm(u[0:3])/o.u_max*100} %')
            u[0:3] = u[0:3] / np.linalg.norm(u[0:3]) * o.u_max * 0.9
        if n > 3 and np.linalg.norm(u[3:6]) > o.a_pid_max:
            o.my_print(f'attention: acceleration reduction: {np.linalg.norm(u[3:6])/o.a_pid_max*100} %')
            u[3:6] = u[3:6] / np.linalg.norm(u[3:6]) * o.a_pid_max * 0.9

        if np.linalg.norm(dr) > tol * 1e2 and i_iteration == shooting_amount:
            mu_ipm *= 10
            i_iteration -= 1
            T_max *= 1.1 if T_max < o.T_max_hard_limit else 1
    if o.method in ['shooting+pd', 'shooting+imp']:
        return u, tol < o.d_to_grab*0.98
    return u

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
        v = np.array([get_v(np.exp(random.uniform(np.log(search_domain[0]), np.log(search_domain[1]))),
                            random.uniform(0, 2 * np.pi),
                            random.uniform(- np.pi / 2, np.pi / 2)) for _ in range(n_vec)])
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
