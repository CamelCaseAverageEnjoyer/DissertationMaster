from colorama import Fore, Style
import matplotlib.pyplot as plt
from p_tqdm import p_map
import numpy as np
import random
import copy
from mylibs.tiny_functions import *

def f_to_capturing(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, return_to_shooting_method = args
    u = o.cases['repulse_vel_control'](u)
    dr, _, w, V, R, j, _ = calculation_motion(o=o, u=u, T_max=T_max, id_app=id_app, interaction=interaction)
    reserve_rate = 1.5
    w = abs(w - o.w_twist)
    e_w = 0. if (o.w_max/reserve_rate - w) > 0 else abs(w - o.w_max/reserve_rate)
    e_V = 0. if (o.V_max/reserve_rate - V) > 0 else abs(V - o.V_max/reserve_rate)
    e_R = 0. if (o.R_max/reserve_rate - R) > 0 else abs(R - o.R_max/reserve_rate)
    e_j = 0. if (o.j_max/reserve_rate - j) > 0 else abs(j - o.j_max/reserve_rate)
    anw = dr - \
        o.mu_ipm * (np.log(o.w_max - w + e_w) + np.log(o.V_max - V + e_V) +
                    np.log(o.R_max - R + e_R) + np.log(o.j_max - j + e_j)) + \
        o.mu_e * (e_w/o.w_max + e_V/o.V_max + e_R/o.R_max + e_j/o.j_max)
    if return_to_shooting_method:
        return anw, np.linalg.norm(dr)
    return np.linalg.norm(anw)

def f_to_detour(u, *args):
    from mylibs.calculation_functions import calculation_motion
    if len(args) == 1:
        args = args[0]
    o, T_max, id_app, interaction, return_vec = args
    print(f"скорость: {np.linalg.norm(u)}")
    u = o.cases['repulse_vel_control'](u)
    dr, dr_average, w, V, R, j, n_crashes = calculation_motion(o, u, T_max, id_app, interaction=interaction)
    w = abs(w - o.w_twist)

    if (o.w_max - w > 0) and (o.V_max - V > 0) and (o.R_max - R > 0) and (o.j_max - j > 0):  # Constraint's fulfillment
        anw = dr - o.mu_ipm * (np.log(o.w_max - w) + np.log(o.V_max - V) +
                               np.log(o.R_max - R) + np.log(o.j_max - j)) + \
              np.array([1e3, 1e3, 1e3]) * n_crashes
    else:
        anw = np.array([1e2, 1e2, 1e2]) * \
              (clip(10*(w - o.w_max), 0, 1) + clip(10*(V - o.V_max), 0, 1) +
               clip(1e-1 * (R - o.R_max), 0, 1) + clip(1e-2 * (j - o.j_max), 0, 1)) + \
              np.array([1e3, 1e3, 1e3]) * n_crashes
    if return_vec:
        return anw
    return np.linalg.norm(anw)

def calc_shooting(o, id_app, r_right, interaction=True, u0=None):
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
    # T_max = o.T_max
    # T_max_hard_limit = T_max*16 if o.X_app.flag_fly[id_app] else o.T_max_hard_limit
    u = o.u_min * tmp / np.linalg.norm(tmp) if u0 is None else np.array(u0)

    # Метод пристрелки
    mu_ipm = o.mu_ipm
    i_iteration = 0
    while i_iteration < shooting_amount:
        du = 0.00001
        u_i = np.eye(3) * du
        mu_ipm /= 5
        mu_e /= 2
        i_iteration += 1

        dr, tol = f_to_capturing(u, o, T_max, id_app, interaction, True)
        dr_x, _ = f_to_capturing(u + u_i[:, 0], o, T_max, id_app, interaction, True)
        dr_y, _ = f_to_capturing(u + u_i[:, 1], o, T_max, id_app, interaction, True)
        dr_z, _ = f_to_capturing(u + u_i[:, 2], o, T_max, id_app, interaction, True)
        o.my_print(f"Точность пристрелки {tol}, разрешённое время {T_max} секунд", mode=None)
        if o.d_to_grab is not None and tol < o.d_to_grab*0.98:  # and not (c or cx or cy or cz):
            break

        Jacobian = np.array([[dr_x[0] - dr[0], dr_y[0] - dr[0], dr_z[0] - dr[0]],
                             [dr_x[1] - dr[1], dr_y[1] - dr[1], dr_z[1] - dr[1]],
                             [dr_x[2] - dr[2], dr_y[2] - dr[2], dr_z[2] - dr[2]]]) / du
        if np.linalg.det(Jacobian):
            Jacobian_1 = np.linalg.inv(Jacobian)
            correction = Jacobian_1 @ dr
            correction = correction/np.linalg.norm(correction)*clip(np.linalg.norm(correction), 0, o.u_max/2)
        else:
            correction = u * np.array([random.uniform(-1, 1)] * 3) * np.linalg.norm(u) * o.k_u
        u = u - correction

        # Ограничение по модулю скорости
        if np.linalg.norm(u) > o.u_max:
            if o.if_any_print:
                print(Style.RESET_ALL + f'attention: shooting speed reduction: {np.linalg.norm(u)/o.u_max*100} %')
            u = u / np.linalg.norm(u) * o.u_max * 0.9

        if np.linalg.norm(dr) > tol * 1e2:
            mu_ipm *= 10
            i_iteration -= 1
            T_max *= 1.1 if T_max < o.T_max_hard_limit else 1
    if o.method == 'shooting+pd':
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


def diff_evolve(func: any, search_domain: list, *args, **kwargs):
    """Функция дифференциальной эволюции.
    :param func: целевая функция
    :param search_domain: 2-мерный список разброса вектора аргументов: [[v[0]_min, v[0]_max], [v[1]_min, v[1]_max],...]
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
