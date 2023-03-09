# Standard libraries
from datetime import datetime
from p_tqdm import p_map
import random
import scipy

# Local libraries
from mylibs.plot_functions import *
from mylibs.calculation_functions import calculation_motion


def dr_gradient(u0, *args):
    o, T_max, id_app, interaction, mu_IPM, mu_e = args
    dr_, _, w_, V_, R_, j_, _ = calculation_motion(o=o, u=u0, T_max=T_max, id_app=id_app, interaction=interaction)
    reserve_rate = 1.5
    e_w = 0. if (o.w_max/reserve_rate - w_) > 0 else abs(w_ - o.w_max/reserve_rate)
    e_V = 0. if (o.V_max/reserve_rate - V_) > 0 else abs(V_ - o.V_max/reserve_rate)
    e_R = 0. if (o.R_max/reserve_rate - R_) > 0 else abs(R_ - o.R_max/reserve_rate)
    e_j = 0. if (o.j_max/reserve_rate - j_) > 0 else abs(j_ - o.j_max/reserve_rate)
    cnd = (e_w + e_V + e_R + e_j > 0)
    anw = dr_ - \
        mu_IPM * (np.log(o.w_max - w_ + e_w) + np.log(o.V_max - V_ + e_V) +
                  np.log(o.R_max - R_ + e_R) + np.log(o.j_max - j_ + e_j)) + \
        mu_e * (e_w/o.w_max + e_V/o.V_max + e_R/o.R_max + e_j/o.j_max), \
        cnd, np.linalg.norm(dr_)
    return np.linalg.norm(anw)


def dr_non_gradient(u, *args):
    o, T_max, id_app, interaction = args
    mu_IPM = 0.01
    dr_p, dr_m_p, w_p, V_p, R_p, j_p, N_crashes_p = calculation_motion(o, u, T_max, id_app, interaction=interaction)

    # Нет проблем с ограничениями
    if (o.w_max - w_p > 0) and (o.V_max - V_p > 0) and (o.R_max - R_p > 0) and (o.j_max - j_p > 0):
        anw = dr_p - mu_IPM * (np.log(o.w_max - w_p) + np.log(o.V_max - V_p) +
                               np.log(o.R_max - R_p) + np.log(o.j_max - j_p)) + \
            np.array([1e50, 1e50, 1e50]) * N_crashes_p
    # Нарушение ограничений
    else:
        anw = np.array([1e5, 1e5, 1e5]) * \
            (clip(10*(w_p - o.w_max), 0, 1) + clip(10*(V_p - o.V_max), 0, 1) +
             clip(1e-1*(R_p - o.R_max), 0, 1) + clip(1e-2*(j_p - o.j_max), 0, 1)) + \
            np.array([1e50, 1e50, 1e50]) * N_crashes_p
    return np.linalg.norm(anw)


def diff_evolve_sample(o, T_max, id_app, interaction, convex, tmp1x, tmp1y, tmp1z, tmp2x, tmp2y, tmp2z,
                       tmp3x, tmp3y, tmp3z, ajx, ajy, ajz, a_recjx, a_recjy, a_recjz, iter_decrease):
    """- Зачем я существую?         \n
    - Ты - условная единица парпрога"""
    tmp1 = np.array([tmp1x, tmp1y, tmp1z])
    tmp2 = np.array([tmp2x, tmp2y, tmp2z])
    tmp3 = np.array([tmp3x, tmp3y, tmp3z])
    aj = np.array([ajx, ajy, ajz])
    a_recj = np.array([a_recjx, a_recjy, a_recjz])
    u_p = o.X_app.v[id_app]
    mu_IPM = o.mu_ipm / 2 ** iter_decrease

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
            dr_p, _, _ = dr_gradient(o=o, u0=u_pre, T_max=T_max, id_app=id_app, interaction=interaction,
                                     mu_IPM=mu_IPM, mu_e=0.1)
        else:
            dr_p = dr_non_gradient(u_pre, o, T_max, id_app, interaction)
    else:
        dr_p = a_recj

    if np.linalg.norm(u_p) == 0.:
        du = u_pre - np.array(o.X_app.v[id_app])
        du_m = np.linalg.norm(du)
        dr_p = dr_p - mu_IPM * np.log([o.du_impulse_max - du_m, o.du_impulse_max - du_m, o.du_impulse_max - du_m]) \
            if o.du_impulse_max > du_m else np.array([1e10, 1e10, 1e10])

    if convex:
        dr_n, _, real_dr = dr_gradient(o=o, u0=u_new, T_max=T_max, id_app=id_app, interaction=interaction,
                                       mu_IPM=mu_IPM, mu_e=o.mu_e)
    else:
        dr_n = dr_non_gradient(u_new, o, T_max, id_app, interaction)

    if np.linalg.norm(u_p) == 0.:
        du = u_new - np.array(o.X_app.v[id_app])
        du_m = np.linalg.norm(du)
        dr_n = dr_n - mu_IPM * np.log([o.du_impulse_max - du_m, o.du_impulse_max - du_m, o.du_impulse_max - du_m]) \
            if o.du_impulse_max > du_m else np.array([1e10, 1e10, 1e10]) * clip(10 * (du_m - o.du_impulse_max), 0, 1)
    a_recj = dr_n if np.linalg.norm(dr_n) < np.linalg.norm(dr_p) else dr_p
    aj = tmp if np.linalg.norm(dr_n) < np.linalg.norm(dr_p) else aj
    return [aj[0], aj[1], aj[2], a_recj[0], a_recj[1], a_recj[2]]


def diff_evolve(o, T_max, id_app, interaction=True, convex=False):
    """Функция дифференциальной эволюции - неградиентный поиск скорости отталкивания/импульса   \n
    Input:                                                                                      \n
    -> o - AllObjects класс;                                                                    \n
    -> T_max - время эпизода                                                                    \n
    -> id_app - номер аппарата                                                                  \n
    -> interaction - происходит ли при импульсе отталкивание аппарата от конструкции            \n
    Output:                                                                                     \n
    -> u - вектор отталкивания/импульса (ССК при interaction=True, ОСК иначе)"""
    n = o.diff_evolve_vectors
    lst_errors = []
    u = o.X_app.v[id_app]
    a = [[np.exp(random.uniform(np.log(o.u_min), np.log(o.u_max))),
          random.uniform(-np.pi, np.pi),
          random.uniform(-np.pi / 2, np.pi / 2)] for _ in range(n)]
    a_rec = [[None, None, None] for _ in range(n)]
    tmp = [[0. for _ in range(6)] for _ in range(n)]
    tmp1_ind = [-1 for _ in range(n)]
    tmp2_ind = [-1 for _ in range(n)]
    tmp3_ind = [-1 for _ in range(n)]
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
                        [o for _ in range(n)],
                        [T_max for _ in range(n)],
                        [id_app for _ in range(n)],
                        [interaction for _ in range(n)],
                        [convex for _ in range(n)],
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
                        [n for _ in range(n)])
        else:
            for ii in range(n):
                tmp[ii] = diff_evolve_sample(o, T_max, id_app, interaction, convex,
                                             a[tmp1_ind[ii]][0], a[tmp1_ind[ii]][1], a[tmp1_ind[ii]][2],
                                             a[tmp2_ind[ii]][0], a[tmp2_ind[ii]][1], a[tmp2_ind[ii]][2],
                                             a[tmp3_ind[ii]][0], a[tmp3_ind[ii]][1], a[tmp3_ind[ii]][2],
                                             a[ii][0], a[ii][1], a[ii][2], a_rec[ii][0], a_rec[ii][1], a_rec[ii][2], i)

        a = [[tmp[ii][0], tmp[ii][1], tmp[ii][2]] for ii in range(n)]
        a_rec = [[tmp[ii][3], tmp[ii][4], tmp[ii][5]] for ii in range(n)]
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
