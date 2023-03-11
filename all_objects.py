# Standard libraries
import matplotlib.pyplot as plt
from scipy import spatial
from datetime import datetime
import numpy as np
import pandas as pd
import random
from tqdm import tqdm

# Local libraries
from mylibs.calculation_functions import *
from mylibs.construction_functions import *
from mylibs.plot_functions import *
from mylibs.tiny_functions import *
from mylibs.control_function import *
from mylibs.im_sample import *


class AllProblemObjects(object):
    """Класс содержит в себе абсолютно все нужные параметры и некоторые методы"""
    def __init__(self,
                 if_impulse_control=False,              # Управление импульсное
                 if_PID_control=False,                  # Управление ПД-регулятором
                 if_LQR_control=False,                  # Управление ЛКР
                 if_avoiding=False,                     # Исскуственное избежание столкновения

                 N_apparatus=1,                         # Количество аппаратов
                 diff_evolve_vectors=5,                 # Количество проб дифф. эволюции
                 diff_evolve_times=3,                   # Количество эпох дифф. эволюции
                 shooting_amount_repulsion=15,          # Шаги пристрелки отталкивания
                 shooting_amount_impulse=10,            # Шаги пристрелки импульсного управления

                 diff_evolve_F=0.8,                     # Гиперпараметр дифф. эволюции
                 diff_evolve_chance=0.5,                # Гиперпараметр дифф. эволюции
                 mu_IPM=0.01,                           # Гиперпараметр дифф. эволюции
                 mu_e=0.1,

                 T_total=100000.,                       # Необязательное ограничение по времени на строительство
                 T_max=500.,                            # Максимальное время перелёта
                 T_max_hard_limit=4000.,                # Максимальное время перелёта при близости нарушении ограничений
                 dt=1.0,                                # Шаг по времени
                 t_reaction=10.,                        # Время между обнаружением цели и включением управления
                 time_to_be_busy=10.,                   # Время занятости между перелётами
                 u_max=0.2,                             # Максимальная скорость отталкивания
                 du_impulse_max=0.4,                    # Максимальная скорость импульса при импульсном управлении
                 w_max=0.0015,      # Максимально допустимая скорость вращения станции (искуственное ограничение)
                 V_max=0.1,         # Максимально допустимая поступательная скорость станции (искуственное ограничение)
                 R_max=9.,          # Максимально допустимое отклонение станции (искуственное ограничение)
                 j_max=30.,         # Максимально допустимый след матрицы поворота S (искуственное ограничение)
                 a_pid_max=0.001,   # Максимальное ускорение при непрерывном управлении

                 is_saving=False,               # Сохранение vedo-изображений
                 save_rate=1,                   # Итерации между сохранением vedo-изображений
                 coordinate_system='orbital',   # Система координат vedo-изображения

                 choice='3',                    # Тип конструкции
                 choice_complete=False,         # Уже собранная конструкция (для отладки)

                 if_talk=True,                  # Мне было скучно
                 if_multiprocessing=True,       # Многопроцессорность
                 if_testing_mode=False,         # Лишние принтпоинты
                 if_any_print=True,             # Любые принтпоинты

                 Radius_orbit=6800e3,           # Радиус орбиты
                 mu=5.972e24 * 6.67408e-11,     # Гравитационный параметр Земли
                 d_to_grab=0.5,                 # Расстояние захвата до цели
                 d_crash=0.1,                   # Расстояние соударения до осей стержней

                 k_p=1e-4,                      # Коэффициент ПД-регулятора
                 s=None,                        # готовый объект класса Structure
                 c=None,                        # готовый объект класса Container
                 a=None,                        # готовый объект класса Apparatus
                 La=np.array(q_dot([1/np.sqrt(2), 1/np.sqrt(2), 0., 0.], [1/np.sqrt(2), 0., 1/np.sqrt(2), 0.]))):

        # Инициализация переменных
        self.survivor = True            # Зафиксирован ли проход "через текстуры" / можно сделать вылет программы
        self.warning_message = None     # Если где-то проблема, вместо вылета программы я обозначаю её сообщениями
        self.t_flyby = T_max * 0.95     # Время необходимости облёта
        self.if_talk = if_talk
        self.if_multiprocessing = if_multiprocessing
        self.if_testing_mode = if_testing_mode
        self.if_any_print = if_any_print
        self.flag_impulse = True
        self.collision_foo = None

        self.d_crash = np.double(d_crash)
        self.if_impulse_control = if_impulse_control
        self.if_PID_control = if_PID_control
        self.if_LQR_control = if_LQR_control
        self.if_avoiding = if_avoiding
        self.control = True if (if_impulse_control or if_PID_control or if_LQR_control) else False

        self.diff_evolve_vectors = diff_evolve_vectors
        self.diff_evolve_times = diff_evolve_times
        self.shooting_amount_repulsion = shooting_amount_repulsion
        self.shooting_amount_impulse = shooting_amount_impulse

        self.diff_evolve_F = np.double(diff_evolve_F)
        self.diff_evolve_chance = np.double(diff_evolve_chance)
        self.mu_ipm = np.double(mu_IPM)
        self.mu_e = mu_e

        self.T_total = T_total
        self.T_max = T_max
        self.T_max_hard_limit = T_max_hard_limit
        self.t = np.double(0.)
        self.iter = 0
        self.dt = np.double(dt)
        self.time_to_be_busy = time_to_be_busy
        self.t_reaction = t_reaction
        self.t_reaction_counter = t_reaction
        self.t_flyby_counter = self.t_flyby
        self.u_max = np.double(u_max)
        self.u_min = np.double(u_max / 10)
        self.du_impulse_max = np.double(du_impulse_max)
        self.w_max = np.double(w_max)
        self.V_max = np.double(V_max)
        self.R_max = np.double(R_max)
        self.j_max = np.double(j_max)
        self.a_pid_max = np.double(a_pid_max)

        self.is_saving = is_saving
        self.save_rate = save_rate
        self.coordinate_system = coordinate_system

        self.Radius_orbit = np.double(Radius_orbit)
        self.Re = np.double(Radius_orbit)
        self.mu = np.double(mu)
        self.d_to_grab = np.double(d_to_grab)

        self.k_p = np.double(k_p)
        self.k_d = np.double(kd_from_kp(k_p))
        self.La = np.double(La)

        # self.X, self.X_cont = get_construction(choice, choice_complete)
        # self.X_app = get_apparatus(self.X, N_apparatus)
        self.choice = choice
        self.X, self.X_cont, self.X_app = get_all_components(choice)
        self.N_nodes = self.X.n_nodes
        self.N_beams = self.X.n_beams
        self.N_cont_beams = len(self.X_cont.mass)
        self.N_app = N_apparatus
        self.t_start = np.zeros(self.N_app + 1)
        self.M = np.double(np.sum(self.X.mass) + np.sum(self.X_cont.mass))

        self.w_hkw = np.sqrt(mu / Radius_orbit ** 3)
        self.W_hkw = np.double(np.array([0., 0., self.w_hkw]))
        self.U, self.S, self.A, self.R_e = self.call_rotation_matrix(La, 0.)
        self.J, self.r_center = call_inertia(self, [], app_y=0)
        self.R = np.double(np.zeros(3))
        self.V = np.double(np.zeros(3))
        self.J_1 = np.double(np.linalg.inv(self.J))
        for i in range(self.N_app):
            print(self.X_app.r)
        self.line_app = [[self.X_app.r[i][0], self.X_app.r[i][1], self.X_app.r[i][2]] for i in range(self.N_app)]
        self.line_app_orf = [[self.X_app.r[i][0], self.X_app.r[i][1], self.X_app.r[i][2]] for i in range(self.N_app)]
        self.line_str = self.R
        self.taken_beams = np.array([])
        self.taken_beams_p = np.array([])
        self.tmp_numer_frame = 0

        self.C_R = get_C_hkw(self.R, np.zeros(3), self.w_hkw)
        self.C_r = [get_C_hkw(self.X_app.r[i], self.X_app.v[i], self.w_hkw) for i in range(self.N_app)]

        self.v_p = [np.zeros(3) for _ in range(N_apparatus)]
        self.dr_p = [np.zeros(3) for _ in range(N_apparatus)]
        self.a_orbital = [np.double(np.zeros(3)) for _ in range(N_apparatus)]  # Ускорение
        self.A_orbital = np.double(np.zeros(3))  # Ускорение
        self.a_self = [np.double(np.zeros(3)) for _ in range(self.N_app)]  # Ускорение
        self.w = np.double(np.zeros(3))
        self.Om = np.double(self.U.T @ self.w + self.W_hkw)
        self.e = np.double(np.zeros(3))
        self.tg_tmp = np.double(np.array([100., 0., 0.]))
        self.flag_vision = [False for _ in range(self.N_app)]

        # Выбор значений в зависимости от аргументов
        self.cases = dict({'diff_vel_control': lambda a, cnd: ((a if np.linalg.norm(a) < self.u_max else
                                                                a / np.linalg.norm(a) * self.u_max * 0.95)
                                                               if np.linalg.norm(a) > self.u_min else
                                                               a / np.linalg.norm(a) * self.u_min * 1.05) if cnd
                                                                                                          else a})
        self.cases = dict({'repulse_vel_control': lambda v: (v if np.linalg.norm(v) < self.u_max else
                                                             v / np.linalg.norm(v) * self.u_max * 0.95)
                                                            if np.linalg.norm(v) > self.u_min else
                                                            v / np.linalg.norm(v) * self.u_min * 1.05})

    def get_discrepancy(self, id_app, vector=False):
        """Возвращает невязку аппарата с целью"""
        discrepancy = self.X_app.r[id_app] - self.b_o(self.X_app.target[id_app])
        return discrepancy if vector else np.linalg.norm(discrepancy)

    def get_angular_momentum(self):
        return np.linalg.norm(self.J @ self.S @ self.w)

    def get_kinetic_energy(self):
        """Возвращет кинетическую энергию вращения станции"""
        return self.w.T @ self.S.T @ self.J @ self.S @ self.w / 2

    def get_potential_energy(self):
        tmp = 0
        J_orf = self.b_i(self.J)
        gamma = self.R_e / self.Radius_orbit
        for i in range(3):
            for j in range(3):
                tmp += 3 * J_orf[i][j] * gamma[i] * gamma[j]
        return 1 / 2 * self.mu / self.Radius_orbit ** 3 * (tmp + np.trace(self.J))  # self.mu/self.Radius_orbit + ()

    def w_update(self):
        self.w = self.U @ np.double(self.Om - self.W_hkw)

    def Om_update(self):
        self.Om = np.double(self.U.T @ self.w + self.W_hkw)

    def call_rotation_matrix(self, La=None, t=None):
        """Функция подсчёта матриц поворота из кватернионов; \n
        На вход подаются кватернионы Lu,Ls и скаляр \n
        Заодно считает вектор от центра Земли до центра масс системы в ИСК."""
        La = self.La if (La is None) else La
        t = self.t if (t is None) else t
        A = quart2dcm(La)
        U = np.array([[0., 1., 0.],
                      [0., 0., 1.],
                      [1., 0., 0.]]) @ \
            np.array([[np.cos(t * self.w_hkw), np.sin(t * self.w_hkw), 0],
                      [-np.sin(t * self.w_hkw), np.cos(t * self.w_hkw), 0],
                      [0, 0, 1]])
        S = np.double(A @ U.T)
        R_e = np.double(U.T @ np.array([0, 0, self.Radius_orbit]))
        return U, S, A, R_e

    def orbital_acceleration(self, rv):
        return np.array([-2 * self.w_hkw * rv[5],
                         -self.w_hkw ** 2 * rv[1],
                         2 * self.w_hkw * rv[3] + 3 * self.w_hkw ** 2 * rv[2]])

    def rv_right_part(self, rv, a):
        a_orbit = self.orbital_acceleration(rv)
        return np.array([rv[3], rv[4], rv[5], a_orbit[0] + a[0], a_orbit[1] + a[1], a_orbit[2] + a[2]])

    def rk4_acceleration(self, r, v, a):
        rv = np.append(r, v)
        k1 = self.rv_right_part(rv, a)
        k2 = self.rv_right_part(rv + k1 * self.dt / 2, a)
        k3 = self.rv_right_part(rv + k2 * self.dt / 2, a)
        k4 = self.rv_right_part(rv + k3 * self.dt, a)
        rv = self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return np.double(rv[0:3] + r), np.double(rv[3:6] + v)

    def lw_right_part(self, LaOm, t, J):
        """Функция правых частей для угловой скорости; \n
        Используется в методе Рунге-Кутты."""
        La, Om = LaOm[0:4], LaOm[4:7]
        dLa = np.double(1 / 2 * q_dot([0, Om[0], Om[1], Om[2]], La))
        U, S, A, R_e = self.call_rotation_matrix(La=La + dLa, t=t)
        J1 = np.linalg.inv(J)
        M_external = np.double(3 * self.mu * my_cross(A @ R_e, J @ A @ R_e) / self.Radius_orbit ** 5)
        return np.append(dLa, A.T @ J1 @ (M_external - my_cross(A @ Om, J @ A @ Om)))

    def rk4_w(self, La, Om, J, t):
        """Интегрирование уравнение Эйлера методом Рунге-Кутты 4 порядка; \n
        Используется в виде: \n
        La, Om = rk4_w(La, Om, J, t)."""
        LaOm = np.append(La, Om)  # Запихиваем в один вектор
        k1 = self.lw_right_part(LaOm, t, J)
        k2 = self.lw_right_part(LaOm + k1 * self.dt / 2, t + self.dt / 2, J)
        k3 = self.lw_right_part(LaOm + k2 * self.dt / 2, t + self.dt / 2, J)
        k4 = self.lw_right_part(LaOm + k3 * self.dt, t + self.dt, J)
        LaOm = self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return (LaOm[0:4] + La) / np.linalg.norm(LaOm[0:4] + La), LaOm[4:7] + Om

    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt

        # Вращение конструкции по Эйлеру
        self.La, self.Om = self.rk4_w(self.La, self.Om, self.J, self.t)
        self.U, self.S, self.A, self.R_e = self.call_rotation_matrix()
        tmp = np.double(self.w.copy())
        self.w_update()
        self.e = np.double((self.w - tmp) / self.dt)

        # Поступательное движение конструкции
        self.R = r_HKW(self.C_R, self.mu, self.w_hkw, self.t - self.t_start[self.N_app])
        self.V = v_HKW(self.C_R, self.mu, self.w_hkw, self.t - self.t_start[self.N_app])
        self.A_orbital = self.orbital_acceleration(np.append(self.R, self.V))

        # Поступательное движение аппаратов
        for id_app in self.X_app.id:
            if self.X_app.flag_fly[id_app] == 1:
                if self.X_app.flag_hkw[id_app]:
                    r = r_HKW(self.C_r[id_app], self.mu, self.w_hkw, self.t - self.t_start[id_app])
                    v = v_HKW(self.C_r[id_app], self.mu, self.w_hkw, self.t - self.t_start[id_app])
                else:
                    r = self.X_app.r[id_app]
                    v = self.X_app.v[id_app]
                    r, v = self.rk4_acceleration(r, v, self.a_self[id_app])
            else:
                r = self.b_o(self.X_app.target[id_app])
                v = np.zeros(3)
            self.X_app.r[id_app] = r
            self.X_app.v[id_app] = v
            self.a_orbital[id_app] = self.orbital_acceleration(np.append(r, v))

    def control_step(self, id_app):
        """ Функция ускорения бортового двигателя / подачи импульса двигателя.
        Вызывается на каждой итерации.
        :param id_app: id-номер аппарата
        :return: None
        """
        if control_condition(o=self, id_app=id_app):
            if self.if_impulse_control:
                impulse_control(o=self, id_app=id_app)
            if self.if_PID_control and self.t_reaction_counter < 0:
                pd_control(o=self, id_app=id_app)
            if self.if_LQR_control and self.t_reaction_counter < 0:
                lqr_control(o=self, id_app=id_app)
            if self.if_avoiding:
                self.a_self[id_app] += avoiding_force(self, id_app)
                if np.linalg.norm(self.a_self[id_app]) > self.a_pid_max:
                    self.a_self[id_app] *= self.a_pid_max / np.linalg.norm(self.a_self[id_app])

    def i_o(self, a, U=None):
        """Инерциальная -> Орбитальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        if len(a_np.shape) == 1:
            return np.double(U @ a_np - np.array([0, 0, self.Re]))
        if len(a_np.shape) == 2:
            return np.double(U @ a_np @ U.T)
        raise Exception("Put vector or matrix")

    def o_i(self, a, U=None):
        """Орбитальная -> Инерциальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        if len(a_np.shape) == 1:
            return np.double(U.T @ (a_np + np.array([0, 0, self.Re])))
        if len(a_np.shape) == 2:
            return np.double(U.T @ a_np @ U)
        raise Exception("Put vector or matrix")

    def o_b(self, a, S=None, R=None, r_center=None):
        """Орбитальная -> Связная"""
        a_np = np.array(a)
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return np.double(S @ (a_np - R) + r_center)
        if len(a_np.shape) == 2:
            return np.double(S @ a_np @ S.T)
        raise Exception("Put vector or matrix")

    def b_o(self, a, S=None, R=None, r_center=None):
        """Связная -> Орбитальная"""
        a_np = np.array(a)
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return np.double(S.T @ (a_np - r_center) + R)
        if len(a_np.shape) == 2:
            return np.double(S.T @ a_np @ S)
        raise Exception("Put vector or matrix")

    def i_b(self, a, U=None, S=None, R=None, r_center=None):
        """Инерциальная -> Связная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return np.double(S @ ((U @ a_np - np.array([0, 0, self.Re])) - R) + r_center)
        if len(a_np.shape) == 2:
            return np.double(S @ U @ a_np @ U.T @ S.T)
        raise Exception("Put vector or matrix")

    def b_i(self, a, U=None, S=None, R=None, r_center=None):
        """Связная -> Инерциальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return np.double(U.T @ ((S.T @ (a_np - r_center) + R) + np.array([0, 0, self.Re])))
        if len(a_np.shape) == 2:
            return np.double(U.T @ S.T @ a_np @ S @ U)
        raise Exception("Put vector or matrix")

    def file_save(self, f, txt):
        if self.is_saving:
            f.write(txt)

    def my_print(self, txt, mode=None, test=False):
        if (self.if_any_print and not test) or (self.if_testing_mode and test):
            if mode is None:
                print(Style.RESET_ALL + txt)
            if mode == "b":
                print(Fore.BLUE + txt)
            if mode == "g":
                print(Fore.GREEN + txt)
            if mode == "y":
                print(Fore.YELLOW + txt)
            if mode == "r":
                print(Fore.RED + txt)
            if mode == "c":
                print(Fore.CYAN + txt)

    def copy(self):
        slf = AllProblemObjects(choice=self.choice)

        slf.survivor = self.survivor
        slf.warning_message = False
        slf.t_flyby = self.t_flyby
        slf.if_talk = False
        slf.if_multiprocessing = self.if_multiprocessing
        slf.if_testing_mode = self.if_testing_mode
        slf.if_any_print = self.if_any_print
        slf.flag_impulse = self.flag_impulse
        slf.collision_foo = None

        slf.d_crash = self.d_crash
        slf.if_impulse_control = self.if_impulse_control
        slf.if_PID_control = self.if_PID_control
        slf.if_LQR_control = self.if_LQR_control
        slf.if_avoiding = self.if_avoiding
        slf.control = self.control

        slf.diff_evolve_vectors = self.diff_evolve_vectors
        slf.diff_evolve_times = self.diff_evolve_times
        slf.shooting_amount_repulsion = self.shooting_amount_repulsion
        slf.shooting_amount_impulse = self.shooting_amount_impulse

        slf.diff_evolve_F = self.diff_evolve_F
        slf.diff_evolve_chance = self.diff_evolve_chance
        slf.mu_ipm = self.mu_ipm
        slf.mu_e = self.mu_e

        slf.T_total = self.T_total
        slf.T_max = self.T_max
        slf.T_max_hard_limit = self.T_max_hard_limit
        slf.iter = self.iter
        slf.t = self.t
        slf.dt = self.dt
        slf.time_to_be_busy = self.time_to_be_busy
        slf.t_reaction = self.t_reaction
        slf.t_reaction_counter = self.t_reaction_counter
        slf.t_flyby_counter = self.t_flyby
        slf.u_max = self.u_max
        slf.u_min = self.u_min
        slf.du_impulse_max = self.du_impulse_max
        slf.w_max = self.w_max
        slf.V_max = self.V_max
        slf.R_max = self.R_max
        slf.j_max = self.j_max
        slf.a_pid_max = self.a_pid_max

        slf.is_saving = self.is_saving
        slf.save_rate = self.save_rate
        slf.coordinate_system = self.coordinate_system

        slf.Radius_orbit = self.Radius_orbit
        slf.Re = self.Re
        slf.mu = self.mu
        slf.d_to_grab = self.d_to_grab

        slf.k_p = self.k_p
        slf.k_d = self.k_d
        slf.La = self.La.copy()

        slf.X = self.X.copy()
        slf.X_cont = self.X_cont.copy()
        slf.N_nodes = self.N_nodes
        slf.N_beams = self.N_beams
        slf.N_cont_beams = self.N_cont_beams
        for i in range(self.N_app):
            self.X_app.r = np.array(self.X_app.r)
            self.X_app.v = np.array(self.X_app.v)
        slf.X_app = self.X_app.copy()
        slf.N_app = self.N_app
        slf.t_start = self.t_start.copy()
        slf.M = self.M

        slf.w_hkw = self.w_hkw
        slf.W_hkw = self.W_hkw.copy()
        slf.U = self.U.copy()
        slf.S = self.S.copy()
        slf.A = self.A.copy()
        slf.R_e = self.R_e.copy()

        slf.J = self.J.copy()
        slf.r_center = self.r_center.copy()
        slf.R = self.R.copy()
        slf.V = self.V.copy()
        slf.J_1 = self.J_1.copy()

        slf.line_str = slf.R
        slf.taken_beams = self.taken_beams.copy()
        slf.taken_beams_p = self.taken_beams_p.copy()

        slf.C_R = self.C_R.copy()
        slf.C_r = self.C_r.copy()

        slf.v_p = self.v_p.copy()
        slf.dr_p = self.dr_p.copy()
        slf.a_orbital = self.a_orbital.copy()
        slf.A_orbital = self.A_orbital.copy()
        slf.a_self = self.a_self.copy()
        slf.w = self.w.copy()
        slf.Om = self.Om.copy()
        slf.e = self.e.copy()
        slf.tg_tmp = self.tg_tmp.copy()
        slf.flag_vision = self.flag_vision.copy()

        return slf
