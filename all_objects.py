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
                 k_u=1e-1,                      # Коэффициент разброса скорости
                 s=None,                        # готовый объект класса Structure
                 c=None,                        # готовый объект класса Container
                 a=None,                        # готовый объект класса Apparatus
                 file_reset=False,
                 La=np.array(q_dot([1/np.sqrt(2), 1/np.sqrt(2), 0., 0.], [1/np.sqrt(2), 0., 1/np.sqrt(2), 0.]))):

        # Init
        self.file_name = 'storage/main.txt'
        self.main_numerical_simulation = True if s is None else False
        self.survivor = True            # Зафиксирован ли проход "через текстуры" / можно сделать вылет программы
        self.warning_message = True     # Если где-то проблема, вместо вылета программы я обозначаю её сообщениями
        self.t_flyby = T_max * 0.95     # Время необходимости облёта
        self.if_talk = if_talk
        self.if_multiprocessing = if_multiprocessing
        self.if_testing_mode = if_testing_mode
        self.if_any_print = if_any_print
        self.flag_impulse = True
        self.collision_foo = None

        self.d_crash = d_crash
        self.if_impulse_control = if_impulse_control
        self.if_PID_control = if_PID_control
        self.if_LQR_control = if_LQR_control
        self.if_avoiding = if_avoiding
        self.control = True if (if_impulse_control or if_PID_control or if_LQR_control) else False

        self.diff_evolve_vectors = diff_evolve_vectors
        self.diff_evolve_times = diff_evolve_times
        self.shooting_amount_repulsion = shooting_amount_repulsion
        self.shooting_amount_impulse = shooting_amount_impulse

        self.diff_evolve_F = diff_evolve_F
        self.diff_evolve_chance = diff_evolve_chance
        self.mu_ipm = mu_IPM
        self.mu_e = mu_e

        self.T_total = T_total
        self.T_max = T_max
        self.T_max_hard_limit = T_max_hard_limit
        self.t = 0.
        self.iter = 0
        self.dt = dt
        self.time_to_be_busy = time_to_be_busy
        self.t_reaction = t_reaction
        self.t_reaction_counter = t_reaction
        self.t_flyby_counter = self.t_flyby
        self.u_max = u_max
        self.u_min = u_max / 10
        self.du_impulse_max = du_impulse_max
        self.w_max = w_max
        self.V_max = V_max
        self.R_max = R_max
        self.j_max = j_max
        self.a_pid_max = a_pid_max

        self.is_saving = is_saving
        self.save_rate = save_rate
        self.coordinate_system = coordinate_system

        self.Radius_orbit = Radius_orbit
        self.Re = Radius_orbit
        self.mu = mu
        self.d_to_grab = d_to_grab

        self.k_p = k_p
        self.k_d = kd_from_kp(k_p)
        self.k_u = k_u
        self.La = La

        self.choice = choice
        if s is None:
            self.s, self.c, self.a = get_all_components(choice=choice, complete=choice_complete, n_app=N_apparatus)
            if file_reset:
                f = open(self.file_name, 'w')
                f.close()
            else:
                f = open(self.file_name, 'a')
                f.write("--------------------------------------\n")
                f.close()
        else:
            self.s = s
            self.c = c
            self.a = a
        self.N_nodes = self.s.n_nodes
        self.N_beams = self.s.n_beams
        self.N_cont_beams = len(self.c.mass)
        self.N_app = N_apparatus
        self.t_start = np.zeros(self.N_app + 1)
        self.M = np.sum(self.s.mass) + np.sum(self.c.mass)

        self.w_hkw = np.sqrt(mu / Radius_orbit ** 3)
        self.W_hkw = np.array([0., 0., self.w_hkw])
        self.U, self.S, self.A, self.R_e = self.call_rotation_matrix(La, 0.)
        self.J, self.r_center = call_inertia(self, [], app_y=0)
        self.R = np.zeros(3)
        self.V = np.zeros(3)
        self.J_1 = np.linalg.inv(self.J)
        self.line_app = [[] for _ in range(self.N_app)]
        self.line_app_orf = [[] for _ in range(self.N_app)]
        self.line_str = self.R
        self.taken_beams = np.array([])
        self.taken_beams_p = np.array([])
        self.tmp_numer_frame = 0

        self.C_R = get_c_hkw(self.R, np.zeros(3), self.w_hkw)
        self.C_r = [get_c_hkw(self.a.r[i], self.a.v[i], self.w_hkw) for i in range(self.N_app)]

        self.v_p = [np.zeros(3) for _ in range(N_apparatus)]
        self.dr_p = [np.zeros(3) for _ in range(N_apparatus)]
        self.a_orbital = [np.zeros(3) for _ in range(N_apparatus)]  # Ускорение
        self.A_orbital = np.zeros(3)  # Ускорение
        self.a_self = [np.zeros(3) for _ in range(self.N_app)]  # Ускорение
        self.w = np.zeros(3)
        self.Om = self.U.T @ self.w + self.W_hkw
        self.e = np.zeros(3)
        self.tg_tmp = np.array([100., 0., 0.])
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

    def get_discrepancy(self, id_app: int, vector: bool = False):
        """Возвращает невязку аппарата с целью"""
        discrepancy = self.a.r[id_app] - self.b_o(self.a.target[id_app])
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
        # tmp = 3 * gamma.T @ J_orf @ gamma  # попробовать
        return 1 / 2 * self.mu / self.Radius_orbit ** 3 * (tmp + np.trace(self.J))  # self.mu/self.Radius_orbit + ()

    def w_update(self):
        self.w = self.U @ (self.Om - self.W_hkw)

    def om_update(self):
        self.Om = self.U.T @ self.w + self.W_hkw

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
        S = A @ U.T
        R_e = U.T @ np.array([0, 0, self.Radius_orbit])
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
        return rv[0:3] + r, rv[3:6] + v

    def lw_right_part(self, LaOm, t, J):
        """Функция правых частей для угловой скорости; \n
        Используется в методе Рунге-Кутты."""
        La, Om = LaOm[0:4], LaOm[4:7]
        dLa = 1 / 2 * q_dot([0, Om[0], Om[1], Om[2]], La)
        U, S, A, R_e = self.call_rotation_matrix(La=La + dLa, t=t)
        J1 = np.linalg.inv(J)
        M_external = 3 * self.mu * my_cross(A @ R_e, J @ A @ R_e) / self.Radius_orbit ** 5
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

        # Euler rotation
        self.La, self.Om = self.rk4_w(self.La, self.Om, self.J, self.t)
        self.U, self.S, self.A, self.R_e = self.call_rotation_matrix()
        tmp = self.w.copy()
        self.w_update()
        self.e = (self.w - tmp) / self.dt

        # Translational movement of the structure
        self.R = r_hkw(self.C_R, self.w_hkw, self.t - self.t_start[self.N_app])
        self.V = v_hkw(self.C_R, self.w_hkw, self.t - self.t_start[self.N_app])
        self.A_orbital = self.orbital_acceleration(np.append(self.R, self.V))

        # Translational movement of devices
        for id_app in self.a.id:
            if self.a.flag_fly[id_app] == 1:
                if self.a.flag_hkw[id_app]:
                    r = r_hkw(self.C_r[id_app], self.w_hkw, self.t - self.t_start[id_app])
                    v = v_hkw(self.C_r[id_app], self.w_hkw, self.t - self.t_start[id_app])
                else:
                    r = self.a.r[id_app]
                    v = self.a.v[id_app]
                    r, v = self.rk4_acceleration(r, v, self.a_self[id_app])
            else:
                r = self.b_o(self.a.target[id_app])
                v = np.zeros(3)
            self.a.r[id_app] = r
            self.a.v[id_app] = v
            self.a_orbital[id_app] = self.orbital_acceleration(np.append(r, v))

            if self.warning_message and self.main_numerical_simulation and (self.t - self.t_start[id_app]) > 10:
                if self.survivor:
                    self.survivor = call_crash(o=self, r_sat=r, R=self.R, S=self.S, taken_beams=self.taken_beams)
                    if not self.survivor:
                        self.my_print(get_angry_message(), mode='y')
        if not self.survivor and self.warning_message and self.iter % 200 == 0 and self.if_testing_mode:
            self.my_print(get_angry_message(), mode='r')

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

    def get_masses(self, id_app: int):
        id_beam = self.a.flag_beam[id_app]
        m_beam = 0. if (id_beam is None) else self.s.mass[id_beam]
        M_without = self.M - m_beam
        m_extra = self.a.mass[id_app] + m_beam
        return m_extra, M_without

    def repulsion_change_params(self, id_app: int, u0):
        if len(u0.shape) != 1 or len(u0) != 3:
            raise Exception("Неправильный входной вектор скорости!")
        # if self.k_u is not None and not self.main_numerical_simulation:
        #     u0 = velocity_spread(u0, self.k_u)
        R_p = self.R.copy()
        V_p = self.V.copy()
        J_p, r_center_p = call_inertia(self, self.taken_beams_p, app_y=id_app)
        self.J, self.r_center = call_inertia(self, self.taken_beams, app_n=id_app)
        self.J_1 = np.linalg.inv(self.J)
        R = R_p + self.S.T @ (self.r_center - r_center_p)
        r = np.array(self.a.r[id_app])
        r0 = self.o_b(r, R=R_p, r_center=r_center_p)

        m_extra, M_without = self.get_masses(id_app)
        u_rot = my_cross(self.w, self.S.T @ r0)                         # ORF
        V_rot = my_cross(self.w, self.S.T @ (self.r_center - r_center_p))  # ORF
        V0 = - u0 * m_extra / M_without                           # BRF
        u = V_p + self.S.T @ u0 + u_rot                              # ORF
        V = V_p + self.S.T @ V0 + V_rot                              # ORF
        self.w = self.b_o(self.J_1) @ (
                self.b_o(J_p) @ self.w - m_extra * my_cross(r, u) +
                (M_without + m_extra) * my_cross(R_p, V_p) - M_without * my_cross(R, V))  # ORF
        self.om_update()
        self.C_r[id_app] = get_c_hkw(r, u, self.w_hkw)
        self.C_R = get_c_hkw(R, V, self.w_hkw)
        self.t_start[id_app] = self.t
        self.t_start[self.a.n] = self.t

        self.a.r[id_app] = r
        self.a.v[id_app] = u

    def capturing_change_params(self, id_app: int):
        J_p, r_center_p = call_inertia(self, self.taken_beams_p, app_n=id_app)
        J, r_center = call_inertia(self, self.taken_beams, app_y=id_app)
        R_p = self.R.copy()
        V_p = self.V.copy()
        r_p = self.a.r[id_app].copy()
        v_p = self.a.v[id_app].copy()
        m_extra, M_without = self.get_masses(id_app)
        self.R += self.S.T @ (r_center - r_center_p)
        self.V = (V_p * M_without + v_p * m_extra) / (M_without + m_extra)
        self.w = np.linalg.inv(self.b_o(J)) @ (self.b_o(J_p) @ self.w -
                                               (M_without + m_extra) * my_cross(self.R, self.V) +
                                               m_extra * my_cross(r_p, v_p) + M_without * my_cross(R_p, V_p))
        self.om_update()
        self.C_R = get_c_hkw(self.R, self.V, self.w_hkw)
        self.a.busy_time[id_app] = self.time_to_be_busy
        self.a.flag_hkw[id_app] = True
        self.a.flag_fly[id_app] = False
        self.t_reaction_counter = self.t_reaction
        self.t_start[self.a.n] = self.t
        self.flag_impulse = True
        self.taken_beams_p = self.taken_beams.copy()

    def i_o(self, a, U=None):
        """Инерциальная -> Орбитальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        if len(a_np.shape) == 1:
            return U @ a_np - np.array([0, 0, self.Re])
        if len(a_np.shape) == 2:
            return U @ a_np @ U.T
        raise Exception("Put vector or matrix")

    def o_i(self, a, U=None):
        """Орбитальная -> Инерциальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        if len(a_np.shape) == 1:
            return U.T @ (a_np + np.array([0, 0, self.Re]))
        if len(a_np.shape) == 2:
            return U.T @ a_np @ U
        raise Exception("Put vector or matrix")

    def o_b(self, a, S=None, R=None, r_center=None):
        """Орбитальная -> Связная"""
        a_np = np.array(a)
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return S @ (a_np - R) + r_center
        if len(a_np.shape) == 2:
            return S @ a_np @ S.T
        raise Exception("Put vector or matrix")

    def b_o(self, a, S=None, R=None, r_center=None):
        """Связная -> Орбитальная"""
        a_np = np.array(a)
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return S.T @ (a_np - r_center) + R
        if len(a_np.shape) == 2:
            return S.T @ a_np @ S
        raise Exception("Put vector or matrix")

    def i_b(self, a, U=None, S=None, R=None, r_center=None):
        """Инерциальная -> Связная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return S @ ((U @ a_np - np.array([0, 0, self.Re])) - R) + r_center
        if len(a_np.shape) == 2:
            return S @ U @ a_np @ U.T @ S.T
        raise Exception("Put vector or matrix")

    def b_i(self, a, U=None, S=None, R=None, r_center=None):
        """Связная -> Инерциальная"""
        a_np = np.array(a)
        U = self.U if (U is None) else U
        S = self.S if (S is None) else S
        R = self.R if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return U.T @ ((S.T @ (a_np - r_center) + R) + np.array([0, 0, self.Re]))
        if len(a_np.shape) == 2:
            return U.T @ S.T @ a_np @ S @ U
        raise Exception("Put vector or matrix")

    def file_save(self, txt):
        file = open(self.file_name, 'a')
        file.write(txt + f" {int(self.main_numerical_simulation)}\n")
        file.close()

    def my_print(self, txt, mode=None, test=False):
        if (self.if_any_print and not test) or (self.if_testing_mode and test):
            if mode is None and not test:
                print(Style.RESET_ALL + txt)
            if mode == "b":
                print(Fore.BLUE + txt)
            if mode == "g":
                print(Fore.GREEN + txt)
            if mode == "y" or test and mode is None:
                print(Fore.YELLOW + txt)
            if mode == "r":
                print(Fore.RED + txt)
            if mode == "c":
                print(Fore.CYAN + txt)
            if mode == "m":
                print(Fore.MAGENTA + txt)

    def copy(self):
        slf = AllProblemObjects(choice=self.choice, s=self.s.copy(), c=self.c.copy(), a=self.a.copy())

        slf.main_numerical_simulation = False
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
        slf.La = copy.deepcopy(self.La)

        slf.t_start = copy.deepcopy(self.t_start)
        slf.M = self.M

        slf.w_hkw = self.w_hkw
        slf.W_hkw = copy.deepcopy(self.W_hkw)
        slf.U = copy.deepcopy(self.U)
        slf.S = copy.deepcopy(self.S)
        slf.A = copy.deepcopy(self.A)
        slf.R_e = copy.deepcopy(self.R_e)

        slf.J = copy.deepcopy(self.J)
        slf.r_center = copy.deepcopy(self.r_center)
        slf.R = copy.deepcopy(self.R)
        slf.V = copy.deepcopy(self.V)
        slf.J_1 = copy.deepcopy(self.J_1)

        slf.taken_beams = copy.deepcopy(self.taken_beams)
        slf.taken_beams_p = copy.deepcopy(self.taken_beams_p)

        slf.C_R = copy.deepcopy(self.C_R)
        slf.C_r = copy.deepcopy(self.C_r)

        slf.v_p = copy.deepcopy(self.v_p)
        slf.dr_p = copy.deepcopy(self.dr_p)
        slf.a_orbital = copy.deepcopy(self.a_orbital)
        slf.A_orbital = copy.deepcopy(self.A_orbital)
        slf.w = copy.deepcopy(self.w)
        slf.Om = copy.deepcopy(self.Om)
        slf.tg_tmp = copy.deepcopy(self.tg_tmp)

        return slf
