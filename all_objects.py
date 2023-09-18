from mylibs.calculation_functions import *
from mylibs.construction_functions import *
from mylibs.plot_functions import *
from mylibs.tiny_functions import *
from mylibs.control_function import *
from mylibs.im_sample import *
from mylibs.numerical_methods import *


class AllProblemObjects(object):
    """Класс содержит в себе абсолютно все параметры задачи и некоторые методы"""
    def __init__(self,
                 if_impulse_control=False,              # Управление импульсное
                 if_PID_control=False,                  # Управление ПД-регулятором
                 if_LQR_control=False,                  # Управление ЛКР
                 if_avoiding=False,                     # Исскуственное избежание столкновения

                 N_apparatus=1,                         # Количество аппаратов
                 diff_evolve_vectors=10,                # Количество проб дифф. эволюции
                 diff_evolve_times=3,                   # Количество эпох дифф. эволюции
                 shooting_amount_repulsion=15,          # Шаги пристрелки отталкивания
                 shooting_amount_impulse=10,            # Шаги пристрелки импульсного управления

                 diff_evolve_F=0.8,                     # Гиперпараметр дифф. эволюции
                 diff_evolve_chance=0.5,                # Гиперпараметр дифф. эволюции
                 mu_IPM=0.001,                          # Гиперпараметр дифф. эволюции
                 mu_e=0.1,

                 T_total=1e10,                          # Необязательное ограничение по времени на строительство
                 T_max=500.,                            # Максимальное время перелёта
                 T_max_hard_limit=2000.,                # Максимальное время перелёта при близости нарушении 
                 freetime=50.,                          # Время неучёта столкновения после отталкивания
                 dt=1.0,                                # Шаг по времени
                 t_reaction=10.,                        # Время между обнаружением цели и включением управления
                 time_to_be_busy=100.,                  # Время занятости между перелётами
                 u_max=0.04,                            # Максимальная скорость отталкивания
                 du_impulse_max=0.4,                    # Максимальная скорость импульса при импульсном управлении
                 w_twist=0.,
                 e_max=1e10,        # Относительное максимальная допустимое отклонение энергии (иск огр)
                 w_max=0.001,       # Максимально допустимая скорость вращения станции (искуственное ограничение)
                 V_max=0.04,        # Максимально допустимая поступательная скорость станции (искуственное ограничение)
                 R_max=50.,         # Максимально допустимое отклонение станции (искуственное ограничение)
                 j_max=45,          # Максимально допустимый след матрицы поворота S (искуственное ограничение)
                 a_pid_max=1e-5,    # Максимальное ускорение при непрерывном управлении

                 is_saving=False,               # Сохранение vedo-изображений
                 save_rate=1,                   # Итерации между сохранением vedo-изображений
                 coordinate_system='orbital',   # Система координат vedo-изображения

                 choice='3',                    # Тип конструкции
                 choice_complete=False,         # Уже собранная конструкция (для отладки)
                 floor=5,

                 if_talk=False,                 # Мне было скучно
                 if_multiprocessing=True,       # Многопроцессорность
                 if_testing_mode=False,         # Лишние принтпоинты
                 if_any_print=True,             # Любые принтпоинты

                 Radius_orbit=6800e3,           # Радиус орбиты
                 mu=5.972e24 * 6.67408e-11,     # Гравитационный параметр Земли
                 d_to_grab=0.5,                 # Расстояние захвата до цели
                 d_crash=0.1,                   # Расстояние соударения до осей стержней

                 k_p=3e-4,                      # Коэффициент ПД-регулятора
                 k_u=1e-1,                      # Коэффициент разброса скорости
                 k_av=1e-5,                     # Коэффициент поля отталкивания
                 k_ac=0.,                       # Коэффициент паразитного ускорения
                 level_avoid=2,
                 s=None,                        # готовый объект класса Structure
                 c=None,                        # готовый объект класса Container
                 a=None,                        # готовый объект класса Apparatus
                 file_reset=False,
                 method='trust-const',
                 fons_fluminis=True,  # А что делаешь ты?
                 if_T_in_shooting=False,
                 begin_rotation='xx'):

        # Параметры типа bool
        self.file_name = 'storage/main.txt'
        self.main_numerical_simulation = s is None
        self.survivor = True            # Зафиксирован ли проход "через текстуры" (можно сделать вылет программы)
        self.warning_message = True     # Если где-то проблема, вместо вылета программы я обозначаю её сообщениями
        self.t_flyby = T_max * 0.95     # Время необходимости облёта
        self.if_talk = if_talk
        self.if_multiprocessing = if_multiprocessing
        self.if_testing_mode = if_testing_mode
        self.if_any_print = if_any_print
        self.flag_impulse = True
        self.collision_foo = None

        # Параметры управления
        self.d_crash = d_crash
        self.if_impulse_control = if_impulse_control
        self.if_PID_control = if_PID_control
        self.if_LQR_control = if_LQR_control
        self.if_avoiding = if_avoiding
        self.control = if_impulse_control or if_PID_control or if_LQR_control
        self.diff_evolve_vectors = diff_evolve_vectors
        self.diff_evolve_times = diff_evolve_times
        self.shooting_amount_repulsion = shooting_amount_repulsion
        self.shooting_amount_impulse = shooting_amount_impulse
        self.diff_evolve_F = diff_evolve_F
        self.diff_evolve_chance = diff_evolve_chance
        self.mu_ipm = mu_IPM
        self.mu_e = mu_e
        self.k_p = k_p
        self.k_d = kd_from_kp(k_p)
        self.k_u = k_u
        self.k_av = k_av
        self.k_ac = k_ac
        self.level_avoid = level_avoid

        # Параметры времени
        self.T_total = T_total
        self.T_max = T_max
        self.T_max_hard_limit = T_max_hard_limit
        self.freetime = freetime
        self.t = 0.
        self.iter = 0
        self.dt = dt
        self.time_to_be_busy = time_to_be_busy
        self.t_reaction = t_reaction
        self.t_reaction_counter = t_reaction
        self.t_flyby_counter = self.t_flyby

        # Кинематические параметры и ограничения
        self.u_max = u_max
        self.u_min = u_max / 1e2
        self.du_impulse_max = du_impulse_max
        self.e_max = e_max
        self.w_max = w_max
        self.V_max = V_max
        self.R_max = R_max
        self.j_max = j_max
        self.a_pid_max = a_pid_max
        self.Radius_orbit = Radius_orbit
        self.Re = Radius_orbit
        self.d_to_grab = d_to_grab
        self.mu = mu

        # Параметры рассматриваемых объектов
        if s is None:
            self.s, self.c, self.a = get_all_components(choice=choice, complete=choice_complete, n_app=N_apparatus,
                                                        floor=floor)
            if file_reset:
                f = open(self.file_name, 'w')
                f.write(f"ограничения {self.R_max} {self.V_max} {self.j_max} {self.w_max}\n")
                f.close()
                f = open('storage/repulsions.txt', 'w')
                f.close()
                f = open('storage/iteration_docking.txt', 'w')
                f.close()
        else:
            self.s = s
            self.c = c
            self.a = a
        # self.N_cont_beams = len(self.c.mass)
        # self.N_app = N_apparatus
        self.choice = choice
        self.t_start = np.zeros(self.a.n + 1)
        self.M = np.sum(self.s.mass) + np.sum(self.c.mass)

        q12 = [[1 / np.sqrt(2) * (begin_rotation[j] in i) for i in ['xyz', 'x', 'y', 'z']] for j in range(2)]
        self.w_hkw = np.sqrt(mu / Radius_orbit ** 3)
        self.w_hkw_vec = np.array([0., 0., self.w_hkw])  # ИСК
        self.La = np.array(q_dot(q12[0], q12[1]))
        self.U, self.S, self.A, self.R_e = self.get_matrices(self.La, 0.)
        self.r_ub = np.zeros(3)
        self.v_ub = np.zeros(3)
        self.J, self.r_center = call_inertia(self, [], app_y=0)  # НЕУЧЁТ НЕСКОЛЬКИХ АППАРАТОВ
        self.J_1 = np.linalg.inv(self.J)
        if self.main_numerical_simulation:
            for i in range(self.a.n):
                self.a.r[i] = self.b_o(self.a.target[i])
        self.taken_beams = np.array([])
        self.taken_beams_p = np.array([])

        self.C_R = get_c_hkw(self.r_ub, np.zeros(3), self.w_hkw)
        self.C_r = [get_c_hkw(self.a.r[i], self.a.v[i], self.w_hkw) for i in range(self.a.n)]

        self.v_p = [np.zeros(3) for _ in range(self.a.n)]
        self.dr_p = [np.zeros(3) for _ in range(self.a.n)]
        self.a_orbital = [np.zeros(3) for _ in range(self.a.n)]  # Ускорение
        self.A_orbital = np.zeros(3)  # Ускорение
        self.a_self = [np.zeros(3) for _ in range(self.a.n)]  # Ускорение
        self.a_self_params = [None for _ in range(self.a.n)]
        self.a_wrong = np.random.rand(3)
        self.a_wrong = self.a_pid_max * k_ac * self.a_wrong / np.linalg.norm(self.a_wrong)
        self.w_twist = w_twist
        self.w = np.array([0., w_twist, 0.])
        self.w_diff = 0.
        self.Om = self.U.T @ self.w + self.w_hkw_vec
        self.e = np.zeros(3)
        self.tg_tmp = np.array([100., 0., 0.])
        self.flag_vision = [False for _ in range(self.a.n)]
        self.U_begin = self.get_potential_energy()
        self.T_begin = self.get_kinetic_energy()
        self.T = 0.
        self.E = 0.
        self.E_max = 0.

        self.fons_fluminis = fons_fluminis
        self.method = method
        self.if_T_in_shooting = if_T_in_shooting
        self.repulsion_counters = [0 for _ in range(self.a.n)]

        # Параметры отображения
        self.is_saving = is_saving
        self.save_rate = save_rate
        self.coordinate_system = coordinate_system
        self.frame_counter = 0
        self.line_app_brf = [[] for _ in range(self.a.n)]
        self.line_app_orf = [[] for _ in range(self.a.n)]
        self.line_str_orf = self.r_ub

        # Отображение параметров задачи
        if self.main_numerical_simulation:
            self.my_print(f"Масса стержней {round(float(np.sum(self.s.mass)), 2)}, масса контейнера "
                          f"{round(float(np.sum(self.c.mass)))}, всего {round(self.M, 2)}. Масса аппарата "
                          f"{round(float(self.a.mass[0]), 2)}", mode='m')

        # Выбор значений в зависимости от аргументов
        self.cases = dict({'acceleration_control': lambda v: v if np.linalg.norm(v) < self.a_pid_max else
                                                             v / np.linalg.norm(v) * self.a_pid_max * 0.95,
                           'repulse_vel_control': lambda v: (v if np.linalg.norm(v) < self.u_max else
                                                             v / np.linalg.norm(v) * self.u_max * 0.95)
                                                             if np.linalg.norm(v) > self.u_min else
                                                             v / np.linalg.norm(v) * self.u_min * 1.05,
                           'diff_vel_control': lambda a, cnd: ((a if np.linalg.norm(a) < self.u_max else
                                                                a / np.linalg.norm(a) * self.u_max * 0.95)
                                                               if np.linalg.norm(a) > self.u_min else
                                                               a / np.linalg.norm(a) * self.u_min * 1.05) if cnd
                                                                                                          else a})

    # ----------------------------------------- РАСЧЁТ ПАРАМЕТРОВ
    def get_e_deviation(self):
        if self.E_max > 1e-4:
            return self.T / self.E_max
        else:
            return self.E

    def get_discrepancy(self, id_app: int, vector: bool = False, r=None):
        """Возвращает невязку аппарата с целью"""
        tmp = self.a.r[id_app] if r is None else r
        discrepancy = tmp - self.b_o(self.a.target[id_app])
        return discrepancy if vector else np.linalg.norm(discrepancy)

    def get_angular_momentum(self):
        return np.linalg.norm(self.J @ self.S @ self.w)

    def get_kinetic_energy(self):
        """Возвращет кинетическую энергию вращения станции"""
        return self.w.T @ self.S.T @ self.J @ self.S @ self.w / 2

    def get_potential_energy(self):
        tmp = 0
        J_irf = self.b_i(self.J)
        gamma = self.R_e / self.Radius_orbit
        for i in range(3):
            for j in range(3):
                tmp += 3 * J_irf[i][j] * gamma[i] * gamma[j]
        # tmp = 3 * gamma.T @ J_orf @ gamma  # попробовать
        return 1 / 2 * self.mu / self.Radius_orbit ** 3 * (tmp + np.trace(self.J))  # self.mu/self.Radius_orbit + ()

    def get_matrices(self, La=None, t=None):
        """Функция подсчёта матриц поворота из кватернионов; \n
        На вход подаётся кватернион La и скаляр \n
        Заодно считает вектор от центра Земли до центра масс системы в ИСК."""
        La = self.La if (La is None) else La
        La /= np.linalg.norm(La)
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

    def get_hkw_acceleration(self, rv):
        return np.array([-2 * self.w_hkw * rv[5],
                         -self.w_hkw ** 2 * rv[1],
                         2 * self.w_hkw * rv[3] + 3 * self.w_hkw ** 2 * rv[2]])

    def get_ext_momentum_rigid_body(self, A, J, R_e):
        return 3 * self.mu * my_cross(A @ R_e, J @ A @ R_e) / self.Radius_orbit ** 5

    def get_masses(self, id_app: int):
        """Нет учёта нескольких аппаратов!"""
        id_beam = self.a.flag_beam[id_app]
        m_beam = 0. if (id_beam is None) else self.s.mass[id_beam]
        m_ub = self.M - m_beam
        m_a = self.a.mass[id_app] + m_beam
        return m_a, m_ub

    def get_repulsion(self, id_app: int):
        file = open('storage/repulsions_-1.txt', 'r')
        lcl_counter = self.repulsion_counters[id_app]
        anw = None
        for line in file:
            lst = line.split()
            if int(lst[1]) == id_app:
                if lcl_counter == 0:
                    anw = np.array([float(lst[2+i]) for i in range(3)])
                lcl_counter -= 1
        file.close()
        return anw

    # ----------------------------------------- РУНГЕ-КУТТЫ 4 ПОРЯДКА
    def rv_right_part(self, rv, a):
        return np.array([rv[3], rv[4], rv[5], a[0], a[1], a[2]])

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
        _, _, A, R_e = self.get_matrices(La=La, t=t)
        J1 = np.linalg.inv(J)
        M_external = self.get_ext_momentum_rigid_body(A, J, R_e)
        # M_external = np.zeros(3)
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
        return La + LaOm[0:4], LaOm[4:7] + Om

    # ----------------------------------------- ОБНОВЛЕНИЕ ПЕРЕМЕННЫХ КЛАССА
    def w_update(self):
        self.w = self.U @ (self.Om - self.w_hkw_vec)

    def om_update(self):
        self.Om = self.U.T @ self.w + self.w_hkw_vec

    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt
        self.U = self.get_potential_energy()
        self.T = self.get_kinetic_energy()  # - self.T_begin
        self.E = self.T + self.U  # - self.U_begin - self.T_begin
        self.E_max = max(self.E_max, self.E)

        # Euler rotation
        tmp = self.Om.copy()
        self.La, self.Om = self.rk4_w(self.La, self.Om, self.J, self.t)
        self.U, self.S, self.A, self.R_e = self.get_matrices()
        self.w_update()
        self.w_diff = np.linalg.norm(self.w - np.array([0., self.w_twist, 0.]))
        self.e = (self.Om - tmp) / self.dt

        # Translational movement of the structure
        self.r_ub = r_hkw(self.C_R, self.w_hkw, self.t - self.t_start[self.a.n])
        self.v_ub = v_hkw(self.C_R, self.w_hkw, self.t - self.t_start[self.a.n])
        self.A_orbital = self.get_hkw_acceleration(np.append(self.r_ub, self.v_ub))

        # Translational movement of devices
        for id_app in self.a.id:
            if self.a.flag_fly[id_app] == 1:
                if self.a.flag_hkw[id_app]:
                    r = r_hkw(self.C_r[id_app], self.w_hkw, self.t - self.t_start[id_app])
                    v = v_hkw(self.C_r[id_app], self.w_hkw, self.t - self.t_start[id_app])
                else:
                    r = self.a.r[id_app]
                    v = self.a.v[id_app]
                    self.a_orbital[id_app] = self.get_hkw_acceleration(np.append(r, v))
                    r, v = self.rk4_acceleration(r, v, self.a_self[id_app] + self.a_orbital[id_app] + self.a_wrong)
            else:
                r = self.b_o(self.a.target_p[id_app])
                v = np.zeros(3)
            self.a_orbital[id_app] = self.get_hkw_acceleration(np.append(r, v))
            self.a.r[id_app] = r
            self.a.v[id_app] = v
            
            if self.d_crash is not None:
                if self.warning_message and self.main_numerical_simulation \
                        and (self.t - self.t_start[id_app]) > self.freetime:
                    if self.survivor:
                        self.survivor = not call_crash(o=self, r_sat=r, R=self.r_ub,
                                                       S=self.S, taken_beams=self.taken_beams)
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
        if self.method == 'linear-propulsion' and self.a_self_params[id_app] is not None:
            self.a_self[id_app] = simple_control(self, self.a_self_params[id_app],
                                                 (self.t - self.t_start[id_app]) / self.T_max)
        if self.control and control_condition(o=self, id_app=id_app):
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

    def repulse_app_config(self, id_app: int):
        # Алгоритм выбора цели
        if self.a.flag_start[id_app]:  # IF ON START
            id_beam = self.s.call_possible_transport(self.taken_beams)[0]
            self.taken_beams = np.append(self.taken_beams, id_beam)
            self.my_print(f"Аппарат {id_app} забрал стержень {id_beam}", mode="b")
        else:
            if self.a.flag_beam[id_app] is None:  # GOING HOME
                id_beam = None
                self.my_print(f"Аппарат {id_app} возвращается на базу", mode="b")
            else:  # GOING TO POINT
                id_beam = self.a.flag_beam[id_app]
                self.my_print(f"Аппарат {id_app} летит установливать стержень {id_beam}", mode="b")

        if id_beam is not None:
            r_1 = self.s.r1[id_beam]
        else:
            tmp_id = self.s.call_possible_transport(self.taken_beams)[0]
            r_1 = self.s.r_st[tmp_id] - np.array([self.s.length[tmp_id], 0.0, 0.0])

        self.a.target_p[id_app] = self.a.target[id_app].copy()
        self.a.target[id_app] = r_1
        self.a.flag_beam[id_app] = id_beam
        self.a.flag_start[id_app] = False
        self.a.flag_fly[id_app] = True
        self.a.flag_hkw[id_app] = False if (self.if_PID_control or self.if_LQR_control) else True
        self.a.r_0[id_app] = self.get_discrepancy(id_app)
        self.flag_vision[id_app] = False

    def get_repulsion_change_params(self, id_app: int):
        R_p = self.r_ub.copy()
        V_p = self.v_ub.copy()
        J_p, r_center_p = call_inertia(self, self.taken_beams_p, app_y=id_app)
        r0 = np.array(self.a.target_p[id_app])
        r = self.b_o(self.a.target_p[id_app], r_center=r_center_p)
        J, r_center = call_inertia(self, self.taken_beams, app_n=id_app)
        J_1 = np.linalg.inv(self.J)
        R0c = r_center - r_center_p
        r0c = r0 - r_center_p
        R = R_p + self.S.T @ R0c
        m_extra, M_without = self.get_masses(id_app)

        return m_extra, M_without, J, J_1, J_p, r_center, r_center_p, r, R, r0c, R0c, R_p, V_p

    def repulsion_change_params(self, id_app: int, u0):
        if len(u0.shape) != 1 or len(u0) != 3:
            raise Exception("Неправильный входной вектор скорости!")
        '''R_p = self.R.copy()
        V_p = self.V.copy()
        self.J, self.r_center = call_inertia(self, self.taken_beams_p, app_y=id_app)
        J_p = self.J.copy()
        r_center_p = self.r_center.copy()
        r0 = np.array(self.a.target_p[id_app])
        r = self.b_o(self.a.target_p[id_app])
        self.J, self.r_center = call_inertia(self, self.taken_beams, app_n=id_app)
        self.J_1 = np.linalg.inv(self.J)
        R0 = self.r_center
        R0c = self.r_center - r_center_p
        r0c = r0 - r_center_p
        R = R_p + self.S.T @ R0c

        m_extra, M_without = self.get_masses(id_app)'''

        m_extra, M_without, J, J_1, J_p, r_center, r_center_p, r, R, r0c, R0c, R_p, V_p = \
            self.get_repulsion_change_params(id_app)
        u_rot = my_cross(self.w, self.S.T @ r0c)                        # ORF
        V_rot = my_cross(self.w, self.S.T @ R0c)                        # ORF
        V0 = - u0 * m_extra / M_without                                 # BRF
        u = self.S.T @ u0 + V_p + u_rot                                 # ORF
        V = self.S.T @ V0 + V_p + V_rot                                 # ORF
        self.w = self.b_o(J_1) @ (
                self.b_o(J_p) @ self.w + (M_without + m_extra) * my_cross(R_p, V_p) -
                m_extra * my_cross(r, u) - M_without * my_cross(R, V))  # ORF

        self.om_update()
        self.C_r[id_app] = get_c_hkw(r, u, self.w_hkw)
        self.C_R = get_c_hkw(R, V, self.w_hkw)
        self.t_start[id_app] = self.t
        self.t_start[self.a.n] = self.t
        self.repulsion_counters[id_app] += 1

        self.a.r[id_app] = r
        self.a.v[id_app] = u

    def capturing_change_params(self, id_app: int):
        J_p, r_center_p = call_inertia(self, self.taken_beams_p, app_n=id_app)
        m_extra, M_without = self.get_masses(id_app)
        R_p = self.r_ub.copy()
        V_p = self.v_ub.copy()
        r_p = self.a.r[id_app].copy()
        v_p = self.a.v[id_app].copy()

        id_beam = self.a.flag_beam[id_app]
        if id_beam is not None:
            if np.linalg.norm(np.array(self.s.r1[id_beam]) - np.array(self.a.target[0])) < 1e-2:
                if self.main_numerical_simulation:
                    self.my_print(f"Стержень id:{id_beam} устанавливается", mode="b")
                self.s.flag[id_beam] = np.array([1., 1.])
                self.a.flag_beam[id_app] = None
                self.taken_beams = np.delete(self.taken_beams, np.argmax(self.taken_beams == id_beam))
        else:
            if self.a.target[id_app][0] < -0.6:  # Если "слева" нет промежуточных точек, то окей
                if self.main_numerical_simulation:
                    self.my_print(f'Аппарат id:{id_app} в грузовом отсеке')
                self.a.flag_start[id_app] = True

        J, r_center = call_inertia(self, self.taken_beams, app_y=id_app)
        m, M = self.get_masses(id_app)

        self.r_ub -= self.S.T @ (r_center - r_center_p)
        self.r_ub = np.zeros(3)  # ШАМАНСТВО
        self.v_ub = (V_p * M_without + v_p * m_extra) / (M_without + m_extra)
        self.w = np.linalg.inv(self.b_o(J)) @ (self.b_o(J_p) @ self.w -
                                               (M + m) * my_cross(self.r_ub, self.v_ub) +
                                               m_extra * my_cross(r_p, v_p) + M_without * my_cross(R_p, V_p))
        self.om_update()
        self.C_R = get_c_hkw(self.r_ub, self.v_ub, self.w_hkw)
        '''if self.main_numerical_simulation:
            print(f"R:{R_p}")
            print(f"c:{self.S.T @ r_center} ~ {self.S.T @ r_center_p} = {self.S.T @ (r_center - r_center_p)}")
            print(f"R:{self.R}")
            print(f"V:{self.V}")
            print(f"C:{self.C_R}")'''
        self.a.busy_time[id_app] = self.time_to_be_busy
        self.a.flag_hkw[id_app] = True
        self.a.flag_fly[id_app] = False
        self.t_reaction_counter = self.t_reaction
        self.t_start[self.a.n] = self.t
        self.a.target_p[id_app] = self.a.target[id_app].copy()
        self.flag_impulse = True
        self.taken_beams_p = self.taken_beams.copy()

    # ----------------------------------------- ПЕРЕХОДЫ МЕЖДУ СИСТЕМАМИ КООРДИНАТ
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
        R = self.r_ub if (R is None) else R
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
        R = self.r_ub if (R is None) else R
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
        R = self.r_ub if (R is None) else R
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
        R = self.r_ub if (R is None) else R
        r_center = self.r_center if (r_center is None) else r_center
        if len(a_np.shape) == 1:
            return U.T @ ((S.T @ (a_np - r_center) + R) + np.array([0, 0, self.Re]))
        if len(a_np.shape) == 2:
            return U.T @ S.T @ a_np @ S @ U
        raise Exception("Put vector or matrix")

    # ----------------------------------------- КОСМЕТИКА
    def file_save(self, txt):
        file = open(self.file_name, 'a')
        file.write(txt + f" {int(self.main_numerical_simulation)}\n")
        file.close()

    def repulsion_save(self, txt):
        file = open('storage/repulsions.txt', 'a')
        file.write(txt + f" {int(self.main_numerical_simulation)}\n")
        file.close()

    def my_print(self, txt, mode=None, test=False):
        if (self.if_any_print and not test) or (self.if_testing_mode and test):
            if mode is None and not test:
                print(Style.RESET_ALL + txt)
            if mode == "b":
                print(Fore.BLUE + txt + Style.RESET_ALL)
            if mode == "g":
                print(Fore.GREEN + txt + Style.RESET_ALL)
            if mode == "y" or test and mode is None:
                print(Fore.YELLOW + txt + Style.RESET_ALL)
            if mode == "r":
                print(Fore.RED + txt + Style.RESET_ALL)
            if mode == "c":
                print(Fore.CYAN + txt + Style.RESET_ALL)
            if mode == "m":
                print(Fore.MAGENTA + txt + Style.RESET_ALL)

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
        slf.w_hkw_vec = copy.deepcopy(self.w_hkw_vec)
        slf.U = copy.deepcopy(self.U)
        slf.S = copy.deepcopy(self.S)
        slf.A = copy.deepcopy(self.A)
        slf.R_e = copy.deepcopy(self.R_e)

        slf.J = copy.deepcopy(self.J)
        slf.r_center = copy.deepcopy(self.r_center)
        slf.r_ub = copy.deepcopy(self.r_ub)
        slf.v_ub = copy.deepcopy(self.v_ub)
        slf.J_1 = copy.deepcopy(self.J_1)

        slf.taken_beams = copy.deepcopy(self.taken_beams)
        slf.taken_beams_p = copy.deepcopy(self.taken_beams_p)

        slf.C_R = copy.deepcopy(self.C_R)
        slf.C_r = copy.deepcopy(self.C_r)

        slf.v_p = copy.deepcopy(self.v_p)
        slf.dr_p = copy.deepcopy(self.dr_p)
        slf.a_self = copy.deepcopy(self.a_self)
        slf.a_orbital = copy.deepcopy(self.a_orbital)
        slf.A_orbital = copy.deepcopy(self.A_orbital)
        slf.w = copy.deepcopy(self.w)
        slf.w_twist = self.w_twist
        slf.w_diff = self.w_diff
        slf.Om = copy.deepcopy(self.Om)
        slf.tg_tmp = copy.deepcopy(self.tg_tmp)

        slf.E = self.E
        slf.E_max = self.E_max

        slf.method = self.method
        slf.fons_fluminis = self.fons_fluminis
        slf.if_T_in_shooting = self.if_T_in_shooting

        return slf
