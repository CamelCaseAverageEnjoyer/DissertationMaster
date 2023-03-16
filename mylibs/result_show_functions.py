from all_objects import *
from vedo import *
from datetime import datetime


def pd_control_params_search(dt=1., N=5, T_max=2000., k_min=1e-4, k_max=1e-2):
    filename = 'storage/pid_koeff.txt'
    f = open(filename, 'w')
    f.close()
    k_p_list = np.exp(np.linspace(np.log(k_min), np.log(k_max), N))  # Logarithmic scale
    tolerance_list = []
    start_time = datetime.now()
    tmp_count = 0
    collide = False

    for k_p in k_p_list:
        tmp_count += 1
        id_app = 0
        tolerance = None
        print(Fore.CYAN + f'Подбор ПД-к-в: {tmp_count}/{N}; время={datetime.now() - start_time}' + Style.RESET_ALL)
        o = AllProblemObjects(if_PID_control=True,
                              dt=dt, k_p=k_p,
                              T_max=T_max,
                              if_talk=False,
                              if_any_print=False,
                              choice='3')

        for i_time in range(int(T_max // dt)):
            # Repulsion
            o.a.loc[id_app, 'busy_time'] -= o.dt if o.a.busy_time[id_app] >= 0 else 0
            if o.a.flag_fly[id_app] == 0 and o.a.busy_time[id_app] < 0:
                _ = repulsion(o, id_app, u_a_priori=np.array([0.00749797, 0.00605292, 0.08625441]))

            o.control_step(id_app)

            target_orf = o.b_o(o.a.target[id_app])
            tmp = (np.linalg.norm(target_orf - np.array(o.a.r[id_app])))
            collide = collide or call_crash(o, o.a.r[id_app], o.R, o.S, o.taken_beams)
            if tolerance is None or tmp < tolerance:
                tolerance = -1 if collide else tmp
        tolerance_list.append(tolerance)
        f = open(filename, 'a')
        f.write(f'{k_p} {tolerance}\n')
        f.close()
    plt.title("Подбор коэффициентов ПД-регулятора")
    plt.plot(k_p_list, tolerance_list, c='#009ACD', label='невязка, м')
    plt.scatter(k_p_list, tolerance_list, c='#009ACD')
    plt.legend()
    plt.show()
    return k_p_list[np.argmin(tolerance_list)]


def plot_params_while_main(filename: str):
    f = open('storage/main.txt', 'r')
    o = AllProblemObjects()

    id_max = 0
    for line in f:
        lst = line.split()
        if lst[0] == 'график':
            id_max = max(id_max, 1 + int(lst[1]))
    f.close()

    def params_reset():
        return [[[] for _ in range(id_max)] for _ in range(8)]

    dr, w, j, V, R, t, a, m = params_reset()

    f = open('storage/main.txt', 'r')
    tmp = 0
    for line in f:
        lst = line.split()
        if len(lst) > 0 and tmp % 5 == 0:
            if lst[0] == 'график' and len(lst) == 9:
                id_app = int(lst[1])
                dr[id_app].append(float(lst[2]))
                w[id_app].append(float(lst[3]))
                j[id_app].append(float(lst[4]))
                V[id_app].append(float(lst[5]))
                R[id_app].append(float(lst[6]))
                a[id_app].append(float(lst[7]))
                m[id_app].append(int(lst[8]))
        tmp += 1
    f.close()
    print(Fore.CYAN + f"Аппаратов на графиках: {id_max}" + Style.RESET_ALL)

    fig, axs = plt.subplots(3)
    axs[0].set_xlabel('время, с')
    axs[0].set_ylabel('Невязка, м')
    axs[0].set_title('Параметры в процессе алгоритма')
    axs[1].set_xlabel('время t, с')
    axs[1].set_ylabel('ограниченные величины')
    axs[2].set_xlabel('итерации')
    axs[2].set_ylabel('бортовое ускорение, м/с2')

    clr = ['c', 'indigo', 'm', 'violet', 'teal', 'slategray', 'greenyellow', 'sienna']
    for id_app in range(id_max):
        t[id_app] = np.linspace(0, len(dr[id_app]), len(dr[id_app]))
        for i in range(len(dr[id_app]) - 1):
            axs[0].plot([t[id_app][i], t[id_app][i+1]], np.array([dr[id_app][i], dr[id_app][i+1]]),
                        c=clr[2 * id_app + 2 * m[id_app][i]])
        axs[1].plot(t[id_app], [1 for _ in range(len(t[id_app]))], c='gray')
        axs[2].plot(range(len(a[id_app])), a[id_app], c='c')
        axs[0].plot(t[id_app], np.zeros(len(t[id_app])), c='khaki')
        axs[2].plot(range(len(a[id_app])), np.zeros(len(a[id_app])), c='khaki')
    id_app = 0
    clr = [['skyblue', 'bisque', 'palegreen', 'darksalmon'], ['teal', 'tan', 'g', 'brown']]
    axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(w[id_app][0]) / o.w_max, np.array(w[id_app][1]) /
                                               o.w_max], c=clr[1][0], label='w')
    axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(j[id_app][0]) / o.j_max, np.array(j[id_app][1]) /
                                               o.j_max], c=clr[1][1], label='угол')
    axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(V[id_app][0]) / o.V_max, np.array(V[id_app][1]) /
                                               o.V_max], c=clr[1][2], label='V')
    axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(R[id_app][0]) / o.R_max, np.array(R[id_app][1]) /
                                               o.R_max], c=clr[1][3], label='R')
    for i in range(len(t[id_app]) - 1):
        axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(w[id_app][i]) / o.w_max, np.array(w[id_app][i+1]) /
                                                     o.w_max], c=clr[m[id_app][i]][0])
        axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(j[id_app][i]) / o.j_max, np.array(j[id_app][i+1]) /
                                                     o.j_max], c=clr[m[id_app][i]][1])
        axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(V[id_app][i]) / o.V_max, np.array(V[id_app][i+1]) /
                                                     o.V_max], c=clr[m[id_app][i]][2])
        axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(R[id_app][i]) / o.R_max, np.array(R[id_app][i+1]) /
                                                     o.R_max], c=clr[m[id_app][i]][3])

    axs[1].legend()

    if filename != '[Название]':
        plt.savefig(f"add/{filename}.jpg")
    plt.show()


def plot_a_avoid(x_boards: list = np.array([-10, 5]), z_boards: list = np.array([3, 10])):
    o = AllProblemObjects(choice='3')

    # x_boards: list = [-15, 15], z_boards: list = [1, 10]
    nx = 30
    nz = 15
    x_list = np.linspace(x_boards[0], x_boards[1], nx)
    z_list = np.linspace(z_boards[0], z_boards[1], nz)

    arrows = []
    i = 0
    forces = [np.zeros(3) for i in range(nx*nz)]
    max_force = 0
    for x in x_list:
        for z in z_list:
            tmp = avoiding_force(o, 0, r=[x, 0, z])
            if tmp is not False:
                forces[i] = tmp
                max_force = max(max_force, np.linalg.norm(tmp))
            i += 1
    i = 0
    for x in x_list:
        for z in z_list:
            force = forces[i] / max_force
            i += 1
            if force is not False:
                print(f"force: {force}")
                l1 = [np.array([x, 0, z]), np.array([x, 0, z]) + force]
                l2 = [np.array([x + 0.1*force[2], 0, z + 0.1*force[0]]),
                      np.array([x + 0.1*force[2], 0, z + 0.1*force[0]]) + force]
                farr = FlatArrow(l1, l2, tip_size=1, tip_width=1).c(color='c', alpha=0.9)
                arrows.append(farr)

    arrows.append(plot_iterations_new(o).color("silver"))
    show(arrows, __doc__, viewup="z", axes=1, bg='bb', zoom=1, size=(1920, 1080)).close()

def reader_avoid_field_params_search():
    filename = 'storage/pid_const5_avoiding.txt'
    f = open(filename, 'r')
    k_p = []
    k_a = []
    res = []
    for line in f:
        lst = line.split()
        k_p.append(float(lst[0]))
        k_a.append(float(lst[1]))
        res.append(int(lst[2]))
    f.close()
    plt.title("Подбор коэффициентов ПД-регулятора, поля отталкивания")
    clr = ['lightblue', 'deeppink', 'palegreen']
    state_list = ['преследование', 'столкновение', 'попадание']
    flag = [False, False, False]
    for i in range(len(res)):
        state = int(res[i] + 1)
        if not flag[state]:
            plt.scatter(k_p[i], k_a[i], c=clr[state], label=state_list[state])
            flag[state] = True
        else:
            plt.scatter(k_p[i], k_a[i], c=clr[state])
        print(k_p[i], k_a[i], state)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel('коэффициент ПД-регулятора')
    plt.xlabel('коэффициент поля уклонения')
    plt.legend()
    plt.show()

def plot_avoid_field_params_search(dt=1., N=5, T_max=1000., k_p_min=1e-4, k_p_max=1e-2, k_a_min=1e-5, k_a_max=1e-3):
    """Фунция тестит коэффициенты ПД-регулятора и коэффициент отталкивания на конструкции 5.
    Результат - точки на пространстве {k_PD, k_avoid}, разделящиеся на классы:
    -> попадание в цель res=1
    -> столкновение res=0
    -> преследование res=-1"""
    filename = 'storage/pid_const5_avoiding.txt'
    f = open(filename, 'w')
    f.close()
    k_p_list = np.exp(np.linspace(np.log(k_p_min), np.log(k_p_max), N))  # Logarithmic scale
    k_av_list = np.exp(np.linspace(np.log(k_a_min), np.log(k_a_max), N))  # Logarithmic scale
    start_time = datetime.now()
    tmp_count = 0

    for k_p in k_p_list:
        for k_a in k_av_list:
            res = -1
            tmp_count += 1
            id_app = 0
            print(Fore.CYAN + f'Подбор ПД-к-в: {tmp_count}/{N**2}; время={datetime.now() - start_time}'
                  + Style.RESET_ALL)
            o = AllProblemObjects(if_PID_control=True,
                                  dt=dt, k_p=k_p, k_av=k_a,
                                  T_max=T_max,
                                  if_talk=False,
                                  if_any_print=False,
                                  choice='5')

            for i_time in range(int(T_max // dt)):
                # Repulsion
                if o.a.flag_fly[id_app] == 0:
                    _ = repulsion(o, id_app, u_a_priori=np.array([0.3, 0., 0.]))

                # Control
                o.time_step()
                o.control_step(0)

                # Docking
                res = 0 if call_crash(o, o.a.r[0], o.R, o.S, o.taken_beams) else res
                if not res:
                    break
                if o.get_discrepancy(id_app=0) < o.d_to_grab:
                    res = 1
                    break

            f = open(filename, 'a')
            f.write(f'{k_p} {k_a} {res}\n')
            print(f"{k_p} {k_a} {res}")
            f.close()
    reader_avoid_field_params_search()

def plot_repulsion_error(N: int = 3, dt: float = 5., T_max: float = 1000.):
    start_time = datetime.now()
    filename = 'storage/repulsion_error.txt'
    f = open(filename, 'w')
    k_u_list = [0.01, 0.02, 0.05, 0.1, 0.2]
    errors = [[[] for _ in range(4)] for _ in range(len(k_u_list))]
    clr = ['plum', 'skyblue', 'darkgreen', 'maroon']
    tmp = 0
    for choice in ['1', '2', '3', '4']:
        id_app = 0
        u = None
        for k in range(len(k_u_list)):
            for j in range(N):
                tmp += 1
                print(Fore.CYAN + f"Невязка от погрешности отталкивания: {tmp}/{N*len(k_u_list)*4}; "
                                  f"время={datetime.now() - start_time}")
                f_min = 1e5
                o = AllProblemObjects(dt=dt, T_max=T_max, if_any_print=False, choice=choice, method='shooting',
                                      d_crash=None, d_to_grab=None, if_talk=False)
                for i_time in range(int(T_max // dt)):
                    # Docking
                    f_min = min(f_min, o.get_discrepancy(id_app=0))

                    # Repulsion
                    if o.a.flag_fly[id_app] == 0:
                        u = repulsion(o, id_app) if u is None else \
                            repulsion(o, id_app, u_a_priori=velocity_spread(u, k_u_list[k]))

                    # Control
                    o.control_step(id_app)
                errors[k][int(choice) - 1].append(f_min)
                f.write(f"{choice} {k_u_list[k]} {f_min}")
                if j == 0:
                    plt.scatter(100 * k_u_list[k], f_min, c=clr[int(choice) - 1], label=choice)
                else:
                    plt.scatter(100 * k_u_list[k], f_min, c=clr[int(choice) - 1])

    plt.xlabel('относительный разброс скорости, %')
    plt.xlabel('погрешность, м')
    plt.legend()
    plt.show()
    f.close()

def plot_acceleration_error(dt=1., T_max=1000.):
    k_ac_list = [0.01, 0.02, 0.05, 0.1, 0.2]
