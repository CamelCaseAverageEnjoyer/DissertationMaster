from all_objects import *
from vedo import *
from datetime import datetime
k_ac_list = [0.01, 0.02, 0.05, 0.1, 0.2]


def reader_pd_control_params(name: str = '', eps: float = 1e-1, lng: str = 'ru'):
    global k_ac_list
    filename = 'storage/pid_koeff_' + name + '.txt'
    f = open(filename, 'r')
    k_p, k_a, tol, col = ([], [], [], [])
    box_k_p = [[] for _ in range(5)]
    for line in f:
        lst = line.split()
        k_p += [float(lst[0])]
        k_a += [int(lst[1])]
        tol += [float(lst[2])]
        col += [int(lst[3])]
    f.close()

    tol_line = [[] for _ in range(len(k_ac_list))]
    count = [0 for _ in range(len(k_ac_list))]
    k_p_max = 0.
    k_p_list = []
    for i in range(len(tol)):
        if k_p_max < k_p[i]:
            k_p_max = k_p[i]
            k_p_list += [k_p[i]]
            for k in range(len(k_ac_list)):
                if count[k_a[i]] > 0:
                    tol_line[k][len(k_p_list) - 2] /= count[k]
                tol_line[k].append(0.)
            count = [0 for _ in range(len(k_ac_list))]
        else:
            tol_line[k_a[i]][len(k_p_list) - 1] += tol[i]
            count[k_a[i]] += 1
    for k in range(len(k_ac_list)):
        tol_line[k][len(k_p_list) - 1] /= count[k]

    if lng == 'ru':
        title = "Подбор коэффициентов ПД-регулятора"
        x_label = "Коэффициент k_p"
        y_label = "Точность"
    else:
        title = "PD-controller coeffitient enumeration"
        x_label = "k_p coeffitient"
        y_label = "Tolerance"

    clr = [['aqua', 'violet', 'limegreen', 'gold', 'lightcoral'],
           ['steelblue', 'purple', 'green', 'goldenrod', 'firebrick']]
    plt.title(title)
    for i in range(len(k_ac_list)):
        plt.plot(k_p_list, tol_line[i], c=clr[0][i], label=k_ac_list[i])
        for kp in k_p_list:
            plt.plot([kp * (1 + eps * (i - 2))] * 2, [0., 20.], c=clr[0][i])
    for i in range(len(tol)):
        plt.scatter(k_p[i] * (1 + eps * (k_a[i] - 2)), tol[i], c=clr[col[i]][k_a[i]])
    plt.xscale("log")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()

def pd_control_params_search(name: str = '', dt=0.2, n_p=5, n_a=10, T_max=700., k_min=1e-4, k_max=1e-2):
    """Функция ищет смысл жизни(его нет) / зависимость невязки/столкновения от """
    global k_ac_list
    k_p_list = np.exp(np.linspace(np.log(k_min), np.log(k_max), n_p))  # Logarithmic scale
    k_p_best = 0
    tolerance_best = 1e5

    filename = 'storage/pid_koeff_' + name + '.txt'
    f = open(filename, 'w')
    start_time = datetime.now()
    tmp_count = 0
    collide = False

    for k_p in k_p_list:
        tol_count = 0.
        for k_a_i in range(len(k_ac_list)):
            for _ in range(n_a):
                tmp_count += 1
                id_app = 0
                tolerance = None
                print(Fore.CYAN + f'Подбор ПД-к-в: {tmp_count}/{n_p * n_a * len(k_ac_list)};'
                                  f' время={datetime.now() - start_time}' + Style.RESET_ALL)
                o = AllProblemObjects(if_PID_control=True,
                                      dt=dt, k_p=k_p, k_ac=k_ac_list[k_a_i],
                                      T_max=T_max,
                                      if_talk=False,
                                      if_any_print=False,
                                      choice='3')

                for i_time in range(int(T_max // dt)):
                    # Repulsion
                    o.a.busy_time[id_app] -= o.dt if o.a.busy_time[id_app] >= 0 else 0
                    if o.a.flag_fly[id_app] == 0 and o.a.busy_time[id_app] < 0:
                        _ = repulsion(o, id_app, u_a_priori=np.array([-0.00749797, 0.00605292, 0.08625441]))

                    o.time_step()
                    o.control_step(id_app)

                    tmp = o.get_discrepancy(id_app)
                    collide = call_crash(o, o.a.r[id_app], o.R, o.S, o.taken_beams)
                    if tolerance is None or tolerance > tmp:
                        tolerance = tmp
                    if collide:
                        break
                tol_count += tolerance
                f.write(f'{k_p} {k_a_i} {tolerance} {int(collide)}\n')
        tol_count /= n_a * len(k_ac_list)
        if tol_count < tolerance_best:
            tolerance_best = tol_count
            k_p_best = k_p
    f.close()
    reader_pd_control_params(name=name)
    return k_p_best

def plot_params_while_main(filename: str = "", trial_episodes: bool = False, show_rate: int = 1, limit: int = 1e5,
                           show_probe_episodes=True):
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
    R_max, V_max, j_max, w_max = (1., 1., 1., 1.)

    f = open('storage/main.txt', 'r')
    tmp = 0
    for line in f:
        lst = line.split()
        if len(dr[0]) < limit:
            if len(lst) > 0 and tmp % show_rate == 0:
                if lst[0] == 'ограничения' and len(lst) == 5:
                    print(f"Есть ограничения")
                    R_max = float(lst[1])
                    V_max = float(lst[2])
                    j_max = float(lst[3])
                    w_max = float(lst[4])
                if lst[0] == 'график' and len(lst) == 9 and (show_probe_episodes or bool(int(lst[8]))):
                    id_app = int(lst[1])
                    dr[id_app].append(float(lst[2]))
                    w[id_app].append(float(lst[3]))
                    j[id_app].append(float(lst[4]))
                    V[id_app].append(float(lst[5]))
                    R[id_app].append(float(lst[6]))
                    a[id_app].append(float(lst[7]))
                    m[id_app].append(int(lst[8]))
        else:
            print(Fore.MAGENTA + f"Внимание! Превышен лимит в {limit} точек!" + Style.RESET_ALL)
            break
        tmp += 1
    f.close()
    print(Fore.CYAN + f"Аппаратов на графиках: {id_max}" + Style.RESET_ALL)
    print(Fore.BLUE + f"Точек на графиах: {len(dr[0])}" + Style.RESET_ALL)

    fig, axs = plt.subplots(3)
    axs[0].set_xlabel('время, с')
    axs[0].set_ylabel('Невязка, м')
    axs[0].set_title('Параметры в процессе алгоритма')
    axs[1].set_xlabel('время t, с')
    axs[1].set_ylabel('ограниченные величины')
    axs[2].set_xlabel('итерации')
    axs[2].set_ylabel('бортовое ускорение, м/с2')

    clr = ['c', 'indigo', 'm', 'violet', 'teal', 'slategray', 'greenyellow', 'sienna']
    clr2 = [['skyblue', 'bisque', 'palegreen', 'darksalmon'], ['teal', 'tan', 'g', 'brown']]
    if show_probe_episodes:
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
        axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(w[id_app][0]) / w_max, np.array(w[id_app][1]) /
                                                   w_max], c=clr[1][0], label='w')
        axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(j[id_app][0]) / j_max, np.array(j[id_app][1]) /
                                                   j_max], c=clr[1][1], label='угол')
        axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(V[id_app][0]) / V_max, np.array(V[id_app][1]) /
                                                   V_max], c=clr[1][2], label='V')
        axs[1].plot([t[id_app][0], t[id_app][1]], [np.array(R[id_app][0]) / R_max, np.array(R[id_app][1]) /
                                                   R_max], c=clr[1][3], label='R')
        for i in range(len(t[id_app]) - 1):
            axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(w[id_app][i]) / w_max, np.array(w[id_app][i+1]) /
                                                         w_max], c=clr[m[id_app][i]][0])
            axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(j[id_app][i]) / j_max, np.array(j[id_app][i+1]) /
                                                         j_max], c=clr[m[id_app][i]][1])
            axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(V[id_app][i]) / V_max, np.array(V[id_app][i+1]) /
                                                         V_max], c=clr[m[id_app][i]][2])
            axs[1].plot([t[id_app][i], t[id_app][i+1]], [np.array(R[id_app][i]) / R_max, np.array(R[id_app][i+1]) /
                                                         R_max], c=clr[m[id_app][i]][3])
    else:
        for id_app in range(id_max):
            t[id_app] = np.linspace(0, len(dr[id_app]), len(dr[id_app]))
            axs[0].plot(t[id_app], dr[id_app], c=clr[2 * id_app])
            axs[1].plot(t[id_app], [1 for _ in range(len(t[id_app]))], c='gray')
            axs[2].plot(range(len(a[id_app])), a[id_app], c='c')
            axs[0].plot(t[id_app], np.zeros(len(t[id_app])), c='khaki')
            axs[2].plot(range(len(a[id_app])), np.zeros(len(a[id_app])), c='khaki')
        id_app = 0
        axs[1].plot(t[id_app], np.array(w[id_app]) / w_max, c=clr2[1][0], label='w')
        axs[1].plot(t[id_app], np.array(j[id_app]) / j_max, c=clr2[1][1], label='угол')
        axs[1].plot(t[id_app], np.array(V[id_app]) / V_max, c=clr2[1][2], label='V')
        axs[1].plot(t[id_app], np.array(R[id_app]) / R_max, c=clr2[1][3], label='R')
    axs[1].legend()
    plt.show()


def plot_a_avoid(x_boards: list = np.array([-10, 10]), z_boards: list = np.array([-5, 10])):
    o = AllProblemObjects(choice='0', N_apparatus=0)

    # x_boards: list = [-15, 15], z_boards: list = [1, 10]
    nx = 60
    nz = 60
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
                l2 = [np.array([x - 0.1*force[2], 0, z + 0.1*force[0]]),
                      np.array([x - 0.1*force[2], 0, z + 0.1*force[0]]) + force]
                farr = FlatArrow(l1, l2, tip_size=1, tip_width=1).c(color='c', alpha=0.9)
                arrows.append(farr)

    arrows.append(plot_iterations_new(o).color("silver"))
    show(arrows, __doc__, viewup="z", axes=1, bg='white', zoom=1, size=(1920, 1080)).close()

def reader_avoid_field_params_search(filename: str = '', lng: str = 'ru'):
    # Init
    f = open(filename, 'r')
    k_p, k_a, lvl, res = ([], [], [], [])
    for line in f:
        lst = line.split()
        k_p += [float(lst[0])]
        k_a += [float(lst[1])]
        lvl += [int(lst[2])]
        res += [int(lst[3])]
    f.close()

    # Titles
    if lng == 'ru':
        title = "Подбор коэффициентов ПД-регулятора, поля отталкивания"
        state_list = ['преследование', 'столкновение', 'попадание']
        xlabel = "коэффициент ПД-регулятора"
        ylabel = "коэффициент поля уклонения"
    else:
        title = "?"
        state_list = ['', 'collision', 'reaching']
        xlabel = "PD-controller coeffitient"
        ylabel = "avoiding field coeffitient"

    # Plotting
    clr = ['lightblue', 'deeppink', 'palegreen']
    flag = [True] * 3
    for i in range(len(res)):
        state = int(res[i] + 1)
        if flag[state]:
            plt.scatter(k_p[i], k_a[i], c=clr[state], label=state_list[state])
            flag[state] = False
        else:
            plt.scatter(k_p[i], k_a[i], c=clr[state])
    plt.title(title)
    plt.yscale("log")
    plt.xscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.show()

def plot_avoid_field_params_search(name: str = '', dt=1., N=20, T_max=1000., k_p_min=1e-4, k_p_max=1e-3,
                                   k_a_min=1e-8, k_a_max=1e-1):
    """Фунция тестит коэффициенты ПД-регулятора и коэффициент отталкивания на конструкции 5.
    Результат - точки на пространстве {k_PD, k_avoid}, разделящиеся на классы:
    -> попадание в цель res=1
    -> столкновение res=0
    -> преследование res=-1"""
    k_p_list = np.exp(np.linspace(np.log(k_p_min), np.log(k_p_max), N))  # Logarithmic scale
    k_av_list = np.exp(np.linspace(np.log(k_a_min), np.log(k_a_max), N))  # Logarithmic scale
    start_time = datetime.now()
    tmp_count = 0
    for lvl in [222]:
        filename = 'storage/pid_const5_avoiding_' + name + '_' + str(lvl) + '.txt'
        # f = open(filename, 'a')
        f = open(filename, 'w')

        for k_p in k_p_list:
            for k_a in k_av_list:
                res = -1
                tmp_count += 1
                id_app = 0
                print(Fore.CYAN + f'Подбор ПД-к-в: {tmp_count}/{N**2 * 4}; время={datetime.now() - start_time}'
                    + Style.RESET_ALL)
                o = AllProblemObjects(if_PID_control=True, if_avoiding=True,
                                    dt=dt, k_p=k_p, k_av=k_a,
                                    T_max=T_max, level_avoid=lvl,
                                    if_talk=False,
                                    if_any_print=False,
                                    choice='5')

                for _ in range(int(T_max // dt)):
                    # Repulsion
                    if o.a.flag_fly[id_app] == 0:
                        _ = repulsion(o, id_app, u_a_priori=np.array([-random.uniform(0.01, 0.05), 0., 0.]))

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

                f.write(f'{k_p} {k_a} {lvl} {res}\n')
        # reader_avoid_field_params_search(name=name)
        f.close()

def reader_repulsion_error(name: str = ''):
    k_u_list = [0.01, 0.02, 0.05, 0.1, 0.2]
    k_w_list = [1e-5, 1e-4, 3e-4, 6e-4, 1e-3]
    filename = 'storage/repulsion_error_' + name + '.txt'
    f = open(filename, 'r')
    w, k, tol = ([], [], [])
    for line in f:
        lst = line.split()
        if len(lst) > 1:
            w += [k_w_list[int(lst[0])]]
            k += [int(lst[1])]
            tol += [float(lst[2])]

    tol_line = [[] for _ in range(len(k_u_list))]
    count = [0 for _ in range(len(k_u_list))]
    w_max = 0.
    i_max = 0

    '''for i in range(len(tol)):
        if w_max < w[i]:
            w_max = w[i]
            for k in range(len(k_u_list)):
                if count[k[i]] > 0:
                    tol_line[k][i_max - 2] /= count[k]
                tol_line[k].append(0.)
            count = [0 for _ in range(len(k_u_list))]
        else:
            tol_line[k_a[i]][len(k_p_list) - 1] += tol[i]
            count[k_a[i]] += 1'''

    plt.scatter(np.array(w) + (np.array(k) - 2) * 1e-5, tol, c='c')
    plt.show()
    
    '''for line in f:
        lst = line.split()
        k_p += [float(lst[0])]
        k_a += [int(lst[1])]
        tol += [float(lst[2])]
        col += [int(lst[3])]
    f.close()

    tol_line = [[] for _ in range(len(k_ac_list))]
    count = [0 for _ in range(len(k_ac_list))]
    k_p_max = 0.
    k_p_list = []
    for i in range(len(tol)):
        if k_p_max < k_p[i]:
            k_p_max = k_p[i]
            k_p_list += [k_p[i]]
            for k in range(len(k_ac_list)):
                if count[k_a[i]] > 0:
                    tol_line[k][len(k_p_list) - 2] /= count[k]
                tol_line[k].append(0.)
            count = [0 for _ in range(len(k_ac_list))]
        else:
            tol_line[k_a[i]][len(k_p_list) - 1] += tol[i]
            count[k_a[i]] += 1
    for k in range(len(k_ac_list)):
        tol_line[k][len(k_p_list) - 1] /= count[k]

    if lng == 'ru':
        title = "Подбор коэффициентов ПД-регулятора"
        x_label = "Коэффициент k_p"
        y_label = "Точность"
    else:
        title = "PD-controller coeffitient enumeration"
        x_label = "k_p coeffitient"
        y_label = "Tolerance"

    clr = [['aqua', 'violet', 'limegreen', 'gold', 'lightcoral'],
           ['steelblue', 'purple', 'green', 'goldenrod', 'firebrick']]
    plt.title(title)
    for i in range(len(k_ac_list)):
        plt.plot(k_p_list, tol_line[i], c=clr[0][i], label=k_ac_list[i])
        for kp in k_p_list:
            plt.plot([kp * (1 + eps * (i - 2))] * 2, [0., 20.], c=clr[0][i])
    for i in range(len(tol)):
        plt.scatter(k_p[i] * (1 + eps * (k_a[i] - 2)), tol[i], c=clr[col[i]][k_a[i]])
    plt.xscale("log")
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.show()'''


def plot_repulsion_error(name: str = '', N: int = 40, dt: float = 1., T_max: float = 1000.):
    k_u_list = [0.01, 0.02, 0.05, 0.1, 0.2]
    k_w_list = [1e-5, 1e-4, 3e-4, 6e-4, 1e-3]
    start_time = datetime.now()
    filename = 'storage/repulsion_error_' + name + '.txt'
    f = open(filename, 'w')
    # errors = [[[] for _ in range(4)] for _ in range(len(k_u_list))]
    clr = ['plum', 'skyblue', 'darkgreen', 'maroon']
    tmp = 0
    id_app = 0
    for k in range(len(k_u_list)):
        for w in range(len(k_w_list)):
            u = None
            for j in range(N):
                tmp += 1
                print(Fore.CYAN + f"Невязка от погрешности отталкивания: {tmp}/{N*len(k_u_list)*len(k_w_list)}; "
                                f"время={datetime.now() - start_time}")
                f_min = 1e5
                o = AllProblemObjects(dt=dt, T_max=T_max, if_any_print=False, choice='4', method='shooting',
                                    d_crash=None, d_to_grab=None, if_talk=False, floor=7)
                for j in range(184):
                    o.s.flag[j] = np.array([1, 1])
                o.w = np.array([0., k_w_list[w], 0.])
                o.om_update()
                for i_time in range(int(T_max // dt)):
                    # Docking
                    f_min = min(f_min, o.get_discrepancy(id_app=0))

                    # Repulsion
                    if o.a.flag_fly[id_app] == 0:
                        u = repulsion(o, id_app) if u is None else \
                            repulsion(o, id_app, u_a_priori=velocity_spread(u, k_u_list[k]))

                    # Control
                    o.time_step()
                # errors[k][int(choice) - 1].append(f_min)
                f.write(f"{w} {k} {f_min}\n")
    f.close()


'''current_beam = 0
current_percentage = 0.
go_forw = True
u_crawl = 0.4
o = AllProblemObjects(choice='2', floor=10)

def get_position_along_beam(o, persentage):
    global current_beam
    persentage = clip(persentage, 0, 1)
    return o.b_o(o.s.r1[current_beam] * (1 - persentage) + o.s.r2[current_beam] * persentage)

def local_iteration_func():
    global current_beam, u_crawl, o, current_percentage, go_forw
    o.iter += 1
    o.t = o.iter * o.dt
    id_app = 0
    o.line_str = np.append(o.line_str, o.R)
    o.line_app[id_app] = np.append(o.line_app[id_app], o.o_b(o.a.r[id_app]))
    o.line_app_orf[id_app] = np.append(o.line_app_orf[id_app], o.a.r[id_app])

    # Euler rotation
    o.La, o.Om = o.rk4_w(o.La, o.Om, o.J, o.t)
    o.U, o.S, o.A, o.R_e = o.call_rotation_matrix()
    o.w_update()

    if current_percentage == 1 and go_forw or current_percentage == 0 and not go_forw:
        tmp = o.get_discrepancy(id_app=0)
        if tmp < 1e-1:
            go_forw = not go_forw
            current_percentage = 0.
            o.a.flag_beam[0] = None
            o.s.flag[id_beam] = np.array([1., 1.])
        else:
            id_node = o.s.id_node[current_beam][0] if np.linalg.norm(o.a.r[id_app] - o.b_o(o.s.r1[current_beam])) < 1e-1 else o.s.id_node[current_beam][1]
            for j in range(o.s.n_beams):
                if id_node in o.s.id_node[j]:
                    tmp1 = o.get_discrepancy(id_app=0, r=o.s.r1[j])
                    tmp2 = o.get_discrepancy(id_app=0, r=o.s.r2[j])
                    if tmp > min(tmp1, tmp2):
                        tmp = min(tmp1, tmp2)
                        current_beam = j
                        current_percentage = 0.

    if current_percentage == 1 and not go_forw or current_percentage == 0 and go_forw:
        if np.linalg.norm(o.a.r[id_app] - o.b_o(o.s.r1[0])) < 1e-1:
            go_forw = True
            id_beam = o.s.call_possible_transport(o.taken_beams)[0]
            o.a.flag_beam[0] = id_beam
        id_beam = o.a.flag_beam[0]
        if o.a.flag_beam[0] is None:
            id_beam = 0
        o.a.target[0] = o.s.r1[id_beam]
    percentage_step = o.dt * u_crawl / o.s.length[current_beam]

    r_step_back = get_position_along_beam(o, current_percentage - percentage_step)
    r_step_forw = get_position_along_beam(o, current_percentage + percentage_step)
    d_back = o.get_discrepancy(r=r_step_back, id_app=0)
    d_forw = o.get_discrepancy(r=r_step_forw, id_app=0)
    if d_back < d_forw:
        o.a.r[0] = r_step_back 
        current_percentage = clip(current_percentage - percentage_step, 0, 1)
    else:
        o.a.r[0] = r_step_forw
        current_percentage = clip(current_percentage + percentage_step, 0, 1)
    print(f"b:{current_beam}, p:{current_percentage}, id:{o.a.flag_beam[0]}")

def iteration_timer(eventId=None):
    global o, vedo_picture, fig_view
    if o.t <= o.T_total:
        local_iteration_func()
    fig_view = draw_vedo_and_save(o, o.iter, fig_view, app_diagram=False)

def button_func():
    global timerId
    fig_view.timer_callback("destroy", timerId)
    if "Play" in button.status():
        timerId = fig_view.timer_callback("create")
    button.switch()'''

def crawling(u_crawl: float = 0.01):
    '''global timerId, fig_view, button, evnetId, o
    for j in range(30):
        o.s.flag[j] = np.array([1, 1])
    o.a.target[0] = o.b_o(np.zeros(3))
    o.a.r[0] = o.b_o(np.zeros(3))
    o.warning_message = False

    timerId = 1
    fig_view = Plotter(bg='bb', size=(1920, 1080))
    button = fig_view.add_button(button_func, states=["Play ", "Pause"], size=20,
                                    font='Bongas', bold=True, pos=[0.9, 0.9])
    fig_view.timer_callback("destroy", timerId)
    evnetId = fig_view.add_callback("timer", iteration_timer)

    my_mesh = plot_iterations_new(o).color("silver")
    app_mesh = plot_apps_new(o)
    fig_view.show(__doc__, my_mesh + app_mesh, zoom=0.5)'''
    o = AllProblemObjects(choice='2', floor=10)
    total_length = 0.
    for i in range(o.s.n_beams):
        i_floor = int((o.s.r1[i][0] - o.s.r1[6][0]) // o.s.length[6]) + 1
        # print(f"i:{i} floor:{i_floor}")
        if i_floor:
            total_length += 2 * (o.s.length[0] + o.s.length[6] + o.s.container_length)
        else:
            total_length += 2 * (o.s.length[0] + o.s.container_length)
        print(f"{total_length / u_crawl} секунд: установлен стержень id:{i}")
    print(f"Время сборки: {total_length / u_crawl} секунд")
