from all_objects import *


def def_o(dt=0.1):
    o1 = AllProblemObjects(if_talk=False)
    o1.w = np.array([1e-3, -3e-5, 0])
    o1.om_update()
    return o1


def test_full_energy(order, w=0.001, dt=1., T_max=1000.):
    # Сохранение полной энергии
    result = True
    e = 10**(-order)
    o = AllProblemObjects(if_talk=False, dt=dt)
    o.w = np.array([0, -w, 0])
    o.om_update()
    U0 = o.get_potential_energy()
    T = [o.get_kinetic_energy()]
    U = [o.get_potential_energy() - U0]
    E_0 = o.get_kinetic_energy() + o.get_potential_energy() - U0
    E_list = [E_0]
    t = [0.]
    for i in range(int(T_max/dt)):
        o.time_step()
        if i % 10 == 0:
            T.append(o.get_kinetic_energy())
            U.append(o.get_potential_energy() - U0)
            E = o.get_kinetic_energy() + o.get_potential_energy() - U0
            E_list.append(E)
            t.append(o.dt * i)
            if abs(E - E_0) / E_0 > e:
                result = False
            if abs(E - E_0) / E > e:
                result = False

    T1 = 2*np.pi/o.w_hkw
    for i in range(int(np.floor(T_max / T1))):
        plt.plot([T1 * (i + 1), T1 * (i + 1)], [np.min(U), np.max(T)], c='#5F9EA0')
    plt.title("График энергий")
    plt.plot(t, U, c='#CDC673', label='потенциальная')
    plt.plot(t, T, c='#CD6090', label='вращательная')
    plt.plot(t, E_list, c='#ADFF2F', label='полная')
    plt.legend()
    draw_reference_frames(o, size=10, showing=True)
    return result


def test_rotation(order, o=None, dt=1., w=0.001, T_max=1000.):
    # Правильность задания матриц вращения
    result = True
    e = 10**(-order)
    o = def_o() if o is None else o
    o.w = np.array([0, -w, 0])
    o.om_update()
    for i in range(int(T_max/dt)):
        o.time_step()
        if i % 10 == 0:
            for j in range(3):
                if o.U.T[0][j] > e:
                    if (o.U.T[0][j] - np.linalg.inv(o.U)[0][j]) / o.U.T[0][j] > e:
                        result = False
                if o.A.T[0][j] > e:
                    if (o.A.T[0][j] - np.linalg.inv(o.A)[0][j]) / o.A.T[0][j] > e:
                        result = False
                if np.linalg.inv(o.U)[0][j] > e:
                    if (o.U.T[0][j] - np.linalg.inv(o.U)[0][j]) / np.linalg.inv(o.U)[0][j] > e:
                        result = False
                if np.linalg.inv(o.S)[0][j] > e:
                    if (o.S.T[0][j] - np.linalg.inv(o.S)[0][j]) / np.linalg.inv(o.S)[0][j] > e:
                        result = False
    a = [0, 1, 2]
    c = np.array([0., 1., 2.])
    A = np.array([[0., 1., 2.], [3., 4., 5.], [6., 7., 8.]])
    for _ in range(20):
        Lu = np.random.rand(4)
        Lu /= np.linalg.norm(Lu)
        Ls = np.random.rand(4)
        Ls /= np.linalg.norm(Ls)
        U = quart2dcm(np.double(Lu))
        S = quart2dcm(np.double(Ls))
        if abs(1. - np.linalg.det(U)) > e:
            result = False
        if abs(1. - np.linalg.det(S)) > e:
            result = False
        b = o.i_o(o.o_i(o.i_o(o.o_i(a, U=U), U=U), U=U), U=U)
        for _ in range(50):
            b = o.i_o(o.o_i(o.i_o(o.o_i(b, U=U), U=U), U=U), U=U)
        for i in range(3):
            if abs(1. - np.linalg.norm([U[0][i], U[1][i], U[2][i]])) > e:
                result = False
            if abs(1. - np.linalg.norm([S[0][i], S[1][i], S[2][i]])) > e:
                result = False
            for j in range(3):
                if abs(U.T[j][i] - np.linalg.inv(U)[j][i]) > e:
                    result = False
                if abs(A[j][i] - o.b_i(o.i_b(A))[j][i]) > e:
                    result = False
            if abs(c[i] - o.i_o(o.o_i(c, U=U), U=U)[i]) > e:
                result = False
            if abs(a[i] - o.i_o(o.o_i(b, U=U), U=U)[i]) > e:
                result = False
            if abs(c[i] - o.o_b(o.b_o(c, S=S), S=S)[i]) > e:
                result = False
            if abs(a[i] - o.o_b(o.b_o(a, S=S), S=S)[i]) > e:
                result = False
            if abs(c[i] - o.i_b(o.b_i(c))[i]) > e:
                result = False
            if abs(a[i] - o.i_b(o.b_i(a))[i]) > e:
                result = False
            if abs(o.o_b(o.i_o(c))[i] - o.i_b(c)[i]) > e:
                result = False
            if abs(o.b_o(o.i_b(c))[i] - o.i_o(c)[i]) > e:
                result = False
    return result


def test_runge_kutta(order, o=None, dt=1., w=0.001, T_max=1000.):
    result = True
    e = 10 ** (-order)
    o = def_o() if o is None else o
    o.w = np.array([0, -w, 0])
    o.om_update()
    r = np.array([random.uniform(-10, 10), random.uniform(-10, 10), random.uniform(-10, 10)])
    v = np.array([random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1)])
    c = get_c_hkw(r, v, o.w_hkw)
    for i in range(int(T_max/dt)):
        r, v = o.rk4_acceleration(r, v, [0, 0, 0])
    r_h = r_hkw(c, o.w_hkw, T_max)
    v_h = v_hkw(c, o.w_hkw, T_max)
    for i in range(3):
        if abs(r[i] / r_h[i]) < e:
            result = False
        if abs(r_h[i] / r[i]) < e:
            result = False
        if abs(v[i] / v_h[i]) < e:
            result = False
        if abs(v_h[i] / v[i]) < e:
            result = False
    return result
