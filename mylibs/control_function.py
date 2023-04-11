from mylibs.calculation_functions import *


def avoiding_force(o, id_app, r=None):
    if r is None:
        r = o.o_b(o.a.r[id_app])
    else:
        r = o.o_b(r)
    force = np.zeros(3)

    for i in range(o.N_beams):
        if not(np.any(o.taken_beams == i)):
            if np.sum(o.s.flag[i]) > 0:
                r1 = o.s.r1[i]
                r2 = o.s.r2[i]
            else:
                r1 = [o.s.r_st[i][0], o.s.r_st[i][1], o.s.r_st[i][2]]
                r2 = [o.s.r_st[i][0] - o.s.length[i], o.s.r_st[i][1], o.s.r_st[i][2]]
            tmp = call_crash_internal_func(r, r1, r2, o.d_crash, return_force=True, k_av=o.k_av, level=o.level_avoid)
            if tmp is not False:
                if o.N_app > 0:
                    if np.linalg.norm(tmp) > np.linalg.norm(force) and np.linalg.norm(o.a.target[id_app] - r1) > 0.5:
                        force = tmp.copy()
                else:
                    if np.linalg.norm(tmp) > np.linalg.norm(force):
                        force = tmp.copy()

    for i in range(o.N_cont_beams):
        r1 = o.c.r1[i]
        r2 = o.c.r2[i]
        tmp_point = (np.array(r1) + np.array(r2)) / 2
        tmp_force = call_crash_internal_func(r, r1, r2, o.c.diam[i], return_force=True, k_av=o.k_av, level=o.level_avoid)
        if tmp_force is not False:
            if o.N_app > 0:
                if np.linalg.norm(tmp_force) > np.linalg.norm(force) and np.linalg.norm(o.a.target[id_app] - tmp_point) > 0.5:
                    force = tmp_force.copy()
            else:
                if np.linalg.norm(tmp_force) > np.linalg.norm(force):
                    force = tmp_force.copy()
    return o.S.T @ force

def pd_control(o, id_app):
    o.a.flag_hkw[id_app] = False
    r = np.array(o.a.r[id_app])
    v = np.array(o.a.v[id_app])
    dr = o.get_discrepancy(id_app, vector=True)
    dv = (dr - o.dr_p[id_app]) / o.dt
    o.dr_p[id_app] = dr.copy()
    r1 = o.a.target[id_app] - o.r_center
    a_pid = -o.k_p * dr - o.k_d * dv  # \
    '''         + o.S.T @ (my_cross(o.S @ o.e, r1) + my_cross(o.S @ o.w, my_cross(o.S @ o.w, r1)) +
                        2 * my_cross(o.S @ o.w, o.S @ o.a.v[id_app])) \
             + o.A_orbital - o.a_orbital[id_app]'''
    o.a_self[id_app] = a_pid.copy()

def lqr_control(o, id_app):
    o.a.flag_hkw[id_app] = False
    r = np.array(o.a.r[id_app])
    v = np.array(o.a.v[id_app])
    rv = np.append(r, v)
    r1 = o.a.target[id_app] - o.r_center
    muRe = o.mu / o.Radius_orbit ** 3
    tmp = 3 / o.Radius_orbit ** 4
    '''a = np.array([[0, 0, 0, 1., 0, 0],
                  [0, 0, 0, 0, 1., 0],
                  [0, 0, 0, 0, 0, 1.],
                  [(-muRe*o.A[0][0] + o.Om[1]**2 + o.Om[2]**2 + tmp*o.R_e[0]*o.S[0][2]),
                   (-muRe*o.A[1][0] + o.Om[2] - o.Om[0]*o.Om[1] + tmp*o.R_e[0]*o.S[1][2]),
                   (-muRe*o.A[2][0] - o.Om[1] - o.Om[0]*o.Om[2] + tmp*o.R_e[0]*o.S[2][2]),
                   0, 2*o.Om[2], -2*o.Om[1]],
                  [(-muRe*o.A[0][1] - o.Om[2] - o.Om[1]*o.Om[0] + tmp*o.R_e[1]*o.S[0][2]),
                   (-muRe*o.A[1][1] + o.Om[0]**2 + o.Om[2]**2 + tmp*o.R_e[1]*o.S[1][2]),
                   (-muRe*o.A[2][1] + o.Om[0] - o.Om[1]*o.Om[2] + tmp*o.R_e[1]*o.S[2][2]),
                   -2*o.Om[2], 0, 2*o.Om[0]],
                  [(-muRe*o.A[0][2] + o.Om[1] - o.Om[2]*o.Om[0] + tmp*o.R_e[2]*o.S[0][2]),
                   (-muRe*o.A[1][2] - o.Om[0] - o.Om[2]*o.Om[1] + tmp*o.R_e[2]*o.S[1][2]),
                   (-muRe*o.A[2][2] + o.Om[0]**2 + o.Om[1]**2 + tmp*o.R_e[2]*o.S[2][2]),
                   2*o.Om[1], -2*o.Om[0], 0]])'''
    a = np.array([[0, 0, 0, 1., 0, 0],
                  [0, 0, 0, 0, 1., 0],
                  [0, 0, 0, 0, 0, 1.],
                  [o.Om[1]**2 + o.Om[2]**2, o.e[2] - o.Om[0]*o.Om[1], -o.e[1] - o.Om[0]*o.Om[2],
                   0, 2*(o.Om[2] - o.W_hkw[2]), -2*o.w_hkw - 2*(o.Om[1] - o.W_hkw[1])],
                  [-o.e[2] - o.Om[0]*o.Om[1], o.Om[0]**2 + o.Om[2]**2 - o.w_hkw**2, o.e[0] - o.Om[1]*o.Om[2],
                   -2*(o.Om[2] - o.W_hkw[2]), 0, 2*(o.Om[0] - o.W_hkw[0])],
                  [o.e[1] - o.Om[0]*o.Om[2], -o.e[0] - o.Om[1]*o.Om[2], o.Om[0]**2 + o.Om[1]**2 + 3*o.w_hkw**2,
                   2*o.w_hkw + 2*(o.Om[1] - o.W_hkw[1]), -2*(o.Om[0] - o.W_hkw[0]), 0]])
    # print(f"a:{a}")
    b = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0],
                  [1., 0, 0],
                  [0, 1., 0],
                  [0, 0, 1.]])
    x_rate = 1
    u_rate = 1e9
    q = np.eye(6) * x_rate
    r = np.eye(3) * u_rate
    p = scipy.linalg.solve_continuous_are(a, b, q, r)

    # print(f"k:{np.linalg.inv(r) @ b.T @ p}")
    a_lqr = - np.linalg.inv(r) @ b.T @ p @ rv
    # a_lqr += tmp * o.R_e * (o.R[2] - (o.S.T @ o.a.target[id_app])[2]) - tmp * o.R_e * o.R[2] - tmp * o.U.T @ o.R * o.R[2] + \
    #          tmp * o.U.T @ o.a.r[id_app] + muRe * o.A.T @ o.a.target[id_app]
    a_lqr += my_cross(o.W_hkw, my_cross(o.W_hkw, o.R)) + \
             o.orbital_acceleration(np.append(o.S.T @ (o.a.target[id_app] - o.r_center), o.V + my_cross(o.w, o.a.target[id_app] - o.r_center)))
    print(f"проверка: {np.linalg.norm(a_lqr) / o.a_pid_max * 100}%")
    '''a_lqr = - np.linalg.inv(r) @ b.T @ p @ rv \
            + o.S.T @ (my_cross(o.S @ o.e, r1) + my_cross(o.S @ o.w, my_cross(o.S @ o.w, r1)) +
                       2 * my_cross(o.S @ o.w, o.S @ o.a.v[id_app])) \
            + o.A_orbital - o.a_orbital[id_app]'''
    # a_lqr -= o.a_self[id_app]
    a_lqr *= clip(np.linalg.norm(a_lqr), 0, o.a_pid_max) / np.linalg.norm(a_lqr)
    o.a_self[id_app] = a_lqr.copy()

def impulse_control(o, id_app):
    if not o.flag_vision[id_app]:
        o.t_flyby_counter -= o.dt
        o.t_reaction_counter = o.t_reaction
        if o.t_flyby_counter < 0:  # Облёт конструкции
            u = find_repulsion_velocity(o=o, id_app=id_app, interaction=False)
            o.t_flyby_counter = o.t_flyby
            o.t_start[id_app] = o.t
            talk_flyby(o.if_talk)
            o.C_r[id_app] = get_c_hkw(o.a.r[id_app], u, o.w_hkw)
    else:
        o.t_flyby_counter = o.t_flyby
        o.t_reaction_counter -= o.dt
        if o.t_reaction_counter < 0:  # Точное попадание в цель
            r_right = o.b_o(o.a.target[id_app])
            u, target_is_reached = calc_shooting(o=o, id_app=id_app, r_right=r_right, interaction=False)
            o.t_reaction_counter = o.t_reaction
            o.t_start[id_app] = o.t
            talk_shoot(o.if_talk)
            o.flag_impulse = not target_is_reached
            o.C_r[id_app] = get_c_hkw(o.a.r[id_app], u, o.w_hkw)

def control_condition(o, id_app):
    target_orf = o.b_o(o.a.target[id_app])
    see_rate = 1
    not_see_rate = 5
    if o.flag_vision[id_app]:  # If app have seen target, then app see it due to episode end
        o.t_reaction_counter -= o.dt
        return True
    if (o.a.flag_fly[id_app] == 1) and ((o.flag_vision[id_app] and ((o.iter % see_rate) == 0)) or
                                        ((not o.flag_vision[id_app]) and ((o.iter % not_see_rate) == 0))):
        if ((o.if_impulse_control and o.flag_impulse) or o.if_PID_control) and o.a.flag_fly[id_app]:
            points = 40
            o.flag_vision[id_app] = True
            for j in range(points):
                intermediate = (target_orf * j + np.array(o.a.r[id_app]) * (points - j)) / points
                o.flag_vision[id_app] = False if call_crash(o, intermediate, o.R, o.S, o.taken_beams) else o.flag_vision[id_app]
            o.t_reaction_counter = o.t_reaction_counter - o.dt*see_rate if o.flag_vision[id_app] else o.t_reaction
            return o.flag_vision[id_app]
    else:
        return True
    return False
