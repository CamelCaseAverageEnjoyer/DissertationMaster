from all_objects import *


o_global = AllProblemObjects(dt=10.)


'''def dr_gradient_my(u0, *args):
    o, T_max, id_app, interaction, mu_IPM, mu_e = args
    u0 = o.cases['repulse_vel_control'](u0)
    dr_, _, w_, V_, R_, j_, _ = calculation_motion(o=o, u=u0, T_max=T_max, id_app=id_app, interaction=interaction)
    reserve_rate = 1.5
    e_w = 0. if (o.w_max/reserve_rate - w_) > 0 else abs(w_ - o.w_max/reserve_rate)
    e_V = 0. if (o.V_max/reserve_rate - V_) > 0 else abs(V_ - o.V_max/reserve_rate)
    e_R = 0. if (o.R_max/reserve_rate - R_) > 0 else abs(R_ - o.R_max/reserve_rate)
    e_j = 0. if (o.j_max/reserve_rate - j_) > 0 else abs(j_ - o.j_max/reserve_rate)
    if o.if_testing_mode:
        print(Style.RESET_ALL + f"e: {e_w/o.w_max} | {e_V/o.V_max} | {e_R/o.R_max} | {e_j/o.j_max}")
    cnd = (e_w + e_V + e_R + e_j > 0)
    anw = dr_ - \
        mu_IPM * (np.log(o.w_max - w_ + e_w) + np.log(o.V_max - V_ + e_V) +
                  np.log(o.R_max - R_ + e_R) + np.log(o.j_max - j_ + e_j)) + \
        mu_e * (e_w/o.w_max + e_V/o.V_max + e_R/o.R_max + e_j/o.j_max)
    print(f"bbb {anw}")
    print(f"aaa {np.linalg.norm(anw)}")
    return np.linalg.norm(anw)


u0 = np.array(np.zeros(3))
res = scipy.optimize.minimize(dr_gradient_my, u0, args=(o, o.T_max, 0, True, o.mu_IPM, o.mu_e), tol=0.5)
print(res)'''

help(o_global.control_step)
