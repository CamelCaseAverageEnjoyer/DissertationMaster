"""Assembling general problem solution"""
from all_objects import *

vedo_picture = False
to_save = True
o_global = AllProblemObjects(if_impulse_control=False,
                             if_PID_control=False,
                             if_LQR_control=False,
                             if_avoiding=True,

                             is_saving=vedo_picture and to_save,
                             save_rate=50,
                             if_talk=False,
                             if_testing_mode=True,
                             choice_complete=False,

                             # method='shooting',
                             method='hkw_analytics+pd',  # [-1.06955519e-02 -9.05063216e-03 -5.05769658e-05]
                             # method='2d_analytics+pd',
                             # method='diffevolve+shooting+pd',
                             # method='diffevolve+trust-constr',
                             # method='shooting+imp',
                             # method='const-propulsion',
                             # method='linear-propulsion',
                             # method='linear-angle',
                             if_T_in_shooting=False,
                             begin_rotation='xx',

                             shooting_amount_repulsion=15,
                             diff_evolve_times=1,
                             diff_evolve_vectors=100,

                             dt=1.0, T_max=10000., u_max=0.05,
                             a_pid_max=1e-5, k_p=3e-4, freetime=50,
                             choice='3', floor=7, d_crash=0.2,
                             N_apparatus=1, file_reset=True)
'''for j in range(300):
    o_global.s.flag[j] = np.array([1, 1])'''
print(f"Количество стержней: {o_global.s.n_beams}")

def iteration_func(o):
    o.time_step()
    o.line_str_orf = np.append(o.line_str_orf, o.r_ub)

    for id_app in o.a.id:
        # Repulsion
        o.a.busy_time[id_app] -= o.dt if o.a.busy_time[id_app] >= 0 else 0
        if (not o.a.flag_fly[id_app]) and o.a.busy_time[id_app] < 0:  # [-0.0111501  -0.01204346 -0.00513348]
            print(f"отталкивание из файла {o.get_repulsion(id_app)}")
            # u = repulsion(o, id_app, u_a_priori=np.array([-0.0111501,  -0.01204346, -0.00513348]))
            u = repulsion(o, id_app, u_a_priori=o.get_repulsion(id_app))
            o.file_save(f'отталкивание {id_app} {u[0]} {u[1]} {u[2]}')
            o.repulsion_save(f'отталкивание {id_app} {u[0]} {u[1]} {u[2]}')

        # Motion control
        o.control_step(id_app)

        # Capturing
        discrepancy = o.get_discrepancy(id_app)
        if (discrepancy < o.d_to_grab) and o.a.flag_fly[id_app]:
            capturing(o=o, id_app=id_app)

        # Docking
        o.file_save(f'график {id_app} {discrepancy} {np.linalg.norm(o.w)} '
                    f'{np.linalg.norm(180 / np.pi * np.arccos(clip((np.trace(o.S) - 1) / 2, -1, 1)))} '
                    f'{np.linalg.norm(o.v_ub)} {np.linalg.norm(o.r_ub)} {np.linalg.norm(o.a_self[id_app])}')
        o.line_app_brf[id_app] = np.append(o.line_app_brf[id_app], o.o_b(o.a.r[id_app]))
        o.line_app_orf[id_app] = np.append(o.line_app_orf[id_app], o.a.r[id_app])

        # Stop criteria
        if np.linalg.norm(o.a.r[id_app]) > 1e3:
            o.my_print(f"МИССИЯ ПРОВАЛЕНА! ДОПУЩЕНА ПОТЕРЯ АППАРАТА!", mode='m')
            o.t = 2 * o.T_total
    return o

def iteration_timer(eventId=None):
    global o_global, vedo_picture, fig_view
    if o_global.t <= o_global.T_total:
        o_global = iteration_func(o_global)
        if vedo_picture and o_global.iter % o_global.save_rate == 0:
            fig_view = draw_vedo_and_save(o_global, o_global.iter, fig_view, app_diagram=False)

def button_func():
    global timerId
    fig_view.timer_callback("destroy", timerId)
    if "Play" in button.status():
        timerId = fig_view.timer_callback("create")
    button.switch()


if __name__ == "__main__":
    global timerId, fig_view, button, evnetId
    if vedo_picture:
        timerId = 1
        fig_view = Plotter(bg='white', size=(1920, 1080))
        button = fig_view.add_button(button_func, states=["Play ", "Pause"], size=20,
                                     font='Bongas', bold=True, pos=[0.9, 0.9])
        fig_view.timer_callback("destroy", timerId)
        evnetId = fig_view.add_callback("timer", iteration_timer)

        my_mesh = plot_iterations_new(o_global).color("silver")
        app_mesh = plot_apps_new(o_global)
        fig_view.show(__doc__, my_mesh + app_mesh, zoom=0.5)

    else:
        while True:
            iteration_timer()
