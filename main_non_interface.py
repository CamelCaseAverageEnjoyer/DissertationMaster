"""Assembling general problem solution"""
from all_objects import *

vedo_picture = False
o_global = AllProblemObjects(if_impulse_control=False,
                             if_PID_control=False,
                             if_LQR_control=False,
                             if_avoiding=False,

                             is_saving=False,
                             save_rate=10,
                             if_talk=True,
                             if_testing_mode=True,
                             choice_complete=False,

                             w_twist=0.,
                             w_max=1e5,
                             e_max=1e2,
                             j_max=1e5,
                             R_max=1e5,
                             # method='shooting',
                             # method='shooting+pd',
                             # method='shooting+imp',
                             method='diffevolve+const-propulsion',
                             begin_rotation='xx',
                             shooting_amount_repulsion=30,

                             diff_evolve_times=3,
                             diff_evolve_vectors=10,

                             dt=10.0, T_max=5000., u_max=0.03,
                             choice='3', floor=7, d_crash=0.2,
                             N_apparatus=1, file_reset=True)
'''for j in range(24):
    o_global.s.flag[j] = np.array([1, 1])'''
# o_global.a.r[0] = np.array([-10, 0., 0.])
print(f"Количество стержней: {o_global.s.n_beams}")

def iteration_func(o):
    o.time_step()
    o.line_str = np.append(o.line_str, o.R)

    for id_app in o.a.id:
        # Repulsion
        o.a.busy_time[id_app] -= o.dt if o.a.busy_time[id_app] >= 0 else 0
        if (not o.a.flag_fly[id_app]) and o.a.busy_time[id_app] < 0:
            # u = repulsion(o, id_app, u_a_priori=np.array([-0.01, 0., -0.01]))
            u = repulsion(o, id_app)
            o.file_save(f'отталкивание {id_app} {u[0]} {u[1]} {u[2]}')

        # Motion control
        o.control_step(id_app)

        # Capturing
        discrepancy = o.get_discrepancy(id_app)
        if (discrepancy < o.d_to_grab) and o.a.flag_fly[id_app]:
            capturing(o=o, id_app=id_app)

        # Docking
        o.file_save(f'график {id_app} {discrepancy} {o.get_e_deviation()} '
                    f'{np.linalg.norm(180 / np.pi * np.arccos(clip((np.trace(o.S) - 1) / 2, -1, 1)))} '
                    f'{np.linalg.norm(o.V)} {np.linalg.norm(o.R)} {np.linalg.norm(o.a_self[id_app])}')
        o.line_app[id_app] = np.append(o.line_app[id_app], o.o_b(o.a.r[id_app]))
        o.line_app_orf[id_app] = np.append(o.line_app_orf[id_app], o.a.r[id_app])

        # Stop criteria
        if np.linalg.norm(o.a.r[id_app]) > 1000:
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
        fig_view = Plotter(bg='bb', size=(1920, 1080))
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
