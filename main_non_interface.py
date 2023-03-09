"""General problem solution"""
from all_objects import *

vedo_picture = True
o = AllProblemObjects(if_impulse_control=False,
                      if_PID_control=False,
                      if_LQR_control=False,
                      if_avoiding=False,
                      is_saving=False,
                      save_rate=1,
                      if_talk=False,
                      if_multiprocessing=True,
                      if_testing_mode=True,
                      dt=5.,
                      choice='1',
                      T_max=1000.,
                      u_max=0.1,
                      N_apparatus=1)


def time_is(t, t0):
    print(f'|TIME: model iterative {t} | real program work {datetime.now() - t0}')

def iteration_func(o, f):
    o.time_step()
    o.line_str = np.append(o.line_str, o.R)  # Center mass line

    for id_app in o.X_app.id:
        # Repulsion
        o.X_app.loc[id_app, 'busy_time'] -= o.dt if o.X_app.busy_time[id_app] >= 0 else 0
        if (not o.X_app.flag_fly[id_app]) and o.X_app.busy_time[id_app] < 0:
            # u = repulsion(o, id_app, u_a_priori=np.array([-0.04, 0., 0.00]))
            u = repulsion(o, id_app)
            o.file_save(f, f'отталкивание {id_app} {u[0]} {u[1]} {u[2]}\n')

        # Motion control
        o.control_step(id_app)

        # Capturing
        discrepancy = o.get_discrepancy(id_app)
        # if (discrepancy < o.d_to_grab) and o.X_app.flag_fly[id_app]:
            # capturing(o=o, id_app=id_app)

        # Docking
        o.file_save(f, f'график {id_app} {discrepancy} {np.linalg.norm(o.w)} '
                       f'{np.linalg.norm(180 / np.pi * np.arccos(clip((np.trace(o.S) - 1) / 2, -1, 1)))} '
                       f'{np.linalg.norm(o.V)} {np.linalg.norm(o.R)}\n')
        o.file_save(f, f'управление {id_app} {np.linalg.norm(o.a_self[id_app])}\n')
        # if o.X_app.flag_fly[id_app]:
        o.line_app[id_app] = np.append(o.line_app[id_app], o.o_b(o.X_app.r[id_app]))  # Line of flying app
        o.line_app_orf[id_app] = np.append(o.line_app_orf[id_app], o.X_app.r[id_app])
    return o, f

def iteration_timer(eventId=None):
    global timerId, button, fig_view, o, f, start_time, vedo_picture
    if o.t <= o.T_total:
        o, f = iteration_func(o, f)
    if vedo_picture:
        fig_view = draw_vedo_and_save(o, o.iter, fig_view, app_diagram=False)

def button_func():
    global timerId
    fig_view.timer_callback("destroy", timerId)
    if "Play" in button.status():
        timerId = fig_view.timer_callback("create", dt=1)
    button.switch()

def run_it_all(o, vedo_picture):
    global timerId, fig_view, button, start_time, collision_foo, f, evnetId
    f = open('storage/main.txt', 'a')
    f.write(f'------------------------------------------------------------\n')
    collision_foo = None  # [None, 'Stop', 'Line']
    start_time = datetime.now()
    if vedo_picture:
        timerId = 1
        fig_view = Plotter(bg='bb', size=(1920, 1080))
        button = fig_view.add_button(button_func, states=["Play ", "Pause"], size=20,
                                     font='Bongas', bold=True, pos=[0.9, 0.9])
        fig_view.timer_callback("destroy", timerId)
        evnetId = fig_view.add_callback("timer", iteration_timer)

        my_mesh = plot_iterations_new(o).color("silver")
        app_mesh = plot_apps_new(o)
        fig_view.show(__doc__, my_mesh + app_mesh, zoom=0.5)

    else:
        while True:
            iteration_timer()
    f.close()


if __name__ == "__main__":
    f = open('storage/main.txt', 'w')
    f.close()
    run_it_all(o, vedo_picture)
