from mylibs.result_show_functions import *

# reader_avoid_field_params_search(filename='storage/pid_const5_avoiding__222.txt')
# plot_avoid_field_params_search()
# plot_repulsion_error()

# reader_pd_control_params()
# pd_control_params_search()

# plot_repulsion_error()
# reader_repulsion_error()

# dv_col_noncol_difference(name='twist_0001', w_twist=0.0001, dt=10, t_max=500, u_max=0.2, times=6)
# reader_dv_col_noncol_difference(name='twist_0001', plot_name=', закрутка 0 рад/с')
'''for w_t in ['00001', '0001', '0002', '0005']:
    print(f"Пошло поехало: {w_t}")
    dv_col_noncol_difference(name='twist_' + w_t, w_twist=float('0.' + w_t), dt=10., t_max=4000., u_max=0.002, times=5)'''

# dv_from_w_twist(dt=10., t_max=4000., u_max=0.002, times=5)
# reader_dv_from_w_twist(plot_name='dt=10, t_max=4000, u_max=0.002')

plot_params_while_main(dt=1., show_rate=100, limit=5000, show_probe_episodes=True)

# crawling()
