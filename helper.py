import numpy as np

from mylibs.result_show_functions import *
from mylibs.test_functions import *

'''y1 = [54.5, 1.28, 0.47453]
y2 = [13.998, 16.11, 3.196, 3.21, 3.2618, 3.244, 3.0996, 0.95499, 0.38613]
y3 = [12.74299, 0.4739]
y4 = [17.74797, 1,457, 0.45737]
x = [1, 2, 3]
plt.plot([1, len(y2)], [0.5, 0.5], label='ε')
plt.plot(np.arange(len(y1)) + 1, y1)
plt.plot(np.arange(len(y2)) + 1, y2)
plt.plot(np.arange(len(y3)) + 1, y3)
plt.ylabel("Невязка Δr, м", fontsize=15)
plt.xlabel("Итерации", fontsize=15)
plt.legend()
plt.show()

y1 = [3.2620e+03, 8.6042e+03, 8.2748e+03, 3.0384e+03, 8.4549e+02, 3.4871e+02, 1.6901e+02, 1.1868e+02, 1.1067e+02, 1.1016e+02, 1.1034e+02, 
        1.1053e+02, 1.1063e+02, 1.1059e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1038e+02, 1.1038e+02, 1.1035e+02, 
        1.1035e+02, 1.1034e+02, 1.1019e+02, 1.1001e+02, 1.1001e+02, 1.0996e+02, 1.0982e+02, 1.0958e+02, 1.0846e+02, 1.0783e+02, 1.0670e+02, 
        1.0670e+02, 1.0642e+02, 1.0598e+02, 1.0598e+02, 1.0582e+02, 1.0553e+02, 1.0496e+02, 1.0496e+02, 1.0458e+02, 1.0458e+02, 1.0433e+02, 
        1.0271e+02, 1.0271e+02, 8.6593e+01, 6.1708e+01, 2.6835e+01, 1.3145e+01, 1.1598e+01, 1.1604e+01, 1.1604e+01, 1.1604e+01, 1.1604e+01, 
        7.9902e+00, 4.9761e+00, 3.5065e+00, 3.5065e+00, 2.1466e+00, 1.3726e+00, 9.2869e-01, 9.3122e-01, 7.9038e-01, 7.9038e-01, 7.9038e-01,
        3.7559e-01, 6.7691e-02, 4.0578e-02, 3.9191e-02, 3.6617e-02, 3.6559e-02, 3.6559e-02, 3.6559e-02, 3.6559e-02, 2.4321e-02]
plt.plot([1, len(y1)], [0.5, 0.5], label='ε')
plt.plot(np.arange(len(y1)) + 1, [np.sqrt(y) for y in y1])
plt.ylabel("Невязка Δr, м", fontsize=15)
plt.xlabel("Итерации", fontsize=15)
plt.legend()
plt.show()'''

'''y1 = [3.2620e+03, 8.6042e+03, 8.2748e+03, 3.0384e+03, 8.4549e+02, 3.4871e+02, 1.6901e+02, 1.1868e+02, 1.1067e+02, 1.1016e+02, 1.1034e+02,
        1.1053e+02, 1.1063e+02, 1.1059e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1041e+02, 1.1038e+02, 1.1038e+02, 1.1035e+02,
        1.1035e+02, 1.1034e+02, 1.1019e+02, 1.1001e+02, 1.1001e+02, 1.0996e+02, 1.0982e+02, 1.0958e+02, 1.0846e+02, 1.0783e+02, 1.0670e+02,
        1.0670e+02, 1.0642e+02, 1.0598e+02, 1.0598e+02, 1.0582e+02, 1.0553e+02, 1.0496e+02, 1.0496e+02, 1.0458e+02, 1.0458e+02, 1.0433e+02,
        1.0271e+02, 1.0271e+02, 8.6593e+01, 6.1708e+01, 2.6835e+01, 1.3145e+01, 1.1598e+01, 1.1604e+01, 1.1604e+01, 1.1604e+01, 1.1604e+01,
        7.9902e+00, 4.9761e+00, 3.5065e+00, 3.5065e+00, 2.1466e+00, 1.3726e+00, 9.2869e-01, 9.3122e-01, 7.9038e-01, 7.9038e-01, 7.9038e-01,
        3.7559e-01, 6.7691e-02, 4.0578e-02, 3.9191e-02, 3.6617e-02, 3.6559e-02, 3.6559e-02, 3.6559e-02, 3.6559e-02, 2.4321e-02]
f = open('storage/iteration_docking_scipy.txt', 'w')
for i in range(len(y1)):
        f.write(f"{i} {np.sqrt(y1[i])}\n")
f.close()'''



# plot_avoid_field_params_search()
# reader_avoid_field_params_search(filename='storage/pid_const5_avoiding__2.txt')
# plot_repulsion_error()

# reader_pd_control_params()
# pd_control_params_search()

# plot_repulsion_error()
# reader_repulsion_error()

# dv_col_noncol_difference(name='', w_twist=0, dt=10, t_max=5000, u_max=0.03, times=2)
# reader_dv_col_noncol_difference(name='', plot_name=', закрутка 0 рад/с')
'''for w_t in ['00001', '0001', '0002', '0005']:
    print(f"Пошло поехало: {w_t}")
    dv_col_noncol_difference(name='twist_' + w_t, w_twist=float('0.' + w_t), dt=10., t_max=4000., u_max=0.002, times=5)'''

dt=2.
t_max=5000.
u_max=0.03
times=5
# dv_from_w_twist(dt=dt, t_max=t_max, u_max=u_max, times=times)
# reader_dv_from_w_twist(plot_name=f'dt={dt}, t_max={t_max}, u_max={u_max}')

n_p = 20
n_t = 10
# for u0 in [0.04, 0.05]:
u0 = 0.03
name = str(u0)  # + 'control_Z_-5'
# full_bundle_of_trajectories(name=name, dt=10., t_max=5000, u0=u0, n_p=n_p, n_t=n_t)  # , control=np.array([0., 0., 1e-6]))
# reader_full_bundle_of_trajectories(name=name, n_p=n_p, n_t=n_t)
# full_bundle_of_trajectories_controlled(dt=0.5, t_max=15000, control=1e-5)

n_x = 200
n_y = n_x
name = 'scipy'
# heatmap_function(name=name, n_x=n_x, n_y=n_y, target_toward=False, scipy_meh=True)
# reader_heatmap_function(name=name, max_value=115, n_x=n_x, n_y=n_y)

# test_center_mass(u_max=3e-2, dt=100., T_max=100000.)
# test_stitching_collision_function(n=100)

plot_params_while_main(dt=5., show_rate=10, limit=3000, show_probe_episodes=False, energy_show=True,
                       show_j=False, show_R=True, filename='', propulsion_3d_plot=False)
# get_repilsions("")

# plot_iter_downgrade('_scipy')

'''a = np.linspace(0, 10*np.pi, 100)
c = [np.cos(i) for i in a]
s = [np.sin(i) for i in a]
anw = [my_atan2(c[i], s[i]) for i in range(len(c))]
plt.plot(anw)
plt.show()'''

# crawling()
