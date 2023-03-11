# Standard libraries
import numpy as np
import pandas as pd
import copy


def package_beams(N, h):
    ast = 2 * h / np.sqrt(3)
    x = np.linspace(0., 5., N)
    y = np.linspace(0., 5., N)
    flag = 0
    count = 0
    max_count = 0
    for i in range(N - 1):
        if flag > 6:
            flag = 0
            max_count += 1
        if flag == 0:  # вверх, новый уровень
            x[i + 1] = x[i]
            y[i + 1] = y[i] + ast
            count = 0
        if flag == 1:  # влево вниз
            x[i + 1] = x[i] - h
            y[i + 1] = y[i] - ast / 2
        if flag == 2:  # вниз
            x[i + 1] = x[i]
            y[i + 1] = y[i] - ast
        if flag == 3:  # вправо вниз
            x[i + 1] = x[i] + h
            y[i + 1] = y[i] - ast / 2
        if flag == 4:  # вправо вверх
            x[i + 1] = x[i] + h
            y[i + 1] = y[i] + ast / 2
        if flag == 5:  # вверх
            x[i + 1] = x[i]
            y[i + 1] = y[i] + ast
        if flag == 6:  # влево вверх
            if count > 0:
                x[i + 1] = x[i] - h
                y[i + 1] = y[i] + ast / 2
            else:
                flag = 0
                max_count += 1
                x[i + 1] = (x[i] - h)
                y[i + 1] = (y[i] + ast / 2) + ast

        if count == 0:
            flag += 1
            count = max_count
        else:
            count -= 1
    return x, y, max_count + 2


class Structure(object):
    def __init__(self, choice: str = '1', complete: bool = False, floor: int = 5, mass_per_length: float = 1.):
        if floor < 1:
            raise "Поменяй параметры конструкции: floor должен быть равным 1 или больше"
        self.choice = choice
        self.container_length = 12.
        self.h = 0.15
        self.lvl = 0
        self.x_start = 0.6

        if choice == '0':
            self.n_beams = 1
            self.n_nodes = 2
            self.mass = np.array([10.])
            self.length = np.array([10.])
            self.id_node = np.array([np.array([0, 1])])
            self.r1 = np.array([np.array([-5., 0., 0.])])
            self.r2 = np.array([np.array([5., 0., 0.])])
            self.flag = np.array([np.array([1, 1])])
            self.r_st = np.array(np.zeros(3))

        if choice == '1':
            self.n_beams = 24
            self.n_nodes = 10
            self.id = np.arange(self.n_beams)
            self.mass = np.array([5.] * self.n_beams)
            l_beam = 5.0
            self.h = 0.5
            r_inscribed = l_beam / 2 / np.sqrt(3)
            r_circumscribed = l_beam / np.sqrt(3)
            x_ground_floor = l_beam * np.sqrt(2 / 3)
            ast = 2 * self.h / np.sqrt(3)
            r = np.array([np.array([0., 0., 0.]),
                          np.array([x_ground_floor, 0., r_circumscribed]),
                          np.array([x_ground_floor, -l_beam / 2, -r_inscribed]),
                          np.array([x_ground_floor, l_beam / 2, r_inscribed]),
                          np.array([x_ground_floor + l_beam, 0., r_circumscribed]),
                          np.array([x_ground_floor + l_beam, -l_beam / 2, -r_inscribed]),
                          np.array([x_ground_floor + l_beam, l_beam / 2, -r_inscribed]),
                          np.array([x_ground_floor + 2 * l_beam, 0., r_circumscribed]),
                          np.array([x_ground_floor + 2 * l_beam, -l_beam / 2, -r_inscribed]),
                          np.array([x_ground_floor + 2 * l_beam, l_beam / 2, -r_inscribed])])
            id_1 = [0, 0, 0, 1, 1, 2, 1, 2, 3, 2, 1, 3, 4, 4, 5, 4, 5, 6, 5, 4, 6, 7, 7, 8]
            id_2 = [1, 2, 3, 2, 3, 3, 4, 5, 6, 4, 6, 5, 5, 6, 6, 7, 8, 9, 7, 9, 8, 8, 9, 9]
            flag_1 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            flag_2 = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            x_st = [-self.x_start for _ in range(self.n_beams)]
            y_st = [2 * self.h, 3 * self.h, 3 * self.h, 2 * self.h, 1 * self.h, -1 * self.h, -2 * self.h, -2 * self.h,
                    -2 * self.h, -1 * self.h, 0, 1 * self.h, 2 * self.h,
                    2 * self.h, 2 * self.h, 1 * self.h, 0, -1 * self.h, -1 * self.h, 0, 1 * self.h, 1 * self.h, 0, 0]
            z_st = [-2 * ast, -0.5 * ast, 0.5 * ast, 2 * ast, 2.5 * ast, 1.5 * ast, 1 * ast, 0, -1 * ast,
                    -1.5 * ast, -2 * ast, -1.5 * ast, -1 * ast, 0, 1 * ast, 1.5 * ast, 2 * ast, 0.5 * ast,
                    -0.5 * ast, -1 * ast, -0.5 * ast, 0.5 * ast, 1 * ast, 0]
            self.id_node = np.array([np.array([id_1[i], id_2[i]]) for i in range(self.n_beams)])
            self.r1 = np.array([np.array([r[self.id_node[i][0]][0], r[self.id_node[i][0]][1],
                                          r[self.id_node[i][0]][2]]) for i in range(self.n_beams)])
            self.r2 = np.array([np.array([r[self.id_node[i][1]][0], r[self.id_node[i][1]][1],
                                          r[self.id_node[i][1]][2]]) for i in range(self.n_beams)])
            self.flag = np.array([np.array([flag_1[i], flag_2[i]]) for i in range(self.n_beams)])
            self.r_st = np.array([np.array([x_st[i], y_st[i], z_st[i]]) for i in range(self.n_beams)])
            self.length = np.array([np.linalg.norm(self.r1[i] - self.r2[i]) for i in range(self.n_beams)])

        if choice == '2':
            self.n_beams = 3 + floor * 3 + (floor - 1) * 6
            self.n_nodes = 1 + 3 * floor
            self.container_length = 8.
            self.h = 0.25
            a_beam = 5.0
            r_inscribed = a_beam / 2 / np.sqrt(3)
            r_circumscribed = a_beam / np.sqrt(3)
            x_ground_floor = a_beam * np.sqrt(2 / 3)

            y_st, z_st, self.lvl = package_beams(self.n_beams, self.h)
            r = np.array([np.zeros(3) for _ in range(self.n_nodes)])  # Beam coordinates; r[0,:] = [x,y,z]
            r[1] = np.array([x_ground_floor, 0., r_circumscribed])
            r[2] = np.array([x_ground_floor, -a_beam / 2, -r_inscribed])
            r[3] = np.array([x_ground_floor, a_beam / 2, -r_inscribed])
            r1 = [r[0], r[0], r[0], r[1], r[1], r[2]]
            r2 = [r[1], r[2], r[3], r[2], r[3], r[3]]
            id_node_1 = [0, 0, 0, 1, 1, 2]
            id_node_2 = [1, 2, 3, 2, 3, 3]
            sequence = [[1, 2, 3, 1, 2, 3, 4, 5, 6], [4, 5, 6, 6, 4, 5, 6, 4, 5]]
            for i in range(floor - 1):
                r[3 * i + 4] = np.array([x_ground_floor + (i + 1) * a_beam, 0., r_circumscribed])
                r[3 * i + 5] = np.array([x_ground_floor + (i + 1) * a_beam, -a_beam / 2, -r_inscribed])
                r[3 * i + 6] = np.array([x_ground_floor + (i + 1) * a_beam, a_beam / 2, -r_inscribed])
                for j in range(len(sequence[0])):
                    id_1, id_2 = 3 * i + sequence[0][j], 3 * i + sequence[1][j]
                    r1.append(r[id_1])
                    r2.append(r[id_2])
                    id_node_1.append(id_1)
                    id_node_2.append(id_2)

            self.id = np.arange(self.n_beams)
            self.id_node = np.array([np.array([id_node_1[i], id_node_2[i]]) for i in range(self.n_beams)])
            self.r1 = np.array(r1)
            self.r2 = np.array(r2)
            self.flag = np.array([np.array([int(complete), int(complete)]) for _ in range(self.n_beams)])
            self.length = np.array([np.linalg.norm(np.array(self.r1[i]) - np.array(self.r2[i]))
                                    for i in range(self.n_beams)])
            self.r_st = np.array([np.array([-self.x_start - self.container_length + self.length[i], y_st[i], z_st[i]])
                                  for i in range(self.n_beams)])
            self.mass = self.length * mass_per_length

        if choice == '3':
            length = 5
            self.n_beams = 84  # Beams
            self.n_nodes = 32  # Nodes
            self.x_start = 0.4
            R = length * 6
            a = 2 * np.arcsin(length / 2 / R)
            n_around = 3
            angle = np.linspace(a, a * n_around, n_around)
            r = [[R * np.sin(angle[i]) for i in range(n_around)],
                 [R * (1 - np.cos(angle[i])) for i in range(n_around)],
                 [0. for _ in range(n_around)]]
            for j in range(5):
                r[0] += r[0][0:3]
                r[1] += [r[1][i] * np.cos(2 * np.pi * (j + 1) / 6) for i in range(3)]
                r[2] += [r[1][i] * np.sin(2 * np.pi * (j + 1) / 6) for i in range(3)]
            for j in range(6):
                r[0] += r[0][1:3]
                r[1] += [r[1][i] * np.cos(2 * np.pi * (j + 1 / 2) / 6) for i in [1, 2]]
                r[2] += [r[1][i] * np.sin(2 * np.pi * (j + 1 / 2) / 6) for i in [1, 2]]
            r[0] += [0., 10.]
            r[1] += [0., 0.]
            r[2] += [0., 0.]
            id_node_1 = [30, 30, 30, 30, 30, 30, 15, 0, 3, 6, 9, 12, 15, 15, 0, 0, 0, 3, 3, 3, 6, 6, 6, 9, 9, 9, 12, 12,
                         12, 15, 26, 16, 28, 1, 18, 4, 20, 7, 22, 10, 24, 13, 26, 16, 16, 16, 28, 1, 1, 1, 18, 4, 4, 4,
                         20, 7, 7, 7, 22, 10, 10, 10, 24, 13, 13, 13, 27, 17, 29, 2, 19, 5, 21, 8, 23, 11, 25, 14, 0, 3,
                         6, 9, 12, 15]
            id_node_2 = [0, 3, 6, 9, 12, 15, 0, 3, 6, 9, 12, 15, 16, 28, 28, 1, 18, 18, 4, 20, 20, 7, 22, 22, 10, 24,
                         24, 13, 26, 26, 16, 28, 1, 18, 4, 20, 7, 22, 10, 24, 13, 26, 27, 27, 17, 29, 29, 29, 2, 19, 19,
                         19, 5, 21, 21, 21, 8, 23, 23, 23, 11, 25, 25, 25, 14, 27, 17, 29, 2, 19, 5, 21, 8, 23, 11, 25,
                         14, 27, 31, 31, 31, 31, 31, 31]
            y_st, z_st, self.lvl = package_beams(self.n_beams, self.h)

            self.id = np.arange(self.n_beams)
            self.id_node = np.array([np.array([id_node_1[i], id_node_2[i]]) for i in range(self.n_beams)])
            self.r1 = np.array([np.array([r[0][id_node_1[i]], r[1][id_node_1[i]], r[2][id_node_1[i]]])
                                for i in range(self.n_beams)])
            self.r2 = np.array([np.array([r[0][id_node_2[i]], r[1][id_node_2[i]], r[2][id_node_2[i]]])
                                for i in range(self.n_beams)])

            self.flag = np.array([np.array([int(complete), int(complete)]) for _ in range(self.n_beams)])
            self.length = np.array([np.linalg.norm(self.r1[i] - self.r2[i]) for i in range(self.n_beams)])
            self.mass = self.length * mass_per_length
            self.r_st = np.array([np.array([-self.x_start - self.container_length + self.length[i], y_st[i], z_st[i]])
                                  for i in range(self.n_beams)])

        if choice == '4':
            length = 5.
            beam_length_multiplier = 2
            self.n_beams = 16
            self.h = 0.15
            self.container_length = 20.
            self.id = range(self.n_beams)
            self.id_node_1 = [0, 0, 0, 0,
                              1, 1, 2, 3,
                              0, 0, 0, 0,
                              5 + floor * 4, 5 + floor * 4, 6 + floor * 4, 7 + floor * 4]
            self.id_node_2 = [1, 2, 3, 4,
                              4, 2, 3, 4,
                              5 + floor * 4, 6 + floor * 4, 7 + floor * 4, 8 + floor * 4,
                              8 + floor * 4, 6 + floor * 4, 7 + floor * 4, 8 + floor * 4]
            self.r1 = [np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3),
                       np.array([0., -length / 2, length]), np.array([0., length / 2, length]),
                       np.array([0., length / 2, length]), np.array([0., -length / 2, length]),
                       np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3),
                       np.array([0., -length / 2, -length]), np.array([0., length / 2, -length]),
                       np.array([0., length / 2, -length]), np.array([0., -length / 2, -length])]
            self.r2 = [np.array([0., -length / 2, length]), np.array([0., length / 2, length]),
                       np.array([length, -length / 2, length]), np.array([length, length / 2, length]),
                       np.array([length, -length / 2, length]), np.array([0., -length / 2, length]),
                       np.array([length, length / 2, length]), np.array([length, length / 2, length]),
                       np.array([0., -length / 2, -length]), np.array([0., length / 2, -length]),
                       np.array([length, length / 2, -length]), np.array([length, -length / 2, -length]),
                       np.array([length, -length / 2, -length]), np.array([0., -length / 2, -length]),
                       np.array([length, length / 2, -length]), np.array([length, length / 2, -length])]

            def lcl_f(x: float, y: float, z: float):
                global length, i
                return np.array([x * length, length * (y - 1 / 2), (i + 1 + z) * length])

            for i in range(floor):
                for j in range(12):
                    self.id += [self.n_beams + j + 12 * i + 1]
                self.n_beams += 12
                self.id_node_1 += [k + i * 4 for k in [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 5]]
                self.id_node_2 += [k + i * 4 for k in [5, 6, 7, 8, 8, 5, 6, 7, 6, 7, 8, 8]]
                self.r1 += [lcl_f(0, 0, 0), lcl_f(0, 1, 0), lcl_f(1, 1, 0), lcl_f(1, 0, 0),
                            lcl_f(0, 0, 0), lcl_f(0, 1, 0), lcl_f(1, 1, 0), lcl_f(1, 0, 0),
                            lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 1, 1), lcl_f(0, 0, 1)]
                self.r2 += [lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 0, 1), lcl_f(1, 1, 1),  # в конце (1, 0, 1)?
                            lcl_f(1, 0, 1), lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 1, 1),
                            lcl_f(0, 1, 1), lcl_f(1, 1, 1), lcl_f(1, 0, 1), lcl_f(1, 0, 1)]

                for j in range(12):
                    self.id += [self.n_beams + j + 12 * i + 1 + 12 * floor]
                self.n_beams += 12
                self.id_node_1 += [4 * (floor + 1 + i) + k for k in [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 5]]
                self.id_node_1 += [4 * (floor + 1 + i) + k for k in [5, 6, 7, 8, 8, 5, 6, 7, 6, 7, 8, 8]]
                self.r1 += [lcl_f(0, 0, 0), lcl_f(0, 1, 0), lcl_f(1, 1, 0), lcl_f(1, 0, 0),
                            lcl_f(0, 0, 0), lcl_f(0, 1, 0), lcl_f(1, 1, 0), lcl_f(1, 0, 0),
                            lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 1, 1), lcl_f(0, 0, 1)]
                self.r2 += [lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 1, 1), lcl_f(1, 0, 1),
                            lcl_f(1, 0, 1), lcl_f(0, 0, 1), lcl_f(0, 1, 1), lcl_f(1, 1, 1),
                            lcl_f(0, 1, 1), lcl_f(1, 1, 1), lcl_f(1, 0, 1), lcl_f(1, 0, 1)]
                last_id = 8 + i * 4 + (floor + 1) * 4

            r_big_circle = length * (floor + 1)
            length *= beam_length_multiplier
            big_circle_floors = round(2 * np.pi * floor / beam_length_multiplier)
            phi = -np.pi * 2 / big_circle_floors
            rot_nods = np.array([np.array([0., -length / 2, r_big_circle + length]),
                                 np.array([0., length / 2, r_big_circle + length]),
                                 np.array([0., length / 2, r_big_circle]),
                                 np.array([0., -length / 2, r_big_circle])])
            rotation_matrix = np.array([[np.cos(phi), 0., np.sin(phi)],
                                        [0., 1., 0.],
                                        [-np.sin(phi), 0., np.cos(phi)]])

            for i in range(big_circle_floors):
                for j in range(12):
                    self.id_node += [self.n_beams]
                    self.n_beams += 1
                copies = copy.deepcopy(rot_nods)
                for k in range(4):
                    rot_nods[4] = rotation_matrix @ rot_nods[4]
                self.id_node_1 += [k + 1 + (last_id + 1 + i) * 4 for k in [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 5]]
                self.id_node_2 += [k + 1 + (last_id + 1 + i) * 4 for k in [5, 6, 7, 8, 8, 5, 6, 7, 6, 7, 8, 8]]
                self.r1 += [copies[0], copies[1], copies[2], copies[3], copies[0], copies[1], copies[2], copies[3],
                            rot_nods[0], rot_nods[1], rot_nods[2], rot_nods[0]]
                self.r2 += [rot_nods[0], rot_nods[1], rot_nods[2], rot_nods[3], rot_nods[3], rot_nods[0],
                            rot_nods[1], rot_nods[2], rot_nods[1], rot_nods[2], rot_nods[3], rot_nods[3]]

            y_st, z_st, self.lvl = package_beams(self.n_beams, self.h)
            self.id = np.array(self.id)
            self.r1 = np.array(self.r1)
            self.r2 = np.array(self.r2)
            self.id_node = np.array([np.array([self.id_node_1[i], self.id_node_2[2]]) for i in range(self.n_nodes)])
            self.flag = np.array([np.array([int(complete)] * 2) for i in range(self.n_beams)])
            self.length = np.array([np.linalg.norm(self.r1[i] - self.r2[i]) for i in range(self.n_beams)])
            self.mass = self.length * mass_per_length
            self.r_st = np.array([[-self.x_start - self.container_length + self.length[i], y_st[i], z_st[i]]
                                  for i in range(self.n_beams)])

        if choice == '5':
            direction = ['+2', '-0', '+1', '+0', '-2', '-0', '-1', '-2', '+0', '-1', '+0', '+1', '+2', '-1', '+2',
                         '-0', '-2', '-0', '+2', '+2', '+1', '+0', '+1', '+0', '-1', '-2', '+1', '-2', '+1', '-0',
                         '-2', '-1', '-0', '+1', '+2', '-0', '-1', '+2', '+1', '+0', '+2', '+0', '-2', '+0', '+2']
            r1 = []
            r2 = []
            point = np.array([0., 0., 0.])
            for i in range(len(direction)):
                r1 += [point]
                point[int(direction[i][1])] += 4. if direction[i][0] == '+' else -4.
                r2 += [point]

            self.id = np.array([0])
            self.mass = np.array([5.])
            self.length = np.array([1.])
            self.id_node = np.array([[0, 1]])
            self.r1 = np.array([r2[len(r2) - 1]])
            self.r2 = np.array([r2[len(r2) - 1] + [0., 0., 1.]])
            self.flag = np.array([np.zeros(2)])
            self.r_st = np.array([np.zeros(3)])

    def copy(self):
        s = Structure(choice=self.choice)
        s.n_beams = self.n_beams
        s.x_start = self.x_start
        s.n_nodes = self.n_nodes
        s.mass = self.mass.copy()
        s.length = self.length.copy()
        s.id_node = self.id_node.copy()
        s.r1 = self.r1.copy()
        s.r2 = self.r2.copy()
        s.flag = self.flag.copy()
        s.r_st = self.r_st.copy()
        s.container_length = self.container_length
        return s

    def call_possible_transport(self, taken_beams):
        """ Функция осуществляет последовательность сборки"""
        beams_to_take = np.array([])
        mask_non_fixed = [int(np.sum(self.flag[i]) == 0) for i in range(self.n_beams)]

        needed_number_nodes = np.zeros(self.n_nodes)
        current_number_nodes = np.zeros(self.n_nodes)
        mask_open_nodes = np.zeros(self.n_nodes)
        needed_number_nodes[self.id_node[0][0]] = 1   # Костыль на точку стыковки коснтрукции
        current_number_nodes[self.id_node[0][0]] = 1  # с грузовым контейнером

        for i in range(self.n_nodes):
            needed_number_nodes[self.id_node[i][0]] += 1  # Сколько стержней приходят в узел
            needed_number_nodes[self.id_node[i][1]] += 1
            current_number_nodes[self.id_node[i][0]] += self.flag[i][0]  # Сколько стержней в узле находятся
            current_number_nodes[self.id_node[i][1]] += self.flag[i][1]

        for i in range(self.n_nodes):  # В каких узлах неполное кол-во стержней, но есть хоть один
            if (needed_number_nodes[i] - current_number_nodes[i] > 0) and (current_number_nodes[i] > 0):
                mask_open_nodes[i] = 1  # Основная маска

        for i in range(self.n_beams):
            print(f"YOLO {len(self.id_node)}")
            if mask_non_fixed[i] > 0:  # Нетронутые со склада
                if mask_open_nodes[self.id_node[i][0]] + mask_open_nodes[self.id_node[i][1]] > 0:  # Надобность балки
                    beams_to_take = np.append(beams_to_take, self.id[i])

        i = 0
        print(f"beams_to_take: {beams_to_take}, taken:{taken_beams}")
        while i < len(beams_to_take):  # Удалить те, которые уже взяты
            if beams_to_take[i] in taken_beams:
                beams_to_take = np.delete(beams_to_take, i)
                i -= 1
            i += 1
        print([int(i) for i in beams_to_take])

        return [int(i) for i in beams_to_take]


class Container(object):
    def __init__(self, s: Structure, choice='1'):
        self.choice = choice

        if choice == '0':
            self.n = 1
            self.id = np.arange([0])
            self.mass = np.array([5])
            self.diam = np.array([0.5])
            self.r1 = np.array([np.array([5., 0., 0.])])
            self.r2 = np.array([np.array([10., 0., 0.])])
            self.flag_grab = np.array([True])

        if choice == '1':
            self.n = 7
            self.id = np.arange(self.n)
            self.mass = np.array([150., 5., 5., 5., 5., 5., 5.])
            self.diam = np.array([5.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
            self.r1 = np.array([np.array([0.0, 0.0, 0.0]),
                                np.array([-s.x_start / 2, 1.0, 6.0]),
                                np.array([-s.x_start / 2, 1.0, 4.0]),
                                np.array([-s.x_start / 2, -1.0, 4.0]),
                                np.array([-s.x_start / 2, 1.0, -6.0]),
                                np.array([-s.x_start / 2, 1.0, -4.0]),
                                np.array([-s.x_start / 2, -1.0, -4.0])])
            self.r2 = np.array([np.array([-s.x_start, 0.0, 0.0]),
                                np.array([-s.x_start / 2, -1.0, 6.0]),
                                np.array([-s.x_start / 2, 1.0, 6.0]),
                                np.array([-s.x_start / 2, -1.0, 6.0]),
                                np.array([-s.x_start / 2, -1.0, -6.0]),
                                np.array([-s.x_start / 2, 1.0, -6.0]),
                                np.array([-s.x_start / 2, -1.0, -6.0])])
            self.flag_grab = np.array([False, True, False, False, True, False, False])

        def get_handrails(x: float, sequence: list):
            return [np.array([-x, sequence[0][k], sequence[1][k]]) for k in range(3)] + \
                [np.array([-x, sequence[0][k], -sequence[1][k]]) for k in range(3)] + \
                [np.array([-x, -sequence[1][k], sequence[0][k]]) for k in range(3)] + \
                [np.array([-x, sequence[1][k], sequence[0][k]]) for k in range(3)]

        if choice in ['2', '3', '4']:
            r_container = s.h * s.lvl * 1.5
            self.n = 13
            self.id = list(range(self.n))
            self.mass = [r_container] + [1.] * (self.n - 1)
            self.diam = [5.] + [0.1] * (self.n - 1)
            sequence_1 = [[1., 1., -1.], [6., 4., 4.]]
            sequence_2 = [[-1., 1., -1.], [6., 6., 6.]]
            self.r1 = [np.zeros(3)] + get_handrails(x=s.x_start / 2, sequence=sequence_1)
            self.r2 = [np.array([-s.x_start, 0.0, 0.0])] + get_handrails(x=s.x_start / 2, sequence=sequence_2)
            self.flag_grab = [False, True, False, False, True, False, False, True, False, False, True, False, False]

            r_beams = 0.05
            # r_around = 1.5
            r_around = r_container - 0.5
            n_around = round(2 * np.pi * r_container / r_beams * 0.6)
            for i in range(n_around):
                self.id += [self.n]
                self.mass += [0.5]
                self.diam += [r_beams]
                angle = i / n_around * 2 * np.pi
                self.r1 += [[-s.x_start, r_around * np.cos(angle), r_around * np.sin(angle)]]
                self.r2 += [[-s.x_start - s.container_length, r_around * np.cos(angle), r_around * np.sin(angle)]]
                self.flag_grab += [False]
                self.n += 1

            n_crossbar = int(s.lvl * 2)
            rotation_matrix = np.array([[1., 0., 0.], [0., -1 / 2, -np.sqrt(3) / 2], [0., np.sqrt(3) / 2, -1 / 2]])
            for i in range(3 * n_crossbar):
                self.id += [self.id[self.n - 1]]
                self.mass += [0.5]
                self.diam += [r_beams]
                self.flag_grab += [False]
                self.n += 1
            r1 = [np.array([-s.x_start - s.container_length * 0.92, 0., 0.]) for _ in range(n_crossbar)]
            r2 = [np.array([-s.x_start - s.container_length * 0.92, 0., 0.]) for _ in range(n_crossbar)]
            for i in range(s.lvl):
                r1[i][1] = s.h * (i + 1 / 2)
                r1[i][2] = np.sqrt(r_around ** 2 - r1[i][1] ** 2)
                r2[i][1] = s.h * (i + 1 / 2)
                r2[i][2] = -np.sqrt(r_around ** 2 - r2[i][1] ** 2)
                r1[i + s.lvl][1] = - s.h * (i + 1 / 2)
                r1[i + s.lvl][2] = np.sqrt(r_around ** 2 - r1[i + s.lvl][1] ** 2)
                r2[i + s.lvl][1] = - s.h * (i + 1 / 2)
                r2[i + s.lvl][2] = -np.sqrt(r_around ** 2 - r2[i + s.lvl][1] ** 2)
            self.r1 += r1
            self.r2 += r2
            r1 = [rotation_matrix @ r1[i] for i in range(n_crossbar)]
            r2 = [rotation_matrix @ r2[i] for i in range(n_crossbar)]
            self.r1 += r1
            self.r2 += r2
            r1 = [rotation_matrix @ r1[i] for i in range(n_crossbar)]
            r2 = [rotation_matrix @ r2[i] for i in range(n_crossbar)]
            self.r1 += r1
            self.r2 += r2

            self.id = np.array(self.id)
            self.mass = np.array(self.mass)
            self.diam = np.array(self.diam)
            self.r1 = np.array(self.r1)
            self.r2 = np.array(self.r2)
            self.flag_grab = np.array(self.flag_grab)

        if choice == '5':
            direction = ['+2', '-0', '+1', '+0', '-2', '-0', '-1', '-2', '+0', '-1', '+0', '+1', '+2', '-1', '+2',
                         '-0', '-2', '-0', '+2', '+2', '+1', '+0', '+1', '+0', '-1', '-2', '+1', '-2', '+1', '-0',
                         '-2', '-1', '-0', '+1', '+2', '-0', '-1', '+2', '+1', '+0', '+2', '+0', '-2', '+0', '+2']
            self.n = len(direction) + 1
            self.id = np.arange(self.n)
            self.mass = np.array([50.] * self.n)
            self.diam = np.array([0.5] * (self.n - 1) + [0.05])
            r1 = []
            r2 = []
            point = np.array([0., 0., 0.])
            for i in range(len(direction)):
                r1 += [point]
                point[int(direction[i][1])] += 4. if direction[i][0] == '+' else -4.
                r2 += [point]
            self.r1 = np.array(r1 + [r2[len(r2) - 1]])
            self.r2 = np.array(r2 + [r2[len(r2) - 1] + [0., 0., 1.]])
            self.flag_grab = np.array([False] * (self.n - 1) + [True])

    def copy(self):
        c = Container(choice=self.choice)
        c.n = self.n
        c.id = self.id
        c.mass = self.mass
        c.diam = self.diam
        c.r1 = self.r1
        c.r2 = self.r2
        c.flag_grab = self.flag_grab
        return c


class Apparatus(object):
    def __init__(self, X: Structure, n: int = 1, mass: float = 10.):
        if n > len(X.id):
            raise ValueError('Слишком много аппаратов! Получи ошибку!')
        self.X = X
        self.n = n
        self.id = np.arange(n)
        self.mass = np.array([mass] * n)
        self.flag_fly = np.array([False] * n)
        self.flag_start = np.array([True] * n)
        self.flag_beam = np.array([None] * n)
        self.flag_hkw = np.array([True] * n)
        self.flag_to_mid = np.array([True] * n)
        self.busy_time = np.array([i * 40. for i in range(n)])
        self.target_p = np.array([np.zeros(3) for _ in range(n)])
        self.v = np.array([np.zeros(3) for _ in range(n)])
        id_list = X.call_possible_transport([])
        self.target = np.array([(np.array(X.r_st[id_list[i]]) + np.array([-0.3 - X.length[id_list[i]], 0, 0]))
                                for i in range(n)])
        self.r = np.array([(np.array(X.r_st[id_list[i]]) + np.array([-0.3 - X.length[id_list[i]], 0, 0]))
                           for i in range(n)])
        self.r_0 = np.array([mass] * n)

    def copy(self):
        a = Apparatus(X=self.X, n=self.n, mass=self.mass[0])
        a.flag_fly = self.flag_fly.copy()
        a.flag_start = self.flag_start.copy()
        a.flag_beam = self.flag_beam.copy()
        a.flag_hkw = self.flag_hkw.copy()
        a.flag_to_mid = self.flag_to_mid.copy()
        a.busy_time = self.busy_time.copy()
        a.target_p = self.target_p.copy()
        a.target = self.target.copy()
        a.v = self.v.copy()
        a.r = self.r.copy()
        a.r_0 = self.r_0.copy()
        return a


def get_all_components(choice: str = '1'):
    """Функция инициализирует классы конструкции и аппаратов"""
    s = Structure(choice)
    c = Container(s, choice)
    a = Apparatus(s)
    return s, c, a
