from colorama import Fore, Style
import matplotlib.pyplot as plt
from p_tqdm import p_map
import numpy as np
import random
import copy

def diff_evolve_sample(args, j: int, func: any, v, target_p, comp_index: list, chance: float = 0.5, f: float = 1.):
    mutant = v[comp_index[0]].copy() + f * (v[comp_index[1]] - v[comp_index[2]])
    for i in range(len(mutant)):
        if random.uniform(0, 1) < chance:
            mutant[i] = v[j][i].copy()
    target_p = func(v[j], args) if target_p is None else target_p
    target = func(mutant, args)
    v[j] = mutant.copy() if target < target_p else v[j]
    target_p = target if target < target_p else target_p
    return np.append(v[j], target_p)


def diff_evolve(func: any, search_domain: list, *args, **kwargs):
    """Функция дифференциальной эволюции.
    :param func: целевая функция
    :param search_domain: 2-мерный список разброса вектора аргументов: [[v[0]_min, v[0]_max], [v[1]_min, v[1]_max],...]
    :return: v_best: len_vec-мерный список"""
    chance = 0.5 if 'chance' not in kwargs.keys() else kwargs['chance']
    f = 1. if 'f' not in kwargs.keys() else kwargs['f']
    n_vec = 10 if 'n_vec' not in kwargs.keys() else kwargs['n_vec']
    len_vec = 3 if 'len_vec' not in kwargs.keys() else kwargs['len_vec']
    n_times = 5 if 'n_times' not in kwargs.keys() else kwargs['n_times']
    multiprocessing = True if 'multiprocessing' not in kwargs.keys() else kwargs['multiprocessing']
    print_process = False if 'print_process' not in kwargs.keys() else kwargs['print_process']
    lst_errors = []
    # попробовать tuple(search_domain[i])
    v = np.array([np.array([random.uniform(search_domain[i][0], search_domain[i][1]) for i in range(len_vec)])
                  for _ in range(n_vec)])
    v_record = [copy.deepcopy(v)]
    target_prev = [None for _ in range(n_vec)]
    v_best = None
    for i in range(n_times):
        if print_process:
            print(Fore.CYAN + f"Шаг {i + 1}/{n_times} дифференциальной эволюции" + Style.RESET_ALL)
        comp_index = [[] for _ in range(n_vec)]
        for j in range(n_vec):
            complement = list(range(n_vec))
            complement.remove(j)
            for _ in range(3):
                comp_index[j].append(random.choice(complement))
                complement.remove(comp_index[j][len(comp_index[j]) - 1])
        anw = p_map(diff_evolve_sample,
                    [args for _ in range(n_vec)],
                    [j for j in range(n_vec)],
                    [func for _ in range(n_vec)],
                    [v for _ in range(n_vec)],
                    [target_prev[j] for j in range(n_vec)],
                    [comp_index[j] for j in range(n_vec)],
                    [chance for _ in range(n_vec)],
                    [f for _ in range(n_vec)]) if multiprocessing else \
            [diff_evolve_sample(args, j, func, v, target_prev[j], comp_index[j], chance, f) for j in range(n_vec)]
        v = np.array([np.array(anw[j][0:len_vec]) for j in range(n_vec)])
        v_record += [copy.deepcopy(v)]
        target_prev = [anw[j][len_vec] for j in range(n_vec)]
        lst_errors.append(np.min(target_prev))
        v_best = v[np.argmin(target_prev)]
    print(Fore.MAGENTA + f"Ответ: {v_best}" + Style.RESET_ALL)
    plt.plot(range(len(lst_errors)), lst_errors, c='indigo')
    plt.show()
    return v_best, v_record
