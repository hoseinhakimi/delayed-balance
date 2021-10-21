import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations, permutations
from scipy import sparse

class delay:
    def set_system(self, nodes_count, enmity_ratio, delay, d_sign):
        self.nc, self.er, self.delay, self.d_sign = (nodes_count, enmity_ratio,
                                                     delay, d_sign)
        self.initiate()

    def initiate(self):
        self.create_network()
        self.total_tri_energy()

    def create_network(self):
        links = np.array(list(combinations(np.arange(self.nc), 2)))

        self.lc = len(links)
        signs = np.random.choice(a=[-1, 1],
                                 size=self.lc,
                                 p=[self.er, 1 - self.er])
        adj = sparse.coo_matrix((signs, (links[:, 0], links[:, 1])),
                                (self.nc, self.nc)).toarray()
        self.adj_mat = adj + adj.T

    def cal_rows(self, a):
        e = self.adj_mat[a[0], a[1]] \
            * self.adj_mat[a[1], a[2]] \
            * self.adj_mat[a[2], a[0]]

        val = 0 if e == self.d_sign else np.random.randint(-2*self.delay, -1)
        return val

    def total_tri_energy(self):
        self.tc = len(list(combinations(np.arange(self.nc), 3)))
        tri_combinations = np.array(list(combinations(np.arange(self.nc), 3)))
        timess = np.apply_along_axis(self.cal_rows, 1, tri_combinations)
        self.state_array = np.c_[tri_combinations, timess, timess == 0]
        energy_part = np.matmul(self.adj_mat, self.adj_mat)
        energy_mat = np.multiply(self.adj_mat, energy_part)
        self.tri_energy = -np.sum(energy_mat) / 6

    def permute_nodes(self, array):
        return np.array(list(permutations(array, 3))[:3])

    # dynamic
    def balanceTimer(self, link, time):
        indices = list(range(self.nc))
        indices.remove(link[1])
        indices.remove(link[0])

        for i in indices:
            sign = \
                self.adj_mat[link[0], link[1]] * \
                self.adj_mat[link[0], i] * \
                self.adj_mat[i, link[1]]

            tri_nodes = [link[0], link[1], i]
            tri_nodes.sort()
            row = np.argwhere(
                (self.state_array[:, :-2] == tri_nodes).all(axis=1))[0][0]
            prev = self.state_array[row, 4]

            if (sign == self.d_sign):  # if new triangle sign is equal to desired sign
                self.state_array[row, 3] = time
                self.state_array[row, 4] = 1
            elif prev:  # if new sign is not desired check previous state of triangle
                self.state_array[row, 3] = time - 1
                self.state_array[row, 4] = 0

    def checker(self, link, time):
        linkEnergy = 0
        indices = list(range(self.nc))
        indices.remove(link[1])
        indices.remove(link[0])

        for i in indices:
            tri_nodes = [link[0], link[1], i]
            tri_nodes.sort()
            row = np.argwhere(
                (self.state_array[:, :-2] == tri_nodes).all(axis=1))[0][0]
            if (time - self.state_array[row, 3] <= self.delay
                    or self.state_array[row, 4]):
                linkEnergy -= 1 * self.d_sign
            else:
                linkEnergy += 1 * self.d_sign

        return linkEnergy

    def flipper(self, time):
        link = np.random.choice(range(self.nc), 2, replace=False)
        tri_energy_m = -np.inner(  # energy based on normal dynamics
            self.adj_mat[:, link[0]],
            self.adj_mat[:, link[1]]) * self.adj_mat[link[0], link[1]]
        # energy based on delayed dynamics
        tri_energy = self.checker(link, time)

        if tri_energy > 0:
            self.adj_mat[link[0], link[1]] = self.adj_mat[
                link[1], link[0]] = -1 * self.adj_mat[link[0], link[1]]
            self.balanceTimer(link, time)
            self.tri_energy += -2 * tri_energy_m
