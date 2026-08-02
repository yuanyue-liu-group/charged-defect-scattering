"""Readers for the BerkeleyGW polarizability / inverse-dielectric HDF5 files.

`Epsmat` wraps chimat.h5 / chi0mat.h5 (chi^0) and epsmat.h5 / eps0mat.h5
(eps^-1), exposing the G-vector <-> index and q-vector <-> index maps that
the sums in Eq. (12) and Eq. (13) run over.
"""

import h5py
import numpy as np

class Epsmat(object):

    def __init__(self, name_):
        self.name = name_
    def read_epsh5(self):
        self.feps = h5py.File(self.name, 'r')
        self.mat_diagonal = np.array(self.feps['mats/matrix-diagonal'])  # diagonal term of dielectric matrix
        
        self.qpts = np.array(self.feps['eps_header/qpoints/qpts']) # qpoints coordinates
        self.nmtx = np.array(self.feps['/eps_header/gspace/nmtx']) # Number of matrix elements BekerleyGW actually compute for each q-point.
        #self.mat_diagonal[:,10000, ]
        self.G_vec = np.array(self.feps['/mf_header/gspace/components']) # G-vectors in RHO G-space
        self.gind_eps2rho = np.array(self.feps['/eps_header/gspace/gind_eps2rho']) # convert Epsilon G-space to Rho G-space
        self.gind_rho2eps = np.array(self.feps['/eps_header/gspace/gind_rho2eps']) # convert Rho G-space to Epsilon G-space
        self.mat = np.array(self.feps['mats/matrix']) # Matrix elements

        self.G_vec_tuple_list = []
        for i_ in self.G_vec:
            i_ = tuple(i_)
            self.G_vec_tuple_list.append(i_)
        self.G_ind2vec = dict(enumerate(self.G_vec_tuple_list))
        self.G_vec2ind = {v:k for k,v in self.G_ind2vec.items()}

        self.qlist = []
        for i_ in self.qpts:
            i_ = tuple(i_)
            self.qlist.append(i_)
        self.q_ind2vec = dict(enumerate(self.qlist))
        self.q_vec2ind = {v:k for k,v in self.q_ind2vec.items()}
        
    def read_epsh5_mdfy(self):
        self.feps = h5py.File(self.name, 'r+')
        self.mat_diagonal = np.array(self.feps['mats/matrix-diagonal'])  # diagonal term of dielectric matrix
        
        self.qpts = np.array(self.feps['eps_header/qpoints/qpts']) # qpoints coordinates
        self.nmtx = np.array(self.feps['/eps_header/gspace/nmtx']) # Number of matrix elements BekerleyGW actually compute for each q-point.
        #self.mat_diagonal[:,10000, ]
        self.G_vec = np.array(self.feps['/mf_header/gspace/components']) # G-vectors in RHO G-space
        self.gind_eps2rho = np.array(self.feps['/eps_header/gspace/gind_eps2rho']) # convert Epsilon G-space to Rho G-space
        self.gind_rho2eps = np.array(self.feps['/eps_header/gspace/gind_rho2eps']) # convert Rho G-space to Epsilon G-space
        self.mat = np.array(self.feps['mats/matrix']) # Matrix elements

        self.G_vec_tuple_list = []
        for i_ in self.G_vec:
            i_ = tuple(i_)
            self.G_vec_tuple_list.append(i_)
        self.G_ind2vec = dict(enumerate(self.G_vec_tuple_list))
        self.G_vec2ind = {v:k for k,v in self.G_ind2vec.items()}

        self.qlist = []
        for i_ in self.qpts:
            i_ = tuple(i_)
            self.qlist.append(i_)
        self.q_ind2vec = dict(enumerate(self.qlist))
        self.q_vec2ind = {v:k for k,v in self.q_ind2vec.items()}
