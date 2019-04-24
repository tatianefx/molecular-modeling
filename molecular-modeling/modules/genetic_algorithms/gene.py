#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "tatianefx"

from modules.chemistry.bond import Bond
from modules.common.helper import calculates_angle_a_b_c


class Gene:

    def __init__(self, phi: Bond, psi: Bond):
        self.phi = phi
        self.psi = psi
        self.bond_angle = self.__get_bond_angle()

    """       
                   A
                  . 
         phi ->  .  
                .
               .
              B  ) bond angle
               .
                .
         psi ->  .
                  .
                   C   
    """
    def __get_bond_angle(self):
        a = self.phi.first_atom.position
        b = self.phi.second_atom.position
        c = self.psi.second_atom.position
        return calculates_angle_a_b_c(a, b, c)

