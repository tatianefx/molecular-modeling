from modules.chemistry.bond import Bond
from modules.common.helper import calculeAngleABC


class Gene:

    def __init__(self, phi: Bond, psi: Bond):
        self.phi = phi
        self.psi = psi
        self.bond_angle = self.__get_bond_angle()

    def __get_bond_angle(self):
        a = self.phi.first_atom.position
        b = self.phi.second_atom.position
        c = self.psi.second_atom.position
        return calculeAngleABC(a, b, c)

