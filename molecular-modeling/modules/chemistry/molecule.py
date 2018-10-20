from modules.chemistry.atom import Atom
from modules.chemistry.bond import Bond


class Molecule:

    def __init__(self, atoms: [Atom], bonds: [Bond]):
        self.atoms = atoms
        self.bonds = bonds

