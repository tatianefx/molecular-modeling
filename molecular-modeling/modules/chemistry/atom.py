from enum import Enum


class AtomName(Enum):
    C = 'C'
    H = 'H'
    O = 'O'
    N = 'N'
    S = 'S'


class Atom:

    def __init__(self, acronym, position):
        self.acronym = acronym
        self.position = position

