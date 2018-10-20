from math import sqrt, acos, pi


def calculeDistanceAB(a, b):
    x = (a[0] - b[0])**2
    y = (a[1] - b[1])**2
    z = (a[2] - b[2])**2

    return sqrt(x + y + z)


def calculeAngleABC(a, b, c):
    ab = [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    bc = [c[0] - b[0], c[1] - b[1], c[2] - b[2]]

    abVec = sqrt(ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2])
    bcVec = sqrt(bc[0] * bc[0] + bc[1] * bc[1] + bc[2] * bc[2])

    abNorm = [ab[0] / abVec, ab[1] / abVec, ab[2] / abVec]
    bcNorm = [bc[0] / bcVec, bc[1] / bcVec, bc[2] / bcVec]

    res = abNorm[0] * bcNorm[0] + abNorm[1] * bcNorm[1] + abNorm[2] * bcNorm[2]

    return acos(res) * 180.0 / pi

