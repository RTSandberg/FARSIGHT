from enum import IntEnum
class BoundaryConditions(IntEnum):
    PERIODIC = 0
    OPEN = 1
class SimType(IntEnum):
    WEAK_LD = 1
    STRONG_LD = 2
    STRONG_TWO_STREAM = 3
    COLDER_TWO_STREAM = 4
    FRIEDMAN_BEAM = 5
class Quadrature(IntEnum):
    TRAP = 0
    SIMPSONS = 1

sim_type_to_flim = {SimType.WEAK_LD : (0, 0.44),
                SimType.STRONG_LD : (0, 0.47),
                SimType.STRONG_TWO_STREAM : (0, 0.45),
                SimType.COLDER_TWO_STREAM : (0, 0.8)}

