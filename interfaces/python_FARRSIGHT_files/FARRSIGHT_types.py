from enum import IntEnum
class BoundaryConditions(IntEnum):
    PERIODIC = 0
    OPEN = 1
class ICsType(IntEnum):
    WEAK_LD = 1
    STRONG_LD = 2
    STRONG_TWO_STREAM = 3
    COLDER_TWO_STREAM = 4
    FRIEDMAN_BEAM = 5
    MAXWELL_JUTNER = 6
    RELATIVISTIC_TWO_STREAM = 7
    RELATIVISTIC_WAVE = 8
class ExternalFieldType(IntEnum):
    SINE = 0
    LOGISTIC = 1
class Quadrature(IntEnum):
    TRAP = 0
    SIMPSONS = 1

ics_type_to_flim = {ICsType.WEAK_LD : (0, 0.44),
                ICsType.STRONG_LD : (0, 0.47),
                ICsType.STRONG_TWO_STREAM : (0, 0.45),
                ICsType.COLDER_TWO_STREAM : (0, 0.8)}

