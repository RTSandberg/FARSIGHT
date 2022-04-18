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
    TANH = 0
    SINE = 1
    LOGISTIC = 2
class Quadrature(IntEnum):
    TRAP = 0
    SIMPSONS = 1


ics_type_to_flim_dict = {ICsType.WEAK_LD : (0, 0.44),
                ICsType.STRONG_LD : (0, 0.47),
                ICsType.STRONG_TWO_STREAM : (0, 0.45),
                ICsType.COLDER_TWO_STREAM : (0, 0.8)}

def ics_type_to_flim(ics_type):
    if ics_type in [ICsType.WEAK_LD, ICsType.STRONG_LD, ICsType.STRONG_TWO_STREAM,ICsType.COLDER_TWO_STREAM ]:
        return ics_type_to_flim_dict[ics_type]
    else:
        return (0,1)

