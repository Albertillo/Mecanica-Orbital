###############################################################################################################################################################
# angulomenor2pi.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                #
# La función funciona de la siguiente manera: Dado un ángulo superior a 2pi lo convertirá a un ángulo comprendido entre 0 y 2pi.							  #
###############################################################################################################################################################

from math import *

def angulomenor2pi(ang): #ang en rad.
    if ang>=2*pi:
        k=ang//(2*pi)
        valor=ang-2*pi*k
    elif ang<2*pi and ang>=0:
        valor=ang
    elif ang<0 and ang>=-2*pi:
        valor=ang+2*pi
    else:
        k=abs(ang)//(2*pi)
        valor=ang+2*pi*(k+1)
    return valor
