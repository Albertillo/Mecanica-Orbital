###############################################################################################################################################################
# StumpffC(z).py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                   #
# La función funciona de la siguiente manera: Dado un valor z, calcula la función de Stumpff C(z).                                                            #
# Algoritmo obtenido de:                                                                                                                                      #
# Howard D. Curtis, Orbital Mechanics for Engineering Students, First Edition, Elsevier Butterworth-Heinemann, Oxford (UK), 2005.							  #
###############################################################################################################################################################

from math import *

def C(Z):
    if Z>0:
        valorC=(1-cos(sqrt(Z)))/Z
    elif Z<0:
        valorC=(cosh(sqrt(-Z))-1)/(-Z)
    else:
        valorC=1/2
    return valorC
