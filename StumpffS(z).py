###############################################################################################################################################################
# StumpffS(z).py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                                   #
# La función funciona de la siguiente manera: Dado un valor z, calcula la función de Stumpff S(z).                                                            #
# Algoritmo obtenido de:                                                                                                                                      #
# Howard D. Curtis, Orbital Mechanics for Engineering Students, First Edition, Elsevier Butterworth-Heinemann, Oxford (UK), 2005.							  #
###############################################################################################################################################################

from math import *

def S(Z):
    if Z>0:
        valorS=(sqrt(Z)-sin(sqrt(Z)))/sqrt(Z*Z*Z)
    elif Z<0:
        valorS=(sinh(sqrt(-Z))-sqrt(-Z))/sqrt((-Z)*(-Z)*(-Z))
    else:
        valorS=1/6
    return valorS
