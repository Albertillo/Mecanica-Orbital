######################################################################################################################################
# Julianday.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.            #
# La función funciona de la siguiente manera: Dada una fecha en formato Día, Mes, Año, calculará el día juliano correspondiente.     #
# Algoritmo obtenido de:                                                                                                             #
# E. G. Richards, Calendars, Internet, Loc. URL: http://aa.usno.navy.mil/publications/docs/c15_usb_online.pdf                        #
######################################################################################################################################

def Julianday(D,M,Y): #D=día (puede contener decimales), M=mes, Y=año.
    y=4716
    j=1401
    m=2
    n=12
    r=4
    p=1461
    q=0
    v=3
    u=5
    s=153
    t=2
    w=2
    A=184
    B=274277
    C=-38
    h=M-m
    g=Y+y-(n-h)//n
    f=(h-1+n)%n
    e=(p*g+q)//r+D-1-j
    J=e+(s*f+t)//u-(3*((g+A)//100))//4-C
    return J
