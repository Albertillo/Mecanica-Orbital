####################################################################################################################################################
# Gregoriandate.py: Esta función se ha desarrollado para uso académico. No se recomienda que se utilice con otros propósitos.                      #
# La función funciona de la siguiente manera: Dado un día juliano, calculará la fecha gregoriana correspondiente y la imprimirá como string.       #
# Algoritmo obtenido de:                                                                                                                           #
# E. G. Richards, Calendars, Internet, Loc. URL: http://aa.usno.navy.mil/publications/docs/c15_usb_online.pdf                                      #
####################################################################################################################################################

def Gregoriandate(J): #J es el día juliano.
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
    f=J+j+(((4*J+B)//146097)*3)//4+C
    e=r*f+v
    g=e%p//r
    h=u*g+w
    D=1+(h%s)//u
    M=(h//s+m)%n+1
    Y=e//p-y+(n+m-M)//n
    return str(int(D)) + "/" + str(int(M)) + "/" + str(int(Y))
