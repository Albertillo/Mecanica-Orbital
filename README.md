# Mecanica-Orbital
Software para diferentes cálculos de Mecánica Orbital

Estos códigos en Python se han escrito para los cálculos de mi Trabajo de fin de Grado: Análisis y Optimización de una Misión Espacial al Satélite Europa. Como se han calculado diferentes tipos de trayectorias en el trabajo (Hohmann, no-Hohmann, asistencia gravitatoria múltiple y trayectorias low-thrust), subiré aquí diferentes funciones que pueden ser de ayuda para el cálculo de este tipo de trayectorias. Todos los códigos han sido realizados por mi. Sin embargo, todos los métodos, algoritmos y teoría relacionados se han obtenido de los libros que se mencionan en la introducción que se encuentra en el archivo de cada función.

Las funciones que he subido son las siguientes:\
	· Cálculo de efemérides (resuelve el problema de Kepler, las funciiones son: Kepler, Kepler2 y Kepler3).\
	· Calculadora de día juliano en función de fecha gregoriana (Julianday).\
	· Calculadora de fecha gregoriana en formato string en función de día juliano (Gregoriandate).\
	· Adaptación de ángulos (angulomenor2pi).\
	· Cálculo funciones de Stumpff C(z) y S(z) (StumpffC(z) y StumpffS(z)).\
	· Problema de Lambert (Lambert)\
	· Conversión vector posición y velocidad a elementos orbitales (VectoElem).\
	· Clasificación de planetas con sus elementos orbitales y parámetros de interés (Clasif).\
	· Módulo básico para el cálculo de DeltaV para maniobras tipo Hohmann, no-Hohmann y asistencia gravitatoria con impulso en periapsis. 		(ModBas)\
	· Búsqueda de trayectorias con varias asistencias gravitatorias (BusquedaFlybyMultiple).\
	· Transferencia de Edelbaum (EdelbaumTransfer).\

