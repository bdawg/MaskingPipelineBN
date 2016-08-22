# Pythonised by Paul
# PNS DEC2012
#
#
# 05DEC12 - bugfix Paul


def rad2mas(rad):
    from numpy import pi

    mas=rad*180./pi
    mas=mas*3600.*1000.

    return mas


def mas2rad(mas):
    from numpy import pi

    rad=mas/1000.
    rad=rad/3600.
    rad=rad*pi/180.

    return rad

