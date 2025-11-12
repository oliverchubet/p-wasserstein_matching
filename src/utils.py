import math
from constants import *

def dist(a,b):
    return math.dist(a.coords, b.coords)

def l_infty(a,b):
    return max([abs(x-y) for x,y in zip(a.coords, b.coords)])

def l_1(a,b):
    return sum([abs(x-y) for x,y in zip(a.coords, b.coords)])
    
