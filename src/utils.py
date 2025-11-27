import math
from math import ceil
from constants import *

def dist(a,b):
    return math.dist(a.coords, b.coords)

def rounded_dist(a,b, delta=0.01):
    return ceil(math.dist(a.coords, b.coords)/delta)*delta

def l_infty(a,b):
    return max([abs(x-y) for x,y in zip(a.coords, b.coords)])

def l_1(a,b):
    return sum([abs(x-y) for x,y in zip(a.coords, b.coords)])
    
