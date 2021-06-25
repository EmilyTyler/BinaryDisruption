import csv
import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import stats
plt.rc('font', family='serif')
seed=10
np.random.seed(seed)

def eval_y(x, coeffs):
    y = np.zeros(np.size(x))
    for i in range(np.size(coeffs)):
        y += coeffs[i]*(x)**(np.size(coeffs)-1-i)
    if np.size(x) == 1:
        if (y < 0.0):
            y = 0.0
    else:
        y = np.where(y < 0.0, 0.0, y)
    return y

def x_from_area(u, x_min, y_min, y_max):
    return x_min + u/(y_max-y_min)

def comparison(x, x_min, x_max, y_max):
    if np.size(x) == 1:
        c = 0
    else:
        c = np.zeros(np.size(x))
    for i in range(np.size(x)):
        if (np.size(x) == 1):
            if (x_min < x < x_max):
                c = y_max
        else:
            if (x_min < x[i] < x_max):
                c[i] = y_max
    return c

#Function to draw numbers from a polynomial distribution given the polynomial coefficients
def draw_from_poly(coeffs, x_min, x_max, size=None):
    #Find maximum of distribution function
    xs = np.linspace(x_min, x_max, num=100000)
    ys = eval_y(xs, coeffs)
    y_max = np.max(ys)
    y_min = 0.0
    #Total area under comparison function (rectangle)
    area = (x_max-x_min)*(y_max-y_min)
    result = np.zeros(size)
    for i in range(np.size(result)):
        while(True):
            u = np.random.uniform(0.0, area, size=None)
            x = x_from_area(u, x_min, y_min, y_max)
            y = np.random.uniform(y_min, comparison(x, x_min, x_max, y_max))
            if (y < eval_y(x, coeffs)):
                if size == None:
                    result = x
                else:
                    result[i] = x
                break
    return result
