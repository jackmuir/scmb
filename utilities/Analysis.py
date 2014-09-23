import numpy as np
import matplotlib.pyplot as plt
from pandas.tools.plotting import autocorrelation_plot
from pandas import Series
from ylm import Ylmr as ylmr

topo = np.loadtxt('topo_parameters_l_4.dat')

tomo = np.loadtxt('tomo_parameters_l_4.dat')


def slashandburn(data,burn,thin):
    return data[burn::thin]

def plottrace(data):
    plt.figure()
    ax = plt.subplot(1,1,1)
    ax.plot(data)
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Value')
    plt.show()


def autocorrelation(dat):
    plt.figure()
    autocorrelation_plot(Series(dat))
    plt.show()



def histogram(data):
    s = Series(data)
    plt.figure()
    s.hist(color='k', alpha=0.5, bins=50)
    plt.show()
