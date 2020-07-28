import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import time

def annotate_plot(func):
    def wrapper_annotate_plot(*args, **kwargs):
        func(*args, **kwargs)
        plot_date_time = datetime.strftime(datetime.now(), "%Y/%m/%d %H:%M:%S")
        plt.subplots_adjust(bottom=0.2)
        plt.text(10, 10, 
                f'Generated at {plot_date_time} by {func.__name__}() in {__file__}', 
                transform=None)
    return wrapper_annotate_plot

@annotate_plot
def plot_sine_plt():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    plt.plot(x, y)
    return

@annotate_plot
def plot_sine_subplot():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    fig, ax = plt.subplots()
    ax.plot(x, y)
    return

@annotate_plot
def plot_sine_subplots():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(x, y)
    ax[1].plot(x, y)
    return


if __name__ == "__main__":
    plot_sine_plt()
    time.sleep(1)
    plot_sine_subplot()
    time.sleep(1)
    plot_sine_subplots()
    plt.show()