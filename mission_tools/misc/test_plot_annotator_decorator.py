import matplotlib.pyplot as plt
import numpy as np
import time

from plot_annotator_decorator import annotate_plot

@annotate_plot
def plot_sine_plt():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    plt.plot(x, y)
    plt.title('The greatest sine curve in the land')
    plt.xlabel('x')
    plt.ylabel('y')
    return

@annotate_plot
def plot_sine_subplot():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_title('The greatest sine curve in the land')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    return

@annotate_plot
def plot_sine_subplots():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    fig, ax = plt.subplots(2, 1)
    ax[0].plot(x, y)
    ax[1].plot(x, y)
    ax[0].set_title('The greatest sine curves in the land')
    ax[-1].set_xlabel('x')
    ax[-1].set_ylabel('y')
    return


if __name__ == "__main__":
    plot_sine_plt()
    time.sleep(1)
    plot_sine_subplot()
    time.sleep(1)
    plot_sine_subplots()
    plt.show()