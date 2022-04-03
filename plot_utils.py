import matplotlib.pyplot as plt


def do_plot(xlabel, x, ylabel, y, label, title):
    fig = plt.figure(figsize=(20, 12))
    plt.plot(x, y, label=label)

    plt.title(title)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.legend()
    plt.show()


def do_3d_plot(title, xlabel, x, ylabel, y, zlabel, z):
    fig = plt.figure(figsize=(20, 12))
    ax = plt.axes(projection="3d")

    plt.title(title)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)

    ax.plot3D(x, y, z, "gray")

    plt.show()
