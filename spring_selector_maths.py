import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RangeSlider


# model
# all dimentions mm, Newtons and degreese
r = 10
x = {"a": [0, 0],
     "b": [0, 25.6],
     "c": [0, 190],
     "e": [109, 150]}
theta = [3.3, 15.4]
T = [2460, 312]


def xd(index):
    t = np.deg2rad(theta[index])
    xd = np.matmul(np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]]), x["e"])
    return xd


def vec(start, end, index=0):
    y = dict(x)
    y["d"] = xd(index)
    return np.array(y[end]) - np.array(y[start])


def A(index):
    return 2 * r + np.cross(vec("a", "b", index), vec("b", "d", index)) / np.linalg.norm(vec("b", "d", index))


def B(index):
    return np.cross(vec("a", "c", index), vec("c", "d", index)) / np.linalg.norm(vec("c", "d", index))


def e_0():
    X1 = B(0) * A(1) * T[1]
    X2 = B(1) * A(0) * T[0]
    e_0 = (X1 * np.linalg.norm(vec("c", "d", 0)) - X2 * np.linalg.norm(vec("c", "d", 1))) / (X1 - X2)
    return e_0


def spring_constant(index, e_0):
    K = -T[index] * A(index) / (B(index) * (np.linalg.norm(vec("c", "d", index)) - e_0))
    return K


def Moment_about_a(i, K, solved_e_0):
    p1 = 2 * T[i] * r
    print(f"2TR is: {p1}")
    p2 = T[i] * np.cross(vec("a", "b", i), vec("b", "d", i)) / np.linalg.norm(vec("b", "d", i))
    print(f"T(Vab X Vbd^) is: {p2}")
    p3 = K * (1 - solved_e_0 / (np.linalg.norm(vec("c", "d", i)))) * np.cross(vec("a", "c", i), vec("c", "d", i))
    print(f"K term is is: {p3}")
    print(f"M(a) = {p1 + p2 + p3}")
    return p1 + p2 + p3


# view
# plotting to help with debugging/sanity checking
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-50, 200)
ax.set_ylim(-50, 250)
ax.set_aspect('equal')

# fixed points
for point in [x["a"], x["b"]]:
    circle = plt.Circle(point, r, fill=True)
    ax.add_artist(circle)
x_points = []
x_points += ax.plot(x["c"][0], x["c"][1], "bo")
x_points += ax.plot(x["e"][0], x["e"][1], "bo")

# moving features
for i in range(2):
    x_points += ax.plot(xd(i)[0], xd(i)[1], "o")

x["d1"] = list(xd(0))
x["d2"] = list(xd(1))
annotations = []
for point in x:
    annotations += [ax.annotate(point, (x[point]), color="r")]

solution = plt.text(60, .025, f' ')



# sliders
fig.subplots_adjust(left=0.25, bottom=0.25)
sliders = []
for i, point in enumerate(["b", "c", "e"]):
    axl = fig.add_axes([0.04 + 0.08 * i, 0.25, 0.0225, 0.63])
    print(f"x = : {x[point][0]}")
    sliders += [Slider(
        ax=axl,
        label=f"{point}_x",
        valmin=0,
        valmax=2*x[point][0],
        valinit=x[point][0],
        orientation="vertical")]
    axl = fig.add_axes([0.08 + 0.08 * i, 0.25, 0.0225, 0.63])
    sliders += [Slider(
        ax=axl,
        label=f"{point}_y",
        valmin=0,
        valmax=2*x[point][1],
        valinit=x[point][1],
        orientation="vertical")]
for i, var in enumerate([theta, T]):
    axl = fig.add_axes([0.25, 0.04 + 0.04 * i, 0.65, 0.03])
    sliders += [RangeSlider(
        ax=axl,
        label=f"{['theta', 'T'][i]} range",
        valmin=0,
        valmax=2*max(var),
        valinit=var)]

def update(var):
    global theta, T
    for i, point in enumerate(["b", "c", "e"]):
        # print(sliders[int(2*i)].val)
        x[point][0] = sliders[int(2*i)].val
        x[point][1] = sliders[2*i+1].val
    theta = [sliders[-2].val[0], sliders[-2].val[1]]
    T = [sliders[-1].val[1], sliders[-1].val[0]]
    x["d1"] = list(xd(0))
    x["d2"] = list(xd(1))

    x_points[0].set_data(*x["c"])
    x_points[1].set_data(*x["e"])
    x_points[2].set_data(*x["d1"])
    x_points[3].set_data(*x["d2"])
    for i, point in enumerate(x):
        annotations[i].set_position(x[point])
    e0 = e_0()
    k = spring_constant(0, e0)
    solution.set_text(f'e_0 = {e_0()} \n'
                      f'K = {spring_constant(0, e_0())}\n'
                      f'Extension range > {- e_0() + np.linalg.norm(vec("c", "e"))}')

for slider in sliders:
    slider.on_changed(update)

# plot linearity of result maybe

e0 = e_0()
k = spring_constant(0, e0)
print(e0, k)
k = spring_constant(1, e0)
print(e0, k)
# Moment_about_a(0, k, e0)
update(0)
plt.show()

