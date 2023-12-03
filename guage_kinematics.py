import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RangeSlider
from scipy.optimize import bisect



class TensionGaugeMaths:
    def __init__(self):
        # positions (mm), x_pos, y_pos, radius (if aplic
        self.r_a = 12.1/2
        self.r_b = 16.5/2
        self.x = {"a": [0, 0],
                  "b": [2.22, 25.6],
                  "c": [0, 190],
                  "e": [109, 150.62]}
        # limits (deg)
        self.theta_min = 3.3
        self.theta_max = 15.4
        # spring paramiters
        self.spring_constant = 0.95 # N/mm
        self.e_0 = 39.1  #mm
        # wire paramiters
        self.wire_diameter = 3  # mm
        # Scale properties
        self.scale_thetas = [16, 0]
        self.scale_numbers = [10, 40]

    def xd(self, theta):
        t = np.deg2rad(theta)
        xd = np.matmul(np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]]), self.x["e"])
        return xd

    def vec(self, start, end, theta=3.3):
        y = dict(self.x)
        y["d"] = self.xd(theta)
        return np.array(y[end]) - np.array(y[start])

    def moment_about_a(self, theta, wire_tension):
        p1 = wire_tension * (self.r_a + self.r_b + self.wire_diameter)
        p2 = wire_tension * np.cross(self.vec("a", "b"), self.vec("b", "d", theta)) / np.linalg.norm(self.vec("b", "d", theta))
        p3 = self.spring_constant * (1 - self.e_0 / (np.linalg.norm(self.vec("c", "d", theta)))) * np.cross(self.vec("a", "c"), self.vec("c", "d", theta))
        res = p1 + p2 + p3
        return res

    def tension_2_theta(self, tension):
        theta_solution = bisect(self.moment_about_a, -3, 20, args=tension)
        return theta_solution

    def plot(self):
        tension = 250
        thetas = np.linspace(-3, 20, 100)
        moments = [self.moment_about_a(theta, tension) for theta in thetas]
        plt.plot(thetas, moments, "*")
        plt.show()

class TensionGuage:
    def __init__(self):
        self.maths = TensionGaugeMaths()

    def plot_theta_vs_tension(self):
        tensions = np.linspace(300, 3600, 100)  # newtons
        thetas = [self.maths.tension_2_theta(tension) for tension in tensions]
        print(f"scale: {tensions[50]} \n Angle: {thetas[50]}")
        plt.plot(tensions, thetas)
        plt.xlabel("Wire tension (N)")
        plt.ylabel("Theta (deg)")
        plt.grid(True)
        plt.show()

    def plot_theta_vs_tension_loose_scales(self):
        tensions = np.linspace(300, 3600, 100)  # newtons
        tensions_scale_2_5 = Conversions.newtons_to_scale_2_5mm(tensions)
        tensions_scale_3 = Conversions.newtons_to_scale_3mm(tensions)
        tensions_scale_4 = Conversions.newtons_to_scale_4mm(tensions)
        # 2.5mm
        self.maths.wire_diameter = 2.5
        thetas_2_5 = [self.maths.tension_2_theta(tension) for tension in tensions]
        # 3mm
        self.maths.wire_diameter = 3
        thetas_3 = [self.maths.tension_2_theta(tension) for tension in tensions]
        # 4mm
        self.maths.wire_diameter = 4
        thetas_4 = [self.maths.tension_2_theta(tension) for tension in tensions]
        self.maths.wire_diameter = 3

        plt.plot(tensions_scale_2_5, thetas_2_5, label="2.5mm")
        plt.plot(tensions_scale_3, thetas_3, label="3mm")
        plt.plot(tensions_scale_4, thetas_4, label="4mm")
        plt.plot(self.maths.scale_numbers, self.maths.scale_thetas, label="Tension gauge scale")
        plt.xlabel("Wire tension (scale)")
        plt.ylabel("Theta (deg)")
        plt.legend()
        plt.grid(True)
        plt.show()


class Conversions:
    @staticmethod
    def newtons_to_scale_3mm(tensions):
        # 3mm cable mapping
        scale = np.array([13, 16, 18, 21, 24, 26, 28, 30, 32])
        tension_kgs = np.array([60, 75, 90, 120, 150, 170, 190, 220, 250])
        tension_n = tension_kgs*9.81
        return np.interp(tensions, tension_n, scale)

    @staticmethod
    def newtons_to_scale_2_5mm(tensions):
        # 2.5mm cable mapping
        scale = np.array([5, 8, 10, 13, 16, 18, 21])
        tension_kgs = np.array([33, 50, 58, 70, 90, 110, 140])
        tension_n = tension_kgs*9.81
        return np.interp(tensions, tension_n, scale)

    @staticmethod
    def newtons_to_scale_4mm( tensions):
        # 4mm cable mapping
        scale = np.array([21, 24, 26, 28, 30, 32, 35, 38, 40])
        tension_kgs = np.array([70, 90, 115, 140, 160, 180, 225, 280, 360])
        tension_n = tension_kgs*9.81
        return np.interp(tensions, tension_n, scale)

def main():
    tension_guage = TensionGuage()
    #tension_guage.plot_theta_vs_tension()
    tension_guage.plot_theta_vs_tension_loose_scales()


if __name__ == "__main__":
    main()

