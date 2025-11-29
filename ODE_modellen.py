############################
## Tumor modelling class  ##
## Auteur: Floris M       ##
## Datum: 26-11-2025      ##
############################

from matplotlib import pyplot as plt
import numpy as np
import math
from math import log

class tumorODE:
    """
    Klasse om tumoren te modelleren met behulp van een aantal groeimodellen.
    
    """

    
    def __init__(self, volume, n, delta_t):
        self.volume = volume
        self.n = n
        self.delta_t = delta_t
        
    def __str__(self):
        return f"Start volume: {self.volume}, aantal dagen (n): {self.n}, tijdsstapgrootte: {self.delta_t}"
        

    def linear_model(self, c, mode=None):
        # Dv/Dt = c
        if mode == "equation":
            return lambda t, v: c

        Ts = [0]
        Vs = [self.volume]
        t = 0

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def exponentieel_model(self, c, mode=None):
        # Dv/Dt = c * V
        if mode == "equation":
            return lambda t, v: c * v

        Ts = [0]
        Vs = [self.volume]
        t = 0

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * self.volume * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)

        return Ts, Vs


    def mendelsohn_model(self, c, d, mode=None):
        # Dv/Dt = c * V^d
        if mode == "equation":
            return lambda t, v: c * math.pow(v, d)

        Ts = [0]
        Vs = [self.volume]
        t = 0

        # Werkt nog niet... math range error

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * math.pow(self.volume, d) * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def linear_gelimiteerdegroei_model(self, c, d, mode=None):
        # Dv/Dt = c * V / (V + d)
        if mode == "equation":
            return lambda t, v: c * (v / (v + d))

        Ts = [0]
        Vs = [self.volume]
        t = 0

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (self.volume / (self.volume + d)) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def oppervlak_gelimiteerdegroei_model(self, c, d, mode=None):
        # Dv/Dt = c * V / V^d
        if mode == "equation":
            return lambda t, v: c * (v / (math.pow(v, d)))

        Ts = [0]
        Vs = [self.volume]
        t = 0

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * (self.volume / (math.pow(self.volume,d))) * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)

        return Ts, Vs
    

    def exponentieel_afvlakkendegroei_model(self, c, v_max, mode=None):
        # Dv/Dt = c * (Vmax - V)
        if mode == "equation":
            return lambda t, v: c * (v_max - v)

        Ts = [0]
        Vs = [self.volume]
        t = 0

        # Dv/Dt = c * (Vmax - V)

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (v_max - self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def von_bertalanffy_model(self, c, d, mode=None):
        # Dv/Dt = c * V^2/3 - d * V
        if mode == "equation":
            return lambda t, v: c * math.pow(v, 2/3) - d * v

        Ts = [0]
        Vs = [self.volume]
        t = 0

        # Math domain error wanneer er een tijdsstapgrootte van 1 gebruikt wordt, werkt wel met 0.1.

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = (c * math.pow(self.volume, 2/3) - d * self.volume) * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)

        return Ts, Vs


    def allee_effect_groei(self, c, v_min, v_max, mode=None):
        # Dv/Dt = c * (V - Vmin) * (Vmax - V)
        if mode == "equation":
            return lambda t, v: c * (v - v_min) * (v_max - v)

        Ts = [0]
        Vs = [self.volume]
        t = 0

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * (self.volume - v_min) * (v_max - self.volume) * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def gompertz_groei(self, c, v_max, mode=None):
        # Dv/Dt = c * V * ln(Vmax / V)
        if mode == "equation":
            return lambda t, v: c * v * log(v_max / v)

        Ts = [0]
        Vs = [self.volume]
        t = 0


        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * self.volume * log(v_max / self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs
    

    def solve(self, method, model):
        t_list = [0]
        v_list = [self.volume]
        t = 0
        y = self.volume
        
        if method == "Euler":
            for _ in range(self.n):
                # Update volgens Euler's methode:
                dy = model(t, y)
                t = t + self.delta_t
                y = y + dy * self.delta_t

                t_list.append(t)
                v_list.append(y)

        elif method == "Heun":
            for _ in range(self.n):
                # Eerste update:
                dydt1 = model(t, y)
                t1 = t + self.delta_t
                y1 = y + dydt1 * self.delta_t
                # Tweede update:
                dydt2 = model(t1, y1)
                t = t + self.delta_t
                y = y + (dydt1 + dydt2) / 2.0 * self.delta_t
               
                t_list.append(t)
                v_list.append(y)

        elif method == "Runge-Kutta":
            for _ in range(self.n):
                # Eerste update:
                dydt1 = model(t, y)
                t1 = t + 0.5 * self.delta_t
                y1 = y + 0.5 * dydt1 * self.delta_t
                # Tweede update:
                dydt2 = model(t1, y1)
                t2 = t + 0.5 * self.delta_t
                y2 = y + 0.5 * dydt2 * self.delta_t
                # Derde update:
                dydt3 = model(t2, y2)
                t3 = t + self.delta_t
                y3 = y + dydt2 * self.delta_t
                # Vierde update:
                dydt4 = model(t3, y3)
                t = t + self.delta_t
                y = y + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * self.delta_t

                t_list.append(t)
                v_list.append(y)

    

        else:
            print("Onbekende solver.... probeer: 'Euler', 'Heun' of 'Runge-Kutta'.")

        return t_list, v_list
    
    def compute_curve(self,a, b, y0):

        def ODE(t, y):
            return a * y + b
        
        t, y = 0.0, y0
        ts, ys = [t], [y]
        for _ in range(self.n):
            # Eerste update
            dydt1 = ODE(t, y)
            t1 = t + 0.5 * self.delta_t
            y1 = y + 0.5 * dydt1 * self.delta_t
            # Tweede update
            dydt2 = ODE(t1, y1)
            t2 = t + 0.5 * self.delta_t
            y2 = y + 0.5 * dydt2 * self.delta_t
            # Derde update
            dydt3 = ODE(t2, y2)
            t3 = t + self.delta_t
            y3 = y + dydt2 * self.delta_t
            # Vierde update volgens Runge-Kutta
            dydt4 = ODE(t3, y3)
            t = t + self.delta_t
            y = y + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * self.delta_t
            # Bewaar tussenresultaten
            ts.append(t)
            ys.append(y) 
        return ts, ys


    # To-do: verder afmaken...
    def MSE(self,a, b, y0):
        _, ys_b = self.compute_curve(a, b, y0)
        sum_squared_error = 0.0
        for y_exact, y_model in zip(ys_a, ys_b):
            error = y_exact - y_model
            sum_squared_error += error * error
        mean_squared_error = sum_squared_error / len(ys_b)
        return mean_squared_error


    def plot(Ts, Vs, color, label):
        """
        Geef de lijst met tijdspunten en volume waarden die daar bijhoren om een grafiek te krijgen die deze plot. 
        """
        plt.plot(Ts, Vs, color=color, label=label)


    def show_plot():
        """
        Weergeef de grafieken in een enkel venster:
        """
        plt.gcf().set_size_inches(12.8, 4.8)
        plt.axvline(0.0, color='k')
        plt.xlabel("Tijd")
        plt.ylabel("Volume")
        plt.legend()
        plt.show()