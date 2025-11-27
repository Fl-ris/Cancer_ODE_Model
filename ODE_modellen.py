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
        

    def linear_model(self, c):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c

        for _ in range(self.n):
            t = t + self.delta_t
            delta_volume = c * self.delta_t
            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def exponentieel_model(self, c):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * self.volume * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def mendelsohn_model(self, c, d):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V^d
        # Werkt nog niet... math range error

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * math.pow(self.volume,d) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def linear_gelimiteerdegroei_model(self, c, d):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V / V + d

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (self.volume / (self.volume +d)) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def oppervlak_gelimiteerdegroei_model(self, c, d):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V/V+d

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (self.volume / (math.pow(self.volume,d))) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs
    

    def exponentieel_afvlakkendegroei_model(self, c, v_max):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * (Vmax - V)

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (v_max - self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def von_bertalanffy_model(self, c, d):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V^2/3 - d * V
        # Math domain error wanneer er een tijdsstapgrootte van 1 gebruikt wordt, werkt wel met 0.1.

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * math.pow(self.volume,2/3) - d * self.volume * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def allee_effect_groei(self, c, v_min, v_max):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * (V - Vmin) * (Vmax - V)

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (self.volume - v_min) * (v_max - self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs


    def gompertz_groei(self, c, v_max):
        Ts = []
        Vs = []
        t = 0

        # Dv/Dt = c * V * ln(Vmax / V)

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * self.volume * log(v_max / self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs
    

    def solve(self, method, model):
        t_list = []
        v_list = []
        t = 0
        y = 0
        
        if method == "Euler":
            for _ in range(self.n):
                # Update volgens Euler's methode:
                dy = self.model(t, y)
                t = t + self.delta_t
                y = y + dy * self.delta_t

                t_list.append(t)
                v_list.append(y)

        elif method == "Heun":
            for _ in range(self.n):
                # Eerste update:
                dydt1 = self.model(t, y)
                t1 = t + self.delta_t
                y1 = y + dydt1 * self.delta_t
                # Tweede update:
                dydt2 = self.model(t1, y1)
                t = t + self.delta_t
                y = y + (dydt1 + dydt2) / 2.0 * self.delta_t
               
                t_list.append(t)
                v_list.append(y)

        elif method == "Runge-Kutta":
            for _ in range(self.n):
                # Eerste update:
                dydt1 = self.model(t, y)
                t1 = t + 0.5 * self.delta_t
                y1 = y + 0.5 * dydt1 * self.delta_t
                # Tweede update:
                dydt2 = self.model(t1, y1)
                t2 = t + 0.5 * self.delta_t
                y2 = y + 0.5 * dydt2 * self.delta_t
                # Derde update:
                dydt3 = self.model(t2, y2)
                t3 = t + self.delta_t
                y3 = y + dydt2 * self.delta_t
                # Vierde update:
                dydt4 = self.model(t3, y3)
                t = t + self.delta_t
                y = y + (dydt1 + 2.0 * dydt2 + 2.0 * dydt3 + dydt4) / 6.0 * self.delta_t

                t_list.append(t)
                v_list.append(y)


        else:
            print("Onbekende solver.... probeer: 'Euler', 'Heun' of 'Runge-Kutta'.")


    def plot(Ts, Vs):
        """
        Geef de lijst met tijdspunten en volume waarden die daar bijhoren om een grafiek te krijgen die deze plot. 
        """
        
        plt.figure(figsize=(12.8, 4.8))
        plt.axvline(0.0, color='k')
        plt.plot(Ts, Vs, label="test")
        plt.show() 