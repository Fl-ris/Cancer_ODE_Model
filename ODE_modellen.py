############################
## Tumor modelling class  ##
## Autheur: Floris M      ##
## Datum: 26-11-2025      ##
############################

from matplotlib import pyplot as plt
import numpy as np
import math

class tumorODE:
    
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

        # Dv/Dt = c * V/V+d

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

        # Dv/Dt = c * V/V+d

        for _ in range(self.n):
            t = t + self.delta_t

            delta_volume = c * (v_max - self.volume) * self.delta_t

            self.volume = delta_volume + self.volume

            Ts.append(t)
            Vs.append(self.volume)


        return Ts, Vs

    def plot(Ts, Vs):
        plt.figure(figsize=(12.8, 4.8))
        plt.axvline(0.0, color='k')
        plt.plot(Ts, Vs, label="test")
        plt.show() 