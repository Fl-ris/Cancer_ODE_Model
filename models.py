import numpy as np
class tumor_growth_models:
    """ Klasse met diverse modellen om tumor groei te simuleren.
    """

    def __init__(self, volume, delta_t, n):
        """
        Initialiseer het groeimodel.
        Parameters:
            volume (float): volume
            delta_t (float): delta t, tijdsverschil Δt tussen opeenvolgende tijdstappen
            n (int): aantal tijdstappen
        """
        self.volume = volume
        self.delta_t = delta_t
        self.n = n

    def lineaire_model(self, c):
        """
        Simuleert lineaire groei met constante groeisnelheid.
        Differentialvergelijking:
            dV/dt = c

        Parameters:
            c: groeiconstante

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij volume.self
        """

        # Lijsten voor tijdstippen en bijbehorende volumes
        Ts = [0]            # beginpunt t=0
        Vs = [self.volume]  # beginvolume

        t = 0               # huidige tijd
        V = self.volume     # lokale variabele voor volume, zodat self.volume onveranderd blijft

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t   # update tijd met tijdstap delta_t
            V += c * self.delta_t  # Euler-update van het volume

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def exponentieel_model(self, c):
        """
        Simuleert exponentiële groei met groeisnelheid evenredig aan het huidige volume.
        Differentialvergelijking:
            dV/dt = c * V

        Parameters:
            c : groeiconstante

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij volume.self
        """

        # Lijsten voor tijdstippen en volumes
        Ts = [0]            # beginpunt t=0
        Vs = [self.volume]  # beginvolume

        t = 0               # huidige tijd
        V = self.volume     # lokale variabele voor volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t        # update tijd
            V += c * V * self.delta_t  # Euler-update van V (exponentiële groei)

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def mendelsohn_model(self, c, d):
        """
        Simuleert Mendelsohn-groei met groeisnelheid afhankelijk van V^d.
        Differentialvergelijking:
            dV/dt = c * V^d

        Parameters:
            c : groeiconstante
            d : exponent van het volume

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij volume.self
        """

        # Lijsten voor tijdstippen en volumes
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume

        t = 0  # huidige tijd
        V = self.volume  # lokale variabele voor volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t  # update tijd
            V += c * V ** d * self.delta_t  # Euler-update van V (Mendelsohn-groei)

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def exponentieel_afvlakkend_model(self, c, V_max):
        """
        Simuleert exponentieel afvlakkende groei waarbij de groeisnelheid afneemt
        naarmate het volume dichter bij een maximale waarde V_max komt.
        Differentialvergelijking:
        dV/dt = c * (V_max - V)

        Parameters:
            c : groeiconstante
            V_max : maximale waarde waarheen het model groeit

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij volume.self
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t # update tijd
            V += c * (V_max - V) * self.delta_t  # Euler-update

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def logistisch_model(self, c, V_max):
        """
        Simuleert logistische groei waarbij de groeisnelheid afhankelijk is van
        zowel het huidige volume V als de resterende capaciteit (V_max - V).
        Differentialvergelijking:
        dV/dt = c * V * (V_max - V)

        Parameters:
            c : groeiconstante
            V_max : maximale waarde waarheen het model groeit

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t  # update tijd
            V += c * V * (V_max - V) * self.delta_t  # Euler-update

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def montroll_model(self, c, d, V_max):
        """
        Simuleert Montroll-groei waarbij de groeisnelheid afhankelijk is van
        V en afneemt naarmate het volume dichter bij een maximale waarde V_max komt,
        volgens de niet-lineaire exponent d.
        Differentialvergelijking:
        dV/dt = c * V * (V_max^d - V^d)

        Parameters:
            c : groeiconstante
            d : exponent van V, bepaalt niet-lineariteit van de groei
            V_max : maximale waarde waarheen het model groeit

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]            # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume     # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t # update tijd
            V += c * V * (V_max**d - V**d) * self.delta_t  # Euler-update

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def allee_model(self, c, V_min, V_max):
        """
        Simuleert Allee-effect groei waarbij de groeisnelheid afhangt van
        het verschil met een minimumwaarde V_min en een maximumwaarde V_max.
        Differentialvergelijking:
        dV/dt = c * (V - V_min) * (V_max - V)

        Parameters:
            c : groeiconstante
            V_min : minimumwaarde nodig voor groei (Allee-drempel)
            V_max : maximale waarde waarheen het model groeit

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t  # update tijd
            V += c * (V - V_min) * (V_max - V) * self.delta_t  # Euler-update

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def lineair_gelimiteerd_model(self, c, d):
        """
        Simuleert lineair gelimiteerde groei waarbij de groeisnelheid afneemt
        naarmate V toeneemt door een limiterende factor d.
        Differentialvergelijking:
        dV/dt = c * V / (V + d)

        Parameters:
            c : groeiconstante
            d : limiterende constante die afvlakking bepaalt

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t
            V += c * V / (V + d) * self.delta_t  # Euler-update

            Ts.append(t)  # sla het huidige tijdstip op
            Vs.append(V)  # sla het huidige volume op

        return Ts, Vs

    def oppervlakte_gelimiteerd_model(self, c, d):
        """
        Simuleert oppervlakte-gelimiteerde groei waarbij de groeisnelheid
        afneemt door een wortelfactor in het volume.
        Differentialvergelijking:
        dV/dt = c * V / (V + d)^(1/3)

        Parameters:
            c : groeiconstante
            d : limiterende constante

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t
            V += c * V / (V + d) ** (1 / 3) * self.delta_t  # Euler-update

            Ts.append(t)
            Vs.append(V)

        return Ts, Vs

    def von_bertalanffy_model(self, c, d):
        """
        Simuleert Von Bertalanffy-groei waarbij de aanmaak afhangt van V^(2/3)
        en het verlies proportioneel is aan V.

        Differentialvergelijking:
        dV/dt = c * V^(2/3) - d * V

        Parameters:
            c : groeiconstante voor aanmaak
            d : verliesconstante

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t
            V += (c * V ** (2 / 3) - d * V) * self.delta_t  # Euler-update

            Ts.append(t)
            Vs.append(V)

        return Ts, Vs

    def gompertz_model(self, c, V_max):
        """
        Simuleert Gompertz-groei waarbij de groeisnelheid exponentieel afneemt
        naarmate het volume dichter bij V_max komt.
        Differentialvergelijking:
        dV/dt = c * V * ln(V_max / V)

        Parameters:
            c : groeiconstante
            V_max : maximale waarde waarheen het model groeit

        Returns:
            Ts : lijst van tijdstappen
            Vs : lijst van volume bij elke tijdsstap beginnend bij self.volume
        """
        Ts = [0]  # beginpunt t=0
        Vs = [self.volume]  # beginvolume
        t = 0
        V = self.volume  # lokale kopie van volume

        # Loop over N tijdstappen
        for _ in range(self.n):
            t += self.delta_t
            V += c * V * (np.log(V_max / V)) * self.delta_t  # Euler-update

            Ts.append(t)
            Vs.append(V)

        return Ts, Vs








