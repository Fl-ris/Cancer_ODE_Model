import inspect
import copy
import numpy as np
import math
import matplotlib.pyplot as plt

class tumorODE:
    """
    Klasse voor simulatie van tumor-groei met verschillende ODE-modellen.

    Ondersteunde modellen:
        - lineaire_model
        - exponentieel_model
        - mendelsohn_model
        - exponentieel_afvlakkend_model
        - logistisch_model
        - montroll_model
        - allee_model
        - lineair_gelimiteerd_model
        - oppervlakte_gelimiteerd_model
        - von_bertalanffy_model
        - gompertz_model

    Integratiemethoden:
        - Euler
        - Heun
        - Runge-Kutta 4 (RK4)
    """

    def __init__(self, volume: float, delta_t: float, n: int):
        """
        Initialiseer het groeimodel.

        Parameters:
            volume: Startvolume (V0)
            delta_t: Tijdstap grootte (dt)
            n: Aantal tijdstappen om te simuleren
        """
        self.start_volume = volume
        self.delta_t = delta_t
        self.n = n

    def _step_euler(self, f, V, t, dt):
        """Euler integratie stap."""
        return V + f(V, t) * dt

    def _step_heun(self, f, V, t, dt):
        """Heun methode (Predictor-Corrector)."""
        k1 = f(V, t)
        k2 = f(V + k1*dt, t + dt)
        return V + 0.5 * (k1 + k2) * dt

    def _step_rk4(self, f, V, t, dt):
        """
        Runge-Kutta 4e Orde.
        """
        k1 = f(V, t)
        k2 = f(V + 0.5*k1*dt, t + 0.5*dt)
        k3 = f(V + 0.5*k2*dt, t + 0.5*dt)
        k4 = f(V + k3*dt, t + dt)
        return V + (k1 + 2*k2 + 2*k3 + k4)/6 * dt

    def _simulate(self, f, methode="rk4"):
        """
        Simuleer een ODE-model met de opgegeven integratiemethode, met euler als default.

        Parameters:
            f (callable): functie f(V, t) die dV/dt retourneert
            methode (str): 'euler', 'heun', of 'rk4' met euler als default

        Returns:
            Ts (list[float]): tijdstappen
            Vs (list[float]): volumes bij elke tijdstap
        """
        stepper = {
            "euler": self._step_euler,
            "heun": self._step_heun,
            "rk4": self._step_rk4
        }.get(methode.lower(), self._step_rk4)

        Ts = [0]
        Vs = [self.start_volume]
        V = self.start_volume
        t = 0

        for _ in range(self.n):
            V = stepper(f, V, t, self.delta_t)
            t += self.delta_t
            Ts.append(t)
            Vs.append(V)

        return Ts, Vs


    def lineaire_model(self, c, methode="rk4"):
        """Dv/Dt = c"""
        return self._simulate(lambda V, t: c, methode)

    def exponentieel_model(self, c, methode="rk4"):
        """Dv/Dt = c * V"""
        return self._simulate(lambda V, t: c * V, methode)

    def mendelsohn_model(self, c, d, methode="rk4"):
        """Dv/Dt = c * V^d"""
        # Anders math domain error..
        return self._simulate(lambda V, t: c * math.pow(max(1e-9, V), d), methode)

    def logistisch_model(self, c, V_max, methode="rk4"):
        """Dv/Dt = c * V * (1 - V/Vmax)"""
        return self._simulate(lambda V, t: c * V * (1 - V/V_max), methode)

    def gompertz_model(self, c, V_max, methode="rk4"):
        """Dv/Dt = c * V * ln(Vmax / V)"""
        # Mag geen log(0) zijn...
        func = lambda V, t: c * V * math.log(V_max / V) if V > 1e-9 else 0
        return self._simulate(func, methode)

    def von_bertalanffy_model(self, c, d, methode="rk4"):
        """Dv/Dt = c * V^(2/3) - d * V"""
        return self._simulate(lambda V, t: c * math.pow(max(0, V), 2/3) - d * V, methode)

    def exponentieel_afvlakkend_model(self, c, V_max, methode="rk4"):
        """Dv/Dt = c * (Vmax - V)"""
        return self._simulate(lambda V, t: c * (V_max - V), methode)

    def allee_effect_model(self, c, V_min, V_max, methode="rk4"):
        """Dv/Dt = c * (V - Vmin) * (Vmax - V)"""
        return self._simulate(lambda V, t: c * (V - V_min) * (V_max - V), methode)

    def lineair_gelimiteerd_model(self, c, d, methode="rk4"):
        """Dv/Dt = c * V / (V + d)"""
        return self._simulate(lambda V, t: c * (V / (V + d)), methode)

    def oppervlak_gelimiteerd_model(self, c, d, methode="rk4"):
        """Dv/Dt = c * V / (V + d)^(1/3)"""
        return self._simulate(lambda V, t: c * V / math.pow((V + d), 1/3), methode)


    def MSE(self, model_func, methode, params, data_ts, data_vs):
        """Bereken de Mean Squared Error tussen model en experimentele data."""
        # Filter parameters zodat alleen de benodigde params naar het model gaan
        sig = inspect.signature(model_func)
        gefilterde_params = {k: v for k, v in params.items() if k in sig.parameters}

        # Simuleer model
        model_ts, model_vs = model_func(methode=methode, **gefilterde_params)

        # Interpoleer model resultaten op exact dezelfde tijdstippen als de data
        model_interp = np.interp(data_ts, model_ts, model_vs)

        # Bereken Error
        errors = np.array(data_vs) - model_interp
        return np.mean(errors ** 2)

    def hooke_jeeves(self, model_func, params, data_ts, data_vs, methode="rk4",
                     tol=1e-6, alpha_up=1.2, alpha_down=0.5, max_iter=10000):
        """
        Hooke & Jeeves / Direct Search optimalisatie.
        Zoekt parameters die de MSE minimaliseren.
        """
        sig = inspect.signature(model_func)
        valid_keys = [k for k in sig.parameters.keys() if k != 'methode']
        
        # Stapgrootte initialisatie
        deltas = {k: 0.1 * max(1.0, abs(v)) for k, v in params.items() if k in valid_keys}
        
        huidige_mse = self.MSE(model_func, methode, params, data_ts, data_vs)
        iteratie = 0

        while max(abs(d) for d in deltas.values()) > tol and iteratie < max_iter:
            iteratie += 1
            for key in valid_keys:
                if key not in params: continue
                
                verbeterd = False
                beste_waarde_in_stap = params[key]
                
                # Probeer parameter te verhogen
                params[key] += deltas[key]
                up_mse = self.MSE(model_func, methode, params, data_ts, data_vs)
                
                if up_mse < huidige_mse:
                    huidige_mse = up_mse
                    beste_waarde_in_stap = params[key]
                    deltas[key] *= alpha_up # Versnel
                    verbeterd = True
                else:
                    # Probeer parameter te verlagen (terug naar origineel - delta)
                    params[key] -= 2 * deltas[key] 
                    down_mse = self.MSE(model_func, methode, params, data_ts, data_vs)
                    
                    if down_mse < huidige_mse:
                        huidige_mse = down_mse
                        beste_waarde_in_stap = params[key]
                        deltas[key] *= alpha_up # Versnel
                        verbeterd = True
                    else:
                        # Geen verbetering, reset naar origineel
                        params[key] += deltas[key]

                params[key] = beste_waarde_in_stap
                
                # Als er geen verbetering is gevonden, verklein de stapgrootte
                if not verbeterd:
                    deltas[key] *= alpha_down

        # Return alleen relevante params
        eind_params = {k: v for k, v in params.items() if k in valid_keys}
        return eind_params, huidige_mse

    def informatie_criteria(self, mse, n_data, n_params):
        """
        Bereken AIC, AICc (voor kleine datasets) en BIC.

        Parameters:
            mse      : mean squared error van model
            n_data   : aantal datapunten
            n_params : aantal parameters in model

        Returns:
            aic, aicc, bic
        """

        aic = n_data * np.log(mse) + 2 * n_params
        bic = n_data * np.log(mse) + n_params * np.log(n_data)     
        aicc = aic + (2 * n_params * (n_params + 1)) / (n_data - n_params - 1)

        return aic, aicc, bic

    def fit_and_evaluate(self, model_func, start_params, data_ts, data_vs, methode="rk4"):
        """
        Fit een model op data en retourneer MSE, AIC, en optimale parameters.
        """
        best_params, mse = self.hooke_jeeves(model_func, start_params, data_ts, data_vs, methode=methode)
        n_data = len(data_vs)
        n_params = len(best_params)
        
        aic, aicc, bic = self.informatie_criteria(mse, n_data, n_params)

        return {
            "model_naam": model_func.__name__,
            "functie": model_func,
            "best_params": best_params,
            "mse": mse,
            "AIC": aic,
            "AICc": aicc,
            "BIC": bic
        }


    def plot(self, Ts, Vs, color=None, label=None):
        """Plot een enkele simulatielijn."""

        plt.plot(Ts, Vs, color=color, label=label)


    def show_plot(self, titel="Tumor Groei Simulatie"):
        """Maak de grafiek op en toon deze."""

       
        plt.axvline(0.0, color='k', linestyle='--', alpha=0.5)
        plt.xlabel("Tijd")
        plt.ylabel("Volume")
        plt.title(titel)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()


if __name__ == "__main__":
    

    # Kleine testdataset van: [S.S. Hassan & H.M. Al-Saedi, 2024](https://doi.org/10.1051/bioconf/20249700118)
    Ts = [0, 13, 20, 32, 42, 55, 65, 75, 85, 88, 95, 98, 107, 115, 120]
    Vs = [250, 255, 550, 575, 576, 800, 1050, 1250, 1750, 2000, 2550, 2750, 3000, 3500, 4000]

    modeler = tumorODE(volume=Ts[0], delta_t=1, n=120)

    modellen_lijst = [
        (modeler.lineaire_model, {'c': 10}),
        (modeler.exponentieel_model, {'c': 0.1}),
        (modeler.mendelsohn_model, {'c': 0, 'd': 0}),
        (modeler.logistisch_model, {'c': 0.2, 'V_max': 800}),
        (modeler.gompertz_model, {'c': 0.2, 'V_max': 800}),
        (modeler.von_bertalanffy_model, {'c': 0.5, 'd': 0.1}),
        (modeler.exponentieel_afvlakkend_model, {'c': 0.1, 'V_max': 800}),
        (modeler.lineair_gelimiteerd_model, {'c': 100, 'd': 500}),
    ]

    print(f"Model ,MSE ,AIC")

    resultaten = []

    for model_func, start_params in modellen_lijst:
            res = modeler.fit_and_evaluate(
                model_func, 
                start_params, 
                Ts, 
                Vs
            )
            resultaten.append(res)
            print(f"{res['model_naam']:<35}, {res['mse']:.2f}, {res['AIC']:.2f}")

    resultaten.sort(key=lambda x: x['AIC'])
    beste_resultaat = resultaten[0]


    print(f"Beste model: {beste_resultaat['model_naam']}")
    print(f"Optimale parameters: {beste_resultaat['best_params']}")
    print(f"Beste model: {beste_resultaat['model_naam']}")

    plt.scatter(Ts, Vs, color="black", label="Werkelijke tumordata")


    beste_func = beste_resultaat['functie']
    beste_params = beste_resultaat['best_params']
    ts, vs = beste_func(**beste_params)
    
    modeler.plot(ts, vs, color="red", label=f"Fit: {beste_resultaat['model_naam']}")

    modeler.show_plot()