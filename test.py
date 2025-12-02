import inspect
import copy
import numpy as np
import math

class tumor_growth_models_ode:
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
            volume: Startvolume
            delta_t: Tijdstap Δt
            n: Aantal tijdstappen
        """
        self.volume = volume
        self.delta_t = delta_t
        self.n = n

    # Numerieke integratiemethoden
    def _step_euler(self, f, V, t, dt):
        """Euler integratie stap."""
        return V + f(V, t) * dt

    def _step_heun(self, f, V, t, dt):
        """Heun integratie stap."""
        k1 = f(V, t)
        k2 = f(V + k1*dt, t + dt)
        return V + 0.5 * (k1 + k2) * dt

    def _step_rk4(self, f, V, t, dt):
        """Runge-Kutta 4 integratie stap (4e orde)."""
        k1 = f(V, t)
        k2 = f(V + 0.5*k1*dt, t + 0.5*dt)
        k3 = f(V + 0.5*k2*dt, t + 0.5*dt)
        k4 = f(V + k3*dt, t + dt)
        return V + (k1 + 2*k2 + 2*k3 + k4)/6 * dt


    # Simulatie-functie
    def _simulate(self, f, methode = "euler"):
        """
        Simuleer een ODE-model met de opgegeven integratiemethode, met euler als default.

        Parameters:
            f (callable): functie f(V, t) die dV/dt retourneert
            methode (str): 'euler', 'heun', of 'rk4' met euler als default

        Returns:
            Ts (list[float]): tijdstappen
            Vs (list[float]): volumes bij elke tijdstap
        """
        # Kies de juiste integratie-stap functie
        stepper = {
            "euler": self._step_euler,
            "heun": self._step_heun,
            "rk4": self._step_rk4
        }[methode.lower()]

        Ts = [0]
        Vs = [self.volume]
        V = self.volume
        t = 0

        for _ in range(self.n):
            V = stepper(f, V, t, self.delta_t)  # update volume
            t += self.delta_t
            Ts.append(t)
            Vs.append(V)

        return Ts, Vs

    # groeimodellen
    def lineaire_model(self, c, methode = "euler"):
        """Lineaire groei: dV/dt = c"""
        return self._simulate(lambda V, t: c, methode)

    def exponentieel_model(self, c, methode = "euler"):
        """Exponentiële groei: dV/dt = c * V"""
        return self._simulate(lambda V, t: c * V, methode)

    def mendelsohn_model(self, c, d, methode = "euler"):
        """Mendelsohn model: dV/dt = c * V^d"""
        return self._simulate(lambda V, t: c * V**d, methode)

    def exponentieel_afvlakkend_model(self, c, V_max, methode = "euler"):
        """Exponentiële afvlakkende groei: dV/dt = c * (V_max - V)"""
        return self._simulate(lambda V, t: c * (V_max - V), methode)

    def logistisch_model(self, c, V_max, methode = "euler"):
        """Logistische groei: dV/dt = c * V * (V_max - V)"""
        return self._simulate(lambda V, t: c * V * (V_max - V), methode)

    def montroll_model(self, c, d, V_max, methode = "euler"):
        """Montroll model: dV/dt = c * V * (V_max^d - V^d)"""
        return self._simulate(lambda V, t: c * V * (V_max**d - V**d), methode)

    def allee_model(self, c, V_min, V_max, methode = "euler"):
        """Allee-effect model: dV/dt = c * (V - V_min) * (V_max - V)"""
        return self._simulate(lambda V, t: c * (V - V_min) * (V_max - V), methode)

    def lineair_gelimiteerd_model(self, c, d, methode = "euler"):
        """Lineair gelimiteerd model: dV/dt = c * V / (V + d)"""
        return self._simulate(lambda V, t: c * V / (V + d), methode)

    def oppervlakte_gelimiteerd_model(self, c, d, methode = "euler"):
        """Oppervlakte-gelimiteerd model: dV/dt = c * V / (V + d)^(1/3)"""
        return self._simulate(lambda V, t: c * V / (V + d)**(1/3), methode)

    def von_bertalanffy_model(self, c, d, methode = "euler"):
        """Von Bertalanffy model: dV/dt = c * V^(2/3) - d * V"""
        return self._simulate(lambda V, t: c * V**(2/3) - d * V, methode)

    def gompertz_model(self, c, V_max, methode= "euler"):
        """Gompertz model: dV/dt = c * V * ln(V_max / V)"""
        return self._simulate(lambda V, t: c * V * math.log(V_max / V), methode)


    def MSE(self, model_func, methode, params, data_ts, data_vs):
        """
        Bereken de Mean Squared Error (MSE) tussen model en experimentele data.

        Parameters:
            model_func : functie van het model (bijv. self.exponentieel_model)
            methode    : integratiemethode ('euler', 'heun', 'rk4')
            params     : dict met parameters voor het model
            data_ts    : lijst/array van tijdstippen van data
            data_vs    : lijst/array van gemeten volumes

        Returns:
            mean_squared_error
        """
        # Filter parameters die model accepteert
        sig = inspect.signature(model_func)
        filtered_params = {k: v for k, v in params.items() if k in sig.parameters}

        # Simuleer model
        model_ts, model_vs = model_func(methode=methode, **filtered_params)

        # Interpoleer model op tijdstippen van data
        model_interp = np.interp(data_ts, model_ts, model_vs)

        # Bereken MSE
        errors = np.array(data_vs) - model_interp
        return np.mean(errors ** 2)


    def hooke_jeeves(self, model_func, params, data_ts, data_vs, methode="euler",
                     tol=1e-6, alpha_up=1.2, alpha_down=0.2, max_iter=10000):
        """
        Hooke & Jeeves / Direct Search optimalisatie voor ODE modellen.
        Zoekt parameters die MSE minimaliseren.

        Returns:
            filtered_params : dict van geoptimaliseerde parameters
            mse             : MSE van het optimale model
        """
        sig = inspect.signature(model_func)
        valid_keys = sig.parameters.keys()
        deltas = {k: 0.1 * max(1.0, abs(v)) for k, v in params.items()}
        mse = self.MSE(model_func, methode, params, data_ts, data_vs)
        iteration = 0

        while max(abs(d) for d in deltas.values()) > tol and iteration < max_iter:
            iteration += 1
            for key in params:
                improved = False
                new_params = copy.deepcopy(params)

                # omhoog proberen
                new_params[key] = params[key] + deltas[key]
                new_mse = self.MSE(model_func, methode, new_params, data_ts, data_vs)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= alpha_up
                    improved = True
                    continue

                # omlaag proberen
                new_params[key] = params[key] - deltas[key]
                new_mse = self.MSE(model_func, methode, new_params, data_ts, data_vs)
                if new_mse < mse:
                    params = new_params
                    mse = new_mse
                    deltas[key] *= -alpha_up
                    improved = True
                    continue

                # geen verbetering → stap verkleinen
                if not improved:
                    deltas[key] *= alpha_down

        filtered_params = {k: v for k, v in params.items() if k in valid_keys}
        return filtered_params, mse

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
        if mse <= 0:
            return np.inf, np.inf, np.inf

        aic = n_data * np.log(mse) + 2 * n_params
        bic = n_data * np.log(mse) + n_params * np.log(n_data)
        if n_data / n_params < 40:  # kleine dataset -> AICc
            aicc = aic + (2 * n_params * (n_params + 1)) / (n_data - n_params - 1)
        else:
            aicc = aic
        return aic, aicc, bic

    def fit_and_evaluate(
        self,
        model_func,
        params,
        data_ts,
        data_vs,
        methode="rk4"
    ):
        """
        Fit een model op data en retourneer MSE, AIC, AICc, BIC.

        Parameters:
            model_func : functie van het model
            params     : dict met startwaarden van parameters
            data_ts    : lijst/array van tijdstippen
            data_vs    : lijst/array van gemeten volumes
            methode    : integratiemethode, standaard rk4 Runge-Kutta 4e orde, nauwkeurig numerieke integratiemethode
            (kan ook euler of heun zijn)

        Returns:
            dict met geoptimaliseerde parameters en statistieken
        """
        best_params, mse = self.hooke_jeeves(model_func, params, data_ts, data_vs, methode=methode)
        n_data = len(data_vs)
        n_params = len(best_params)
        aic, aicc, bic = self.information_criteria(mse, n_data, n_params)

        return {
            "params": best_params,
            "mse": mse,
            "AIC": aic,
            "AICc": aicc,
            "BIC": bic
        }
