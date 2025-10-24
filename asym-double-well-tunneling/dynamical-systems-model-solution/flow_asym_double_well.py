from __future__ import print_function, division
import numpy as np
"""
This class implements the four-dimensional dynamical system
governing the time evolution of the mean and variance as derived
in the manuscript.

Last updated: Oct 24th, 2025
"""

class ModelFlow:
    def __init__(self, a, b, c, energy, skewness):
        """Initializes the model with given energy and mass values."""
        self.a = a
        self.b = b
        self.c = c
        self.E = energy
        self.S = skewness
        self.dims = 4  # Number of dimensions in the system

    def flow(self, time, state):
        """Defines the flow of the system for solve_ivp."""
        if not isinstance(state, (list, np.ndarray)):
            raise TypeError(f"Expected state to be iterable, got {type(state)}")
        x, p, v, w = state
        dx = p
        dp = -self.a*x + self.b*(v+x**2) -self.c*(self.S + 3.0*v*x + x**3)
        dv = w
        dw = 4.0*self.E - 2.0*p**2 - self.a*(4.0*v + 2.0*x**2) + self.b*(10.0/3.0*self.S + 8.0*v*x + 4.0/3.0*x**3) -self.c*(9.0*v**2 + 10.0*self.S*x + 12.0*v*x**2 + x**4)
        return [dx, dp, dv, dw]
