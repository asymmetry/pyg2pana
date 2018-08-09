"""
Python Wrapper for Peter Bosted's Fits
======================================

Calculate inclusive inelastic electron scattering cross sections for different
nuclei, using Peter Bosted's fits to structure functions F1 and F2.

See Peter Bosted's page https://userweb.jlab.org/~bosted/ for more references.

    xs -- Calculate non-radiated cross section
    xs_rad -- Calculate radiated cross section
"""

from .pbosted import xs, xs_rad

__all__ = [s for s in dir() if not s.startswith('_')]
