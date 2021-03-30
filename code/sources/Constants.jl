"""
    a

Scale length of the Kuzmin-Toomre disk.
Set to `1.0` by default.
"""
const a = 1.0

"""
    b

Scale length of the Toomre spherical dark halo.
Set to `2.8` by default, see `Reddish et al. (2021)`.
"""
const b = 2.8

"""
    c

Scale length of the Toomre spherical bulge.
Set to `1/20` by default, see `Reddish et al. (2021)`.
"""
const c = 1/20

"""
    error_dis

Error threshold on the convergence distance when determining
physical eigenvalues of the system.
"""
const error_dis = 0.005

"""
    error_var

Error threshold on the convergence variance when determining
physical eigenvalues of the system.
"""
const error_var = 0.005

"""
    threshold_nb_egv

Minimal number of cutoff matrix eigenvalues that must be found
to apply the selection criterion.
"""
const threshold_nb_egv = 10

"""
    threshold_GR

Lower threshold under which we don't compute the growth rate of the system.
Set to `0.1` by default, see `Reddish et al. (2021)`.
"""
const threshold_GR = 0.1
