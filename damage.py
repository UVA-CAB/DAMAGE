# -------------------------------------------------------------------------
# Version: 2.0 (Released: April 2025)
# -------------------------------------------------------------------------
# This script calculates the Diffuse Axonal, Multi-Axis, General Evaluation
# (DAMAGE) metric: Gabler, L.F., Crandall, J.R. & Panzer, M.B. Ann Biomed
# Eng, 2018. https://doi.org/10.1007/s10439-018-02179-9.
# -------------------------------------------------------------------------
# Administrative contact: Ramakrishnan Iyer (zdf5ks@virginia.edu)
# Research contact: Matthew B. Panzer (panzer@virginia.edu)
# -------------------------------------------------------------------------
# Input variable:
# X: Timed head kinematics (nx4), where n is the # of time points.
#    Column 1: Time (s)
#    Column 2: X angular acceleration (rad/s/s)
#    Column 3: Y angular acceleration (rad/s/s)
#    Column 4: Z angular acceleration (rad/s/s)
#    Data entered must be measured with respect to the local head
#    coordinate system defined by SAE J211. Do not enter column headers.
#
# Output variable:
# DAMAGE value
# -------------------------------------------------------------------------


import numpy as np
import pandas as pd


def calculate_damage(T, X):
    """
    Newmark-beta time integration for MDOF system.

    Parameters
    ----------
    T : (n,) array_like
        Time vector
    F : (n, m) array_like
        head rotation acceleration matrix (each row = time step, columns = DOFs)

    Returns
    -------
    x : (n, m) ndarray
        Displacements
    dx : (n, m) ndarray
        Velocities
    d2x : (n, m) ndarray
        Accelerations
    """

    T = np.asarray(T)
    F = np.column_stack((X[:, 0:2], -X[:, 2]))
    F = np.asarray(F)

    n = len(T)
    m = F.shape[1]

    # Initialize response arrays
    x = np.zeros((n, m))    # displacement
    dx = np.zeros((n, m))   # velocity
    d2x = np.zeros((n, m))  # acceleration

    # Newmark-beta parameters (average acceleration method)
    beta1 = 0.5
    beta2 = 0.5

    # System matrices
    M = np.eye(3)

    kxx, kyy, kzz = 32142, 23493, 16935
    kxy, kyz, kxz = 0, 0, 1636.3

    K = np.array([
        [kxx + kxy + kxz, -kxy, -kxz],
        [-kxy, kxy + kyy + kyz, -kyz],
        [-kxz, -kyz, kxz + kyz + kzz]
    ])

    a1 = 5.9148e-3
    C = a1 * K

    # Initial step
    dT = T[1] - T[0]
    Meff = M + dT * beta1 * C + 0.5 * dT**2 * beta2 * K

    d2x[0, :] = np.linalg.solve(Meff, F[0, :])

    # Time stepping
    for i in range(1, n-1):
        dT = T[i] - T[i - 1]

        Meff = M + dT * beta1 * C + 0.5 * dT**2 * beta2 * K
        # print(i)

        # Effective force
        Feff = (
            F[i, :]
            - C @ (dx[i - 1, :] + dT * (1 - beta1) * d2x[i - 1, :])
            - K @ (
                x[i - 1, :]
                + dT * dx[i - 1, :]
                + 0.5 * dT**2 * (1 - beta2) * d2x[i - 1, :]
            )
        )

        # Solve for acceleration
        d2x[i, :] = np.linalg.solve(Meff, Feff)

        # Update velocity and displacement
        d2x1 = (1 - beta1) * d2x[i - 1, :] + beta1 * d2x[i, :]
        d2x2 = (1 - beta2) * d2x[i - 1, :] + beta2 * d2x[i, :]

        dx[i, :] = dx[i - 1, :] + dT * d2x1
        x[i, :] = x[i - 1, :] + dT * dx[i - 1, :] + 0.5 * dT**2 * d2x2
        B = 2.9903 #% (1/m) scale factor
        dmg_i_newmark = B*x
        dmg_t_newmark = np.maximum.accumulate(
            np.sqrt(np.sum(dmg_i_newmark ** 2, axis=1))
        )

        dmg_newmark = max(dmg_t_newmark)

    return dmg_newmark
if __name__ == '__main__':
    # njit(float64[::1](float64[::1], float64, float64[::1, :], float64[::1, :], float64[::1, :], float64[::1, :],
    #                   float64[::1], float64[::1], float64[::1], float64[::1]))(equation_of_motion).inspect_types()

    def main():
        df = pd.read_excel("example_data.xlsx", sheet_name="Sheet1")
        X = df.to_numpy()
        T = X[:, 0]
        Ax= X[:,1]
        Ay = X[:, 2]
        Az = X[:, 3]
        damage = calculate_damage( T,  np.stack((Ax, Ay, Az), axis = 1))
        damage = round (damage, 4)
        print ('DAMAGE = ', damage)

    main()
