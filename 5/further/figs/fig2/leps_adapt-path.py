import numpy as np
import matplotlib.pyplot as plt
plt.style.use('paper')

a = 0.05
b = 0.8
c = 0.0
r0 = 0.742
alpha = 1.942


def q(r, d):
    return 0.5 * d * (1.5 * np.exp(-2*alpha*(r - r0))
                      - np.exp(-alpha*(r - r0)))


def j(r, d):
    return 0.25 * d * (np.exp(-2*alpha*(r - r0))
                       - 6 * np.exp(-alpha*(r - r0)))


def leps(r_ab, r_bc):
    """x := r_AB    y := r_BC"""

    r_ac = r_ab + r_bc

    j_ab = j(r_ab, d=4.746)
    j_bc = j(r_bc, d=4.746)
    j_ac = j(r_ac, d=3.445)

    v = (q(r_ab, d=4.746) / (1 + a) +
         q(r_bc, d=4.746) / (1 + b) +
         q(r_ac, d=3.445) / (1 + c) -
         (j_ab**2 / (1 + a)**2 +
          j_bc**2 / (1 + b)**2 +
          j_ac**2 / (1 + c)**2 -
          j_ab*j_bc / ((1 + a) * (1+b)) -
          j_bc*j_ac / ((1 + b) * (1 + c)) -
          j_ab*j_ac / ((1 + a) * (1 + c)))
         )

    return v


def leps_alt(_x, _y, x0=1.93):
    return leps(_x, _y) - 2 * np.artan(5*(x-x0)) - 2*x


def plot_adaptive_path(h=1E-6, g=50):
    dr_min, dr_max = 0.05, 0.3
    dr_m = dr_max - dr_min
    dr_init = 0.05

    _x, _y = x0 + dr_init, y0 - dr_init
    x_path, y_path = [x0, _x], [y0, _y]

    while _x != xn or _y != yn:

        grad_x = (leps(_x + h, _y) - leps(_x, _y)) / h
        grad_y = (leps(_x, _y + h) - leps(_x, _y)) / h

        if grad_x < 0:
            _x = min(_x + dr_max, xn)
        else:
            _x = min(_x + dr_m*np.exp(-(grad_x / g)**2) + dr_min,
                     xn)

        if grad_y < 0:
            _y = max(_y - dr_max, yn)
        else:
            _y = max(_y - (dr_m*np.exp(-(grad_y / g)**2) + dr_min),
                     yn)

        x_path.append(_x)
        y_path.append(_y)

    plt.plot(x_path, y_path, marker='o', c='white', lw=1, ms=5)

    return None


if __name__ == '__main__':

    x, y = np.meshgrid(np.linspace(0, 2, 20),
                       np.linspace(0, 2, 20))

    x0, y0 = 0.3, 2.0
    xn, yn = 2.0, 0.3

    plt.pcolor(x, y, leps(x, y), vmax=65)
    plot_adaptive_path()
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    plt.tight_layout()
    plt.savefig('tmp.png', dpi=200)
