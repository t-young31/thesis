import numpy as np

m = 18 * 1.66053906660E-27       # 1 amu in kg
kb = 1.380649E-23           # J K-1
T = 298                       # K
h = 6.62607004E-34            # J s
beta = 1.0 / (kb * T)

max_n = 250
nx, ny, nz = np.meshgrid(np.arange(0, max_n),
                         np.arange(0, max_n),
                         np.arange(0, max_n),
                         indexing='ij')

for l_ang in np.linspace(10, 1, num=10):
    l = l_ang*10**(-10)

    q_sum = np.sum(np.exp(-(beta*h**2 /
                            (8.0*m*l**2))
                          * (nx**2 + ny**2 + nz**2)))

    q_integral = ((2.0*np.pi*m*kb*T)/(h**2))**1.5*l**3

    print(100*(np.log(q_sum) - np.log(q_integral))/np.log(q_sum),
          f'l = {l_ang} Å')
