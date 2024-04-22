import numpy as np
import matplotlib.pyplot as plt

z1 = 37
z2 = 95

beta = np.linspace(0, 45*np.pi/180, 20)

m_n = np.linspace(8, 11, 7)

plt.figure()
for i, m_n in enumerate(m_n):
    # print(m_n)
    mt = m_n / np.cos(beta)
    a = mt * (z1 + z2) / 2
    # print(a[0])
    plt.plot(beta*180/np.pi, a, label=f"m_n = {m_n}")
plt.axhline(y=725, color='r', linestyle='-', label="y=725")



# import tikzplotlib
# tikzplotlib.save("HandIn/Sized.tex")

plt.legend()
plt.show()
