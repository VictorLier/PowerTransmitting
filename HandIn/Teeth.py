import numpy as np
import matplotlib.pyplot as plt


i_goal = 292 / 750

z1 = np.arange(20,151,1)
z2 = np.rint(z1 / i_goal)

# Trim z2 vector such that the highest number is under 150
z2 = z2[z2 <= 150]
# Trim z1 vector to match the length of z2
z1 = z1[:len(z2)]

actuali = z1 / z2


alpha_w = 725
beta = 15* np.pi / 180
possible_beta = np.arange(5, 20+5, 5)*np.pi/180
possible_module = np.array([6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18])
alpha_n = 20 * np.pi / 180

color = ['r', 'g', 'b', 'y', 'c', 'm', 'k', 'w']

plt.figure()
for j, beta in enumerate(possible_beta):
    profile_shift = np.zeros((len(possible_module), len(z1)))
    for i, m_t in enumerate(possible_module):
        for n in range(len(z1)):
            a = m_t * (z1[n] + z2[n]) / 2
            alpha_t = np.arctan(np.tan(alpha_n) / np.cos(beta))
            invalpha_t = np.tan(alpha_t) - alpha_t
            alpha_wt = np.arccos(a * np.cos(alpha_t) / alpha_w)
            invalpha_wt = np.tan(alpha_wt) - alpha_wt
            profile_shift[i, n] = (z1[n] + z2[n]) / 2 * (invalpha_wt - invalpha_t) / np.tan(alpha_t)
    profile_shift[profile_shift < 0] = np.nan
    profile_shift[profile_shift > 1.5] = np.nan

    for n in range(len(profile_shift[0,:])):
        for i in range(len(profile_shift[:,0])):
            if not np.isnan(profile_shift[i,n]):
                print(z1[n], possible_module[i])
                plt.scatter(z1[n], possible_module[i])
                # print(z1[n], profile_shift[i,n])
                plt.scatter(z1[n], profile_shift[i,n], color=color[j])
    print("\n")
plt.show()



print("Stop")