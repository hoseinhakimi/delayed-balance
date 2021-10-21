import delayed_balance as blc
import random
import matplotlib.pyplot as plt
size = 30
delays = [0,1, 20000]
for i in delays:
    nblc = blc.TriangularBalance(size, -1, i)
    nblc.TriadDynamics(3)
    # plt.plot(nblc.TimeLine[:, 1], nblc.TimeLine[:, 0], '.')
    plt.plot(nblc.TimeLine[:, 1], '.', label=str(i),markersize=1)
    # plt.savefig("s"+str(i))
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.savefig("final4")
# plt.show()
