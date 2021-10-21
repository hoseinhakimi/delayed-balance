import _thread
import time
import delay as delay
import matplotlib.pyplot as plt

nodes_count = 20
time_steps = nodes_count ** 3
max_pow = 12
delays = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
# results = []

print("sd")
def lifeTime(start):
  time.sleep(start*5)
  for i in range(int(start * max_pow/4), int((start + 1) * max_pow/4)):
    average_life = 0
    T = 2 ** i
    for _ in range(10):
      # trngl_energy_ts = []
      dl = delay.delay()
      dl.set_system(nodes_count, 0.5, T, 1)
      for t in range(time_steps):
        dl.flipper(t)
        e = dl.tri_energy / dl.tc
        # trngl_energy_ts.append(e)
        if (e == -1):
          average_life += t
          break
    print(average_life/10)
    # results.append([T, average_life/10])


# try:
# lifeTime(1)
for i in range(4):
  _thread.start_new_thread( lifeTime, (i,))
# except:
#    print("Error: unable to start thread")


# print(results)
