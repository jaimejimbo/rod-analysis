import methods

file_ = open('num_rods.data', 'r')
times = []
data = []
i = 0
for line in file_:
    i += 1
    times.append(i)
    data_ = line.split(' ')
    data.append(int(data_[2]))

import matplotlib.pyplot as plt
plt.figure()
plt.plot(times, data)
plt.savefig('num_rods.png')
