import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import statistics as stats

myLables = ['256', '1024', '4096', '16384', '65536', '262144', '1048576']

# convert file double into List array double 
def makeArray(filename):
    temp = np.empty((7, 0)).tolist()
    index = 0
    with open(filename, 'r') as data:
        datalines = (line.rstrip() for line in data) 
        for line in datalines:
            # print(line)
            temp[index].append(math.log(float(line)))
            index = (index+1) %7
    return temp

def draw(plt, listt, c, my_label, add):

    avg = []
    mux = []
    for arr in listt:
        arr = (np.asarray(arr).astype(np.float))
        mux.append(arr)
        avg.append(stats.median(arr))

    plt.plot((np.array(range(len(mux)))*5 + add), np.asarray(
    avg).astype(np.float),color=c, label=my_label)

    bp = plt.boxplot(mux, positions=(np.array(
        range(len(mux)))*5 + add), sym='', widths=0.6, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor(c)

    

for ite in [16,36,49,64]:
    m_one = []
    m_two = []
    m_three = []
    m_one = makeArray('M1_data_' + str(ite) + '.txt')
    m_two = makeArray('M2_data_'+ str(ite)+'.txt')
    m_three = makeArray('M3_data_'+ str(ite)+'.txt')
    
    plt.figure()
    draw(plt, m_one, 'r', 'Multiple MPI_Sends(1 Element transmit)', -0.7)
    draw(plt, m_two, 'b', 'MPI_Pack/Unpack (Multiple element)', 0)
    draw(plt, m_three, 'g', 'MPI_Send using Derived datatypes(Vector)', 0.7)

    plt.title('Total Process(P) :: ' + str(ite))
    plt.xlabel('[ BoxPlot ] Total cell in matrix (data point) (N*N)')
    plt.ylabel('Time Scale [log(seconds)]')
    plt.xticks(np.arange(0, 7*5, 5), myLables)
    plt.legend()

    plt.savefig('plot'+str(ite)+'.jpg')

