import csv
import matplotlib.pyplot as plt

our_data = []
pouyan_data = []
with open('compare_cluster_emd_ratios_to_lahnetal') as file:
    reader = csv.reader(file, delimiter='\t',)
    for row in reader:
        our_data.append(row)

with open('pouyan_cluster_emd_ratio_data') as file:
    reader = csv.reader(file, delimiter='\t',)
    for row in reader:
        pouyan_data.append(row)

our_samples = [ 2*int(row[0]) for row in our_data[1:] ]
pouyan_samples = [ 2*int(row[0]) for row in pouyan_data[1:] ]

our_ratio = [ float(row[9]) for row in our_data[1:6] ]
pouyan_ratio = []
for k in range((len(pouyan_data)-1)//5):
    val = sum([float(pouyan_data[5*k + 1 + i][8]) for i in range(5)])/5
    pouyan_ratio.append(val)

#pouyan_ratio = [ sum([float(row[8]) for row in [ pouyan_data[5*k+1:5*(k+1)+1]])/5 for k in range((len(pouyan_data)-1)/5)]]
ratio = [ pouyan_ratio[k]/our_ratio[k] for k in range(len(our_ratio)) ]

#plt.scatter(our_samples,our_ratio, marker='x',label='this paper')
#plt.scatter(pouyan_samples,pouyan_ratio, marker='o',label='Lahn et al')
plt.plot(our_samples[0:5], ratio,marker='o')
plt.xlabel('Sample Size')
plt.ylim(0, 5.5)
plt.title('Ratio of Min Cost with HST vs. Cluster Distances')
plt.ylabel('Ratio')
plt.legend()
plt.savefig('ratio_of_min_cost.png')
plt.show()
