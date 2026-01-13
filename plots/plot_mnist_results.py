import matplotlib.pyplot as plt

results = [
(250,1.2024157219131104),
(500,1.210241753228835),
(750,1.2351571916674182),
(1500,1.256847574740919),
(2000,1.251889585855637)]

x = [2*r[0] for r in results]
y = [r[1] for r in results]

plt.ylim(0,2)
plt.plot(x,y,marker='o')
plt.title('Approximation Ratio on Mnist Data')
plt.xlabel('Sample Size')
plt.ylabel('Approximation Ratio (Our cost / min cost)')
plt.show()
