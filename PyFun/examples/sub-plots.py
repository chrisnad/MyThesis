import matplotlib.pyplot as plt
f1 = plt.figure(1)
ax1 = f1.add_subplot(111)
ax1.plot(range(0,20), label="Bar 1")
plt.xlabel("ooo")
plt.grid(True)
plt.legend(fancybox=True).get_frame().set_alpha(0.5)
#f1.savefig('/home/cnader/Desktop/point-cet-aprem/oulala.png')

f2 = plt.figure(2)
ax2 = f2.add_subplot(111)
ax2.plot(range(10,20), label="asdasda")
plt.xlabel("olaaoo")
plt.legend(fancybox=True).get_frame().set_alpha(0.5)



plt.show()
