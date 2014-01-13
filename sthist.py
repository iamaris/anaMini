import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


#f = open('st.txt','r')
#print f.readlines()
#f.close()
st = []
x1 = []
x2 = []
x3 = []
x4 = []
x5 = []
x6 = []
x7 = []
x8 = []
a = np.genfromtxt("st.txt")
for i in a:
  st.append(i[0])
  x1.append(i[1])
  x2.append(i[2])
  x3.append(i[3])
  x4.append(i[4])
  x5.append(i[5])
  x6.append(i[6])
  x7.append(i[7])
  x8.append(i[8])

s1 = np.array(x1)
s2 = np.array(x2)
s3 = np.array(x3)
s4 = np.array(x4)
s5 = np.array(x5)
s6 = np.array(x6)
s7 = np.array(x7)
s8 = np.array(x8)

s23 = s2+s3
y4 = np.divide(s4,s23)
y5 = np.divide(s5,s23)
y6 = np.divide(s6,s23)
y7 = np.divide(s7,s23)
y8 = np.divide(s8,s23)
e4 = np.sqrt(np.divide(np.sqrt(s4),s4)*np.divide(np.sqrt(s4),s4)+np.divide(np.sqrt(s23),s23)*np.divide(np.sqrt(s23),s23))
e5 = np.sqrt(np.divide(np.sqrt(s5),s5)*np.divide(np.sqrt(s5),s5)+np.divide(np.sqrt(s23),s23)*np.divide(np.sqrt(s23),s23))
e6 = np.sqrt(np.divide(np.sqrt(s6),s6)*np.divide(np.sqrt(s6),s6)+np.divide(np.sqrt(s23),s23)*np.divide(np.sqrt(s23),s23))
e7 = np.sqrt(np.divide(np.sqrt(s7),s7)*np.divide(np.sqrt(s7),s7)+np.divide(np.sqrt(s23),s23)*np.divide(np.sqrt(s23),s23))
e8 = np.sqrt(np.divide(np.sqrt(s8),s8)*np.divide(np.sqrt(s8),s8)+np.divide(np.sqrt(s23),s23)*np.divide(np.sqrt(s23),s23))

#z23 = plt.errorbar(st, s23, yerr=np.sqrt(s23),fmt='o')
#z4 = plt.errorbar(st, y4, yerr=e4,fmt='o')
#z5 = plt.errorbar(st, y5, yerr=e5,fmt='o')
#z6 = plt.errorbar(st, y6, yerr=e6,fmt='o')
#z7 = plt.errorbar(st, y7, yerr=e7,fmt='o')
##z8 = plt.errorbar(st, y8, yerr=e8,fmt='o')

plt.errorbar(st, y4, fmt='o')
plt.errorbar(st, y5, fmt='o')
plt.errorbar(st, y6, fmt='o')
plt.errorbar(st, y7, fmt='o')
plt.errorbar(st, y8, fmt='o')

plt.xlabel('St')
plt.ylabel('# Event')
plt.title('St')
# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()
