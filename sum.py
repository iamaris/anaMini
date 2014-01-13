f = open('lumi.txt','r')
total = []
for i in f.readlines():
  total.append(float(i))
f.close()

ans = 0.
for i in total:
  ans = ans + i

print ans

