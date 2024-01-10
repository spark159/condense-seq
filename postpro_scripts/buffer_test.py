import time
import math

start = time.time()
f = open("test.txt", 'w')
s = ""
buffer_size = 100000000
count = 0
for i in xrange(500000000):
    #print >> f, 'This is a speed test'
    #f.write('This is a speed test\n')
    s += 'This is a speed test\n'
    count += 1
    if count > buffer_size:
        count = 0
        f.write(s)
        s = ""
f.write(s)
f.close()
end = time.time()
print(end - start)
