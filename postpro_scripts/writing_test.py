f = open("t.txt", "w")

print >> f, "hello"

for i in range(10):
    f.write(str(i+1))

print >> f
    
f.close()
