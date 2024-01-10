import multiprocessing as mp
import os

def make_file (fname, s):
    f = open(fname, 'w')
    print >> f, s
    f.close()

pool = mp.Pool(processes=4)
for i in range(4):
    pool.apply(make_file, args=(str(i)+"_temp.txt", str(i)))

for i in range(4):
    os.system("cat " + str(i) + "_temp.txt" + " >> final.txt")
    os.system("rm " + str(i) + "_temp.txt")



