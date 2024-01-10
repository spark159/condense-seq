import os, sys, subprocess, re
 


samtools_cmd = ["module load samtools", "&&", "samtools", "view", "/home/spark159/scratch4-tha4/sangwoo/2022_07_04_GM_NCP/GM-NCP-PEG-1.bam"]
samtools_proc = subprocess.Popen(" ".join(samtools_cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

for line in samtools_proc.stdout:
    print line
    break
#stdout, error = samtools_proc.communicate()

#print stdout, error
#subprocess.call(" ".join(samtools_cmd), shell=True)

#os.system('module load samtools/1.15.1 ')
#subprocess.call("samtools view -h", shell=True)
#subprocess.check_output('module load samtools/1.15.1 && samtools  view -h', shell=True)
#subprocess.call('module load samtools/1.15.1 | samtools  view /home/spark159/scratch4-tha4/sangwoo/2022_07_04_GM_NCP/GM-NCP-PEG-1.bam | head', shell=True)


#out = subprocess.check_output(['samtools','view', '-h'], shell=True, universal_newlines=True)
#samtools_cmd = ["samtools",  "view", "-h"]
#subprocess.call("PATH=$PATH:.", shell=True)
#subprocess.call(". shell_scrpt.sh", shell=True)
#subprocess.call(samtools_cmd)
#subprocess.call("ml samtools | samtools view -h", shell=True)
#subprocess.call("ml samtools | samtools view /home/spark159/scratch4-tha4/sangwoo/2022_07_04_GM_NCP/GM-NCP-PEG-1.bam | head", shell=True)
#subprocess.call("ls -l", shell=True)

