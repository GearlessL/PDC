import os
import sys
import signal
import time
import subprocess
import threading

def shellGetOutput(str1) :
  lock = threading.Lock()
  with lock:
    try:
      result = subprocess.run(str1, stdout=subprocess.PIPE,shell=True, timeout=1e5)
      # Output results
      print(result.stdout.decode())
      appendToFile(result.stdout.decode(), out_filename)
    except subprocess.TimeoutExpired:
      print("The program has been running for longer than the specified time and has been terminated.")
  


def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

alg_names = ["", "Peeling", "AC", "SC", "ParPeel", "ParPeel+", "Shell-PDC"]


files = ["em"];
num_threads = [32]
algs = [1, 2, 3, 4, 5, 6]
read_dir = "./materials/"
for filename in files:
  for nw in num_threads:
    for al in algs:
      out_filename = "./result/" + filename + "/res.out"
      ss = ("./dcore" + " -t "+ str(nw) + " -f " + read_dir  + filename +".txt" + " -a " + str(al))
      
      info = "file: " + filename + "; in: " + read_dir
      appendToFile(info, out_filename)
      info = "algs: " + str(alg_names[al]) + "; threads: " + str(nw) + ";\n"
      appendToFile(info, out_filename)
      
      shellGetOutput(ss)
      
      appendToFile("\n", out_filename)
            
            
            
