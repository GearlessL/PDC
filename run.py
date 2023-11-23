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
      # result = os.system(str1)
      # 输出结果
      print(result.stdout.decode())
      appendToFile(result.stdout.decode(), out_filename)
    except subprocess.TimeoutExpired:
      print("程序运行时间超过指定的时间，已被终止。")
  


def appendToFile(out, filename):
  with open(filename, "a+") as out_file:
    out_file.writelines(out)

alg_names = ["", "Peeling", "AC", "SC", "PDC", "PDC+", "Fast-PDC"]


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
            
            
            