#/usr/bin/python2

import sys

lines = open(sys.argv[1])

avg_time = None 
count = 0

for line in lines:
  if line.startswith("SCALE"):
    if avg_time is not None:
      print(avg_time / count)
    avg_time = 0.0
    count = 0
    print(line.split())
  else:
    avg_time += float(line)
    count += 1
print(avg_time / count)
