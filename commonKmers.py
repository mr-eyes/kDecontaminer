import sys
import kProcessor as kp

kf_file = sys.argv[1].replace(".mqf","")

kf = kp.kDataFrame.load(kf_file)

it = kf.begin()

size = 0

while it != kf.end():
    size += 1
    it.next()

print(f"total size: {size}")