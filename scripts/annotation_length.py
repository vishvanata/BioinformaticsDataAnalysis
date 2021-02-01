

filename = "textfiles\entpoutu.txt"
fin = open(filename)
fin.readline()
annotation_readlines = fin.readlines()
count = 0
for line in annotation_readlines:
    count += 1
print(count)
