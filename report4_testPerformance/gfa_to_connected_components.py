import subprocess
from tqdm import tqdm
from collections import defaultdict
import sys
import pickle
import os

if len(sys.argv) < 1:
    exit("ex: python gfa_to_connected_components.py <gfa_path>")

gfa_path = sys.argv[1]

lines_number = int(subprocess.getoutput('wc -l ' + gfa_path).split()[0])

def find(x):
    l = leaders[x]
    if l is not None:
        l = find(l)
        leaders[x] = l
        return l
    return x

def union(x, y):
    lx, ly = find(x), find(y)
    if lx != ly:
        leaders[lx] = ly

print ("Parsing GFA file ...")

source = []
target = []
with open(gfa_path, 'r') as file:
    next(file)  # skip header
    for line in tqdm(file, total=lines_number):
        if line[0] == 'L':
            line = line.split()
            node_1 = int(line[1])
            node_2 = int(line[3])
            
            source.append(node_1)
            target.append(node_2)

        elif line[0] == "S":
            node = int(line.split()[1])
            source.append(node)
            target.append(node)

print ("Generating connected components...")
leaders = defaultdict(lambda: None)

for i in range(len(source)):
    union(source[i], target[i])


groups = defaultdict(set)

for x in leaders:
    groups[find(x)].add(x)

# pkl_file_name = os.path.basename(gfa_path).split(".")[0] + "_connected_components.pkl"
# print ("Writing " + pkl_file_name + " ...")

with open(gfa_path + ".components.tsv" ,'w') as comp:
    for k, v in groups.items():
        comp.write(",".join(list(map(str, v))) + '\n')

print(f"number of components: {len(groups)}")

#print(groups)

print("Done" , file = sys.stderr)

# with open(pkl_file_name, 'wb') as pickle_file:
#     pickle.dump(groups, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)
