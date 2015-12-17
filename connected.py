import numpy as np
from scipy.ndimage.measurements import label

print "Loading ScornerData"
a = np.loadtxt('ScornerData')
print "ScornerData Complete"

print "Connected components labelling begin"
labeled_array, num_features = label(a)

b = np.array(labeled_array)
noofcells=0
count=0
max=0
min=100000000
for i in range(1,num_features+1) :
    c = np.where(b==i)
    noofcells += len(c[0])
    if len(c[0])>max:
        max=len(c[0])
    if len(c[0])<min:
        min=len(c[0])
    count+=1

print "Number of connected components :",count
print "Average size of a connected components :",noofcells/count
print "Minimum size of a connected component :",min
print "Maximum size of a connected component :",max

