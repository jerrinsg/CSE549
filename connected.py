import numpy as np
from scipy.ndimage.measurements import label

a = np.loadtxt('ScornerMatrix')

labeled_array, num_features = label(a)
print("number :",num_features)
print(labeled_array)
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
    #print len(c[0])
    count+=1
    #print c

print "Average length ",noofcells/count
print "min ",min
print "max ",max
