##Implementation	of	Arrowhead	Algorithm	for	Domain Annotation

###Input	Data	Used
Input data can be found at http://chromosome.sdsc.edu/mouse/hi-c/download.html

###How to run
```
python arrowhead.py --input_data <inputfile> --is_normal <y/n> --apply_threshold1 <y/n> --t1 <threshold1_val1> --t2 <threshold1_val2> --t3 <threshold1_val3> --apply_threshold2 <y/n> --t4 <threshold2_val1> --t5 <threshold2_val2>
```
--is_normal: Is the input matrix normalized

threshold1_val1:	Value	for	variance	when	applying	threshold1

threshold1_val2:	Value	for	Mean(sgn(Ua,b))	when	applying	threshold1

threshold1_val3:	Value	for	Mean(sgn(La,b))	when	applying	threshold1

threshold2_val1:	Value	for	Mean(sgn(Ua,b))	when	applying	threshold2

threshold2_val1:	Value	for	Mean(sgn(La,b))	when	applying	threshold2

###Sample run
```
python arrowHead.py --input_data=uij.chr1 --is_normal n --t1 0.2 --t2 0.5 --t3 0.5 --apply_threshold1 y
```

For	performance	considerations	the	implementation	calculates	only	upper	triangle	values	of	all	the	matrices.	

On	running	the	script	it	also	produces	the	heat	map	for	Arrow	Head	matrix	as	A.jpg	and	for	the S	corner	matrix	as	Scorner.jpg

After	the	arrowhead.py	script	is	run,	we	get	a	file	called	ScornerData	which	is	the	filtered version	of	Scorner	matrix	and	we	need	to	find	the	connected	components	in	this	matrix.

To	find	the	connected	components,	run	the	script	connected_components.py
```
python	connected_components.py
```
