CE4011 STRUCTURAL ANALYSIS ENGINE AND TESTS

*** The program is working inputting the text file name in the code as an example str_input.txt.
*** The inputs should be in the units: m, m^2, kN, kPa etc.
*** The test codes check directly the structural engine (str_engine_hw3.py). And the text files are read in the engine.

*** Object Oriented Structural Engine Txt File Format:

nNode=4
# trX/trY/rotZ: 1=Fixed, 0=Free
n=1	 x=0.0	y=0.0	trX=1	trY=1	rotZ=0
n=2	 x=0.0	y=3.0	trX=0	trY=0	rotZ=0
n=3	 x=4.0	y=3.0	trX=0	trY=0	rotZ=0
n=4  	 x=4.0  y=0.0	trX=0	trY=1	rotZ=0

nMaterial=2

n=1	area=0.02	inertia=0.08	elasticmod=200000000	
n=2	area=0.01	inertia=0.01	elasticmod=200000000

nMember=4
# type: frame/truss | rs/re: moment releases (1=pinned)
# udl: force per meter | pLoad: point load | pDist: distance from startnode
n=1  startnode=1  endnode=2  matProp=1  type=frame  rs=0  re=0  udl=0.0  pLoad=0.0  pDist=0.0
n=2  startnode=2  endnode=3  matProp=1  type=frame  rs=0  re=0  udl=0.0  pLoad=0.0  pDist=0.0
n=3  startnode=4  endnode=3  matProp=1  type=frame  rs=0  re=0  udl=0.0  pLoad=0.0  pDist=0.0
n=4  startnode=1  endnode=3  matProp=2  type=frame  rs=0  re=0  udl=0.0  pLoad=0.0  pDist=0.0

nLoads=2
n=1	nodeId=2	Fx=10	Fy=-10	Mz=0
n=2	nodeId=3	Fx=10	Fy=-10	Mz=0

***Example File Used for Case a:
# Structure (a): 
nNode=4
n=1 x=0.0 y=0.0 trX=0 trY=1 rotZ=0
n=2 x=0.0 y=5.0 trX=0 trY=0 rotZ=0
n=3 x=6.0 y=5.0 trX=0 trY=0 rotZ=0
n=4 x=6.0 y=0.0 trX=0 trY=1 rotZ=0

nMaterial=1
n=1 area=0.02 inertia=0.0005 elasticmod=200000000

nMember=3
n=1 startnode=1 endnode=2 matProp=1 type=frame rs=0 re=0
n=2 startnode=2 endnode=3 matProp=1 type=frame rs=0 re=0 udl=-10.0
n=3 startnode=3 endnode=4 matProp=1 type=frame rs=0 re=0

nLoads=1
nodeId=2 Fx=15.0 Fy=0.0 Mz=0.0


***Example File Used for Case c:

# Structure (c):
nNode=5
n=1 x=0.0 y=0.0 trX=1 trY=1 rotZ=1
n=2 x=0.0 y=5.0 trX=0 trY=0 rotZ=0
n=3 x=2.0 y=5.0 trX=0 trY=1 rotZ=0
n=4 x=4.0 y=0.0 trX=1 trY=1 rotZ=1
n=5 x=4.0 y=5.0 trX=0 trY=0 rotZ=0

nMaterial=1
n=1 area=0.02 inertia=0.0005 elasticmod=200000000

nMember=3
n=1 startnode=1 endnode=2 matProp=1 type=frame rs=0 re=0 udl=0
n=2 startnode=2 endnode=3 matProp=1 type=frame rs=0 re=0 udl=0
n=3 startnode=4 endnode=5 matProp=1 type=frame rs=0 re=0 udl=10.0

nLoads=2
nodeId=2 Fx=15.0 Fy=0.0 Mz=0.0
nodeId=5 Fx=0.0 Fy=0.0 Mz=0.0
