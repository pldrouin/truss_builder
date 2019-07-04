import numpy as np
import truss_builder as tb

#Simple truss with three members, down load applied on top node

n1 = tb.Node('n1',[0,1,0])
n2 = tb.Node('n2',[-1,0,0],slider=[True,False,True],spinner=[True,True,True])
n3 = tb.Node('n3',[1,0,0],slider=[True,False,True],spinner=[True,True,True])
m1 = tb.Member('m1',n1)
m1.AddNode(n2,hinged=[True,True,True])
m2 = tb.Member('m2',n1)
m2.AddNode(n3,hinged=[True,True,True])
m3 = tb.Member('m3',n2,hinged=[True,True,True])
m3.AddNode(n3,hinged=[True,True,True])
n1.AddLoad([0,-10,0])
#m1.AddDistributedLineLoad([1,2,2],[0.99,1,1],[1.01,1,1])
#m1.AddPointLoad([0,2,0], [0,0,0])
#m1.AddDistributedLineLoad([0,2,0],[-0.0001,0,0],[0.0001,0,0])
#m1.momentSum
a,b = tb.GenTable()
#print(a)
#print(b)
#print("Matrix A shape is ",a.shape)
#print("Number of rows with zeros: ",np.count_nonzero(a,axis=0))
#print("Number of columns with zeros:\n",np.count_nonzero(a,axis=1))
#print("Number of columns with zeros for b:\n",(b!=0)*1)
#print("Matrix rank is ",np.linalg.matrix_rank(a))
s = np.linalg.lstsq(a,b,rcond=None)[0]

tb.PrintResults(s)
np.matmul(a,s)-b
