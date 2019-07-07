import numpy as np
import truss_builder as tb

#King post truss with five members, uniform load on rafters

theta = np.pi / 4
d0 = 12
s = 2 * 48
H = 7 * 12
dK = (s/2)/np.tan(theta)

TCL = 16.13333333333333
BCL = 4

LR1x = 0
LR1y = -TCL * (s/2 + d0)
LR2x = 0
LR2y = -TCL * (s/2 + d0)
WKy = -0
WTy = -BCL * s
WP1y = -0
WP2y = -0

dsq = 3.5/2
ddia = dsq*np.sqrt(2)

tr = tb.Truss("Simple Truss")

nRK  = tr.CreateNode('RK',[0,dsq+dK+ddia,0])
nR1 = tr.CreateNode('R1',[-s/2-d0,-d0+ddia,0])
nR2 = tr.CreateNode('R2',[s/2+d0,-d0+ddia,0])
nRT1 = tr.CreateNode('RT1',[-s/2-ddia,dsq,0],slider=[True,False,True])
nR1div = tr.CreateNode('R1div',[-s/4,s/4,0])
nRT2 = tr.CreateNode('RT2',[s/2+ddia,dsq,0],slider=[True,False,True])
nTK  = tr.CreateNode('TK',[0,dsq,0])
#nPB1 = tr.CreateNode('PB1',[-s/2,-H,0],slider=[True,False,True])
#nPB2 = tr.CreateNode('PB2',[s/2,-H,0],slider=[True,False,True])

#mR1 = tr.CreateMember('Rafter 1 (R1)',nR1)
#mR1.AddNode(nRT1)
#mR1.AddNode(nRK)
#mR1.AddDistributedLineLoad([0,LR1y,0],nR1.position,nRK.position)

mR1a = tr.CreateMember('Rafter 1a (R1a)',nR1)
mR1a.AddNode(nRT1)
mR1a.AddNode(nR1div, hinged=[False,False,False])
mR1a.AddDistributedLineLoad([0,LR1y*(d0+s/4)/(d0+s/2),0],nR1.position,nR1div.position)

mR1b = tr.CreateMember('Rafter 1b (R1b)',nR1div, hinged=[False,False,False])
mR1b.AddNode(nRK)
mR1b.AddDistributedLineLoad([0,LR1y*(s/4)/(d0+s/2),0],nR1div.position,nRK.position)

mR2 = tr.CreateMember('Rafter 2 (R2)',nR2)
mR2.AddNode(nRT2)
mR2.AddNode(nRK)
mR2.AddDistributedLineLoad([0,LR2y,0],nR2.position,nRK.position)

mT1 = tr.CreateMember("Tie Beam 1 (T1)",nRT1,hinged=[True,True,True])
mT1.AddNode(nTK)
mT1.AddDistributedLineLoad([0,WTy/2,0],nRT1.position,nTK.position)

mT2 = tr.CreateMember("Tie Beam 2 (T2)",nRT2,hinged=[True,True,True])
mT2.AddNode(nTK)
mT2.AddDistributedLineLoad([0,WTy/2,0],nRT2.position,nTK.position)

mK = tr.CreateMember("King Post (K)",nRK)
mK.AddNode(nTK)
mK.AddDistributedLineLoad([0,WKy,0],nRK.position,nTK.position)

#mP1 = tr.CreateMember("Post 1 (P1)",nRT1)
#mP1.AddNode(nPB1)
#mP1.AddDistributedLineLoad([0,WP1y,0],nRT1.position,nPB1.position)

#mP2 = tr.CreateMember("Post 2 (P2)",nRT2)
#mP2.AddNode(nPB2)
#mP2.AddDistributedLineLoad([0,WP2y,0],nRT2.position,nPB2.position)

a,b = tr.GenTable()
#print(a)
#print(b)
#print("Matrix A shape is ",a.shape)
#print("Number of rows with zeros: ",np.count_nonzero(a,axis=0))
#print("Number of columns with zeros:\n",np.count_nonzero(a,axis=1))
#print("Number of columns with zeros for b:\n",(b!=0)*1)
#print("Matrix rank is ",np.linalg.matrix_rank(a))
s = np.linalg.lstsq(a,b,rcond=None)[0]

tr.PrintResults(s)
print(np.matmul(a,s)-b)

