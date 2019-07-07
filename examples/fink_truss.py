import numpy as np
import truss_builder as tb

#King post truss with five members, uniform load on rafters

d0 = 12
s = 2 * 96
H = 7 * 12
dK = s/2
rsep = 46.375
tsep = 61.833/2

TCL = 139.6/12
BCL = 40/12

LR1x = 0
LR1y = -TCL * (s/2 + d0)
LR2x = 0
LR2y = -TCL * (s/2 + d0)
WKy = -0
WTy = -BCL * s
WP1y = -0
WP2y = -0

tr = tb.Truss("Simple Truss")

nRK  = tr.CreateNode('RK',[0,dK,0])
nR1 = tr.CreateNode('R1',[-s/2-d0,-d0,0])
nR2 = tr.CreateNode('R2',[s/2+d0,-d0,0])
nRT1 = tr.CreateNode('RT1',[-s/2,0,0],slider=[True,False,True])
nR1div = tr.CreateNode('R1div',[-s/8,-s/8+s/2,0])
nRT2 = tr.CreateNode('RT2',[s/2,0,0],slider=[True,False,True])
nTK  = tr.CreateNode('TK',[0,0,0])
nR1sep = tr.CreateNode('R1sep',[-rsep,-rsep+s/2,0])
nR2sep = tr.CreateNode('R2sep',[rsep,-rsep+s/2,0])
nT1sep = tr.CreateNode('T1sep',[-tsep,0,0])
nT2sep = tr.CreateNode('T2sep',[tsep,0,0])
#nPB1 = tr.CreateNode('PB1',[-s/2,-H,0],slider=[True,False,True])
#nPB2 = tr.CreateNode('PB2',[s/2,-H,0],slider=[True,False,True])

#mR1a = tr.CreateMember('Rafter 1a (R1a)',nR1)
#mR1a.AddNode(nRT1)
#mR1a.AddNode(nR1sep)
#mR1a.AddDistributedLineLoad([0,LR1y*(d0+s/2-rsep)/(d0+s/2),0],nR1.position,nR1sep.position)

#mR1c = tr.CreateMember('Rafter 1c (R1c)',nR1sep)
#mR1c.AddNode(nRK)
#mR1c.AddDistributedLineLoad([0,LR1y*rsep/(d0+s/2),0],nR1sep.position,nRK.position)

mR1a = tr.CreateMember('Rafter 1a (R1a)',nR1)
mR1a.AddNode(nRT1)
mR1a.AddNode(nR1sep)
mR1a.AddDistributedLineLoad([0,LR1y*(d0+s/2-rsep)/(d0+s/2),0],nR1.position,nR1sep.position)

mR1b = tr.CreateMember('Rafter 1b (R1b)',nR1sep)
mR1b.AddNode(nR1div, hinged=[False,False,False])
mR1b.AddDistributedLineLoad([0,LR1y*(rsep-s/8)/(d0+s/2),0],nR1sep.position,nR1div.position)

mR1c = tr.CreateMember('Rafter 1c (R1c)',nR1div, hinged=[False,False,False])
mR1c.AddNode(nRK)
mR1c.AddDistributedLineLoad([0,LR1y*(s/8)/(d0+s/2),0],nR1div.position,nRK.position)

mR2a = tr.CreateMember('Rafter 2a (R2a)',nR2)
mR2a.AddNode(nRT2)
mR2a.AddNode(nR2sep)
mR2a.AddDistributedLineLoad([0,LR2y*(d0+s/2-rsep)/(d0+s/2),0],nR2.position,nR2sep.position)

mR2c = tr.CreateMember('Rafter 2c (R2c)',nR2sep)
mR2c.AddNode(nRK)
mR2c.AddDistributedLineLoad([0,LR2y*rsep/(d0+s/2),0],nR2sep.position,nRK.position)

mT1a = tr.CreateMember("Tie Beam 1a (T1a)",nRT1)
mT1a.AddNode(nT1sep)
mT1a.AddDistributedLineLoad([0,WTy/2*(s/2-tsep)/(s/2),0],nRT1.position,nT1sep.position)

mT = tr.CreateMember("Tie Beam centre (Tcentre)",nT1sep)
mT.AddNode(nT2sep)
mT.AddDistributedLineLoad([0,WTy*(2*tsep)/s,0],nT1sep.position,nT2sep.position)

mT2a = tr.CreateMember("Tie Beam 2a (T2a)",nRT2)
mT2a.AddNode(nT2sep)
mT2a.AddDistributedLineLoad([0,WTy/2*(s/2-tsep)/(s/2),0],nRT2.position,nT2sep.position)

mM1a = tr.CreateMember("Member 1a (1a)",nRK)
mM1a.AddNode(nT1sep)

mM1b = tr.CreateMember("Member 1b (1b)",nT1sep)
mM1b.AddNode(nR1sep)

mM2a = tr.CreateMember("Member 2a (2a)",nRK)
mM2a.AddNode(nT2sep)

mM2b = tr.CreateMember("Member 2b (2b)",nT2sep)
mM2b.AddNode(nR2sep)

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

