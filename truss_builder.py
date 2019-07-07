import warnings
import numpy as np

class Node:
    def __init__(self, name, position, slider=[True,True,True], spinner=[True,True,True], index=0):
        self.name = name
        self.index = index
        self.position = np.array(position)
        self.slider = slider
        self.spinner = spinner
        self.memberList = []
        self.loadSum = np.zeros((3))
        self.momentSum = np.zeros((3))

    def release(self):
        self.memberList = []
        
    def AddLoad(self, load):
        self.loadSum += load
    
    def AddMoment(self, moment):
        self.momentSum += moment

    def InitIndices(self):
        self.mindex = np.full((3),-1) #Row indices if the node transmits moments
        self.loadReactionIndex = np.full((3),-1)
        self.momentReactionIndex = np.full((3),-1)
        

class MemberNodeEntry:
    def __init__(self, node, hinged=[True,True,True]):
        self.node = node
        self.hinged = hinged
        self.loadIndex = -1

    def InitIndices(self):
        self.momentIndex = np.full((3),-1) #Column index for node moments
        

class Member:

    def __init__(self, name, refnode, hinged=[True,True,True]):
        self.name = name
        ne = MemberNodeEntry(refnode, hinged)
        self.nodeList = [ne]
        refnode.memberList.append(self)
        self.loadSum = np.zeros((3))
        self.momentSum = np.zeros((3))

    def release(self):
        for n in self.nodeList: n.node=None
        self.nodeList = []
        
    def AddNode(self, node, hinged=[True,True,True]):
        ne = MemberNodeEntry(node, hinged)
        self.nodeList.append(ne)
        node.memberList.append(self)
        
    def AddPointLoad(self, load, position):
        self.loadSum += load
        self.momentSum += np.cross(position - self.nodeList[0].node.position, load)
        
    def AddDistributedLineLoad(self, totalLoad, p1, p2):
        L1norm = np.linalg.norm(totalLoad)
        
        if(L1norm == 0): return
        self.loadSum += totalLoad
        p1 = p1 - self.nodeList[0].node.position
        p2 = p2 - self.nodeList[0].node.position
        pdiff = p2 - p1
        #print("pdiff ",pdiff)
        u_u = pdiff / np.linalg.norm(pdiff)
        px = p1 + (np.dot(p1,p1) - np.dot(p1,p2)) / np.dot(pdiff,pdiff) * pdiff
        #print("px ",px)
        rv = np.linalg.norm(px)
        
        if((px == rv*u_u).all()):
            #print("Identical")
            px = np.array([0,0,0])
            
            if(u_u[0] == 0):   v_u = np.array([1,0,0])
                
            elif(u_u[1] == 0): v_u = np.array([0,1,0])
            
            elif(u_u[2] == 0): v_u = np.array([0,0,1])
                
            else:    
                v_u = np.array([1,1, -(u_u[0] + u_u[1]) / u_u[2]])
                v_u /= np.linalg.norm(v_u)
        else:
            v_u = px / rv
            
        w_u = np.cross(u_u,v_u)
        
        #print("u_u ",u_u,", v_u ",v_u,", w_u ",w_u)
        
        u1 = np.linalg.norm(p1 - px)
        u2 = np.linalg.norm(p2 - px)
        L1_u = totalLoad / L1norm
        L1_u_w = np.dot(L1_u,w_u)
        
        self.momentSum += L1norm / np.linalg.norm(pdiff) * ((L1_u_w * rv * u_u - np.dot(L1_u,u_u) * rv * w_u) * (u2 - u1) + (-L1_u_w * v_u + np.dot(L1_u,v_u) * w_u) * (u2 * u2 - u1 * u1) / 2)
        

class Truss:

   def __init__(self, name):
       self.name = name
       self.allNodes = []
       self.allMembers = []

   def __del__(self):
       for m in self.allMembers: m.release()

       for n in self.allNodes: n.release()
       del self.allMembers
       del self.allNodes

   def CreateMember(self, name, refnode, hinged=[True,True,True]):
       member = Member(name, refnode, hinged)
       self.allMembers.append(member)
       return member

   def CreateNode(self, name, position, slider=[True,True,True], spinner=[True,True,True]):
       node = Node(name, position, slider, spinner, index=len(self.allNodes))
       self.allNodes.append(node)
       return node
      
   def GenTable(self):
       index = 0
       mindex = 0 #Row index for nodes that transmit moments, relative to nms*6+nns*3.
       nms = len(self.allMembers)
       nns = len(self.allNodes)
       a = np.zeros((nms * 6 + nns * 3,0))
       b = np.zeros(nms * 6 + nns * 3)

       for n in self.allNodes: n.InitIndices()
               
       for i,m in enumerate(self.allMembers):
           print('Member ',m.name)
           a = np.hstack((a,np.zeros((a.shape[0],len(m.nodeList) * 3))))

           for n in m.nodeList: n.InitIndices()
               
           refnode = n = m.nodeList[0]
           print('\tNode ',n.node.name,', slider: ',n.node.slider,', spinner: ',n.node.spinner,', hinged: ',n.hinged)
           n.loadIndex = index
           a[6*i][n.loadIndex] = 1
           a[6*i+1][n.loadIndex+1] = 1
           a[6*i+2][n.loadIndex+2] = 1
           index += 3   
           
           for j,d in enumerate(n.hinged):
               
               if d==False:
                   a = np.hstack((a,np.zeros((a.shape[0],1))))
                   n.momentIndex[j] = index
                   a[6*i+3+j][index] = 1
                   
                   if(n.node.mindex[j] == -1):
                       n.node.mindex[j] = mindex
                       mindex += 1  
                       a = np.vstack((a,np.zeros((1,a.shape[1]))))
                       b = np.hstack((b,-n.node.momentSum[j]))
                   a[6*nms+3*nns+n.node.mindex[j]][index] = -1
                   index += 1
                              
           a[6*nms+3*n.node.index][n.loadIndex] = -1
           a[6*nms+3*n.node.index+1][n.loadIndex+1] = -1
           a[6*nms+3*n.node.index+2][n.loadIndex+2] = -1
           
           for n in m.nodeList[1:]:
               print('\tNode ',n.node.name,', slider: ',n.node.slider,', spinner: ',n.node.spinner,', hinged: ',n.hinged)
               n.loadIndex = index
               #Member load equations (3 per member)
               a[6*i][n.loadIndex] = 1
               a[6*i+1][n.loadIndex+1] = 1
               a[6*i+2][n.loadIndex+2] = 1
               index += 3
               
               relpos = n.node.position - refnode.node.position
               #print('relpos ',relpos)
               #Member moment equations (3 per member)
               a[6*i+3][n.loadIndex+1] = -relpos[2]
               a[6*i+3][n.loadIndex+2] = relpos[1]
               a[6*i+4][n.loadIndex] = relpos[2]
               a[6*i+4][n.loadIndex+2] = -relpos[0]
               a[6*i+5][n.loadIndex] = -relpos[1]
               a[6*i+5][n.loadIndex+1] = relpos[0]
                       
               for j,d in enumerate(n.hinged):
               
                   if d==False:
                       a = np.hstack((a,np.zeros((a.shape[0],1))))
                       n.momentIndex[j] = index
                       a[6*i+3+j][index] = 1
                       
                       if(n.node.mindex[j] == -1):
                           n.node.mindex[j] = mindex
                           mindex += 1
                           a = np.vstack((a,np.zeros((1,a.shape[1]))))
                           b = np.hstack((b,-n.node.momentSum[j]))
                       #Create equations for moment carrying nodes (between 0 and 3 equations per node)
                       a[6*nms+3*nns+n.node.mindex[j]][index] = -1
                       index += 1
               #Node load equations (3 per node)       
               a[6*nms+3*n.node.index][n.loadIndex] = -1
               a[6*nms+3*n.node.index+1][n.loadIndex+1] = -1
               a[6*nms+3*n.node.index+2][n.loadIndex+2] = -1
           
           b[6*i:6*i+3] = -m.loadSum
           b[6*i+3:6*i+6] = -m.momentSum
           
       for n in self.allNodes:
           #print('Node ',n.name)
           b[6*nms+3*n.index:6*nms+3*n.index+3] = -n.loadSum
           
           for j,s in enumerate(n.slider):
               
               if s==False:
                   #print(j,' not slider')
                   n.loadReactionIndex[j] = a.shape[1]
                   a = np.hstack((a,np.zeros((a.shape[0],1))))
                   #Add to node load equations for sliders
                   a[6*nms+3*n.index+j][n.loadReactionIndex[j]] = 1
           
           for j,s in enumerate(n.spinner):
               
               if n.mindex[j]>=0 and s==False:
                   #print(j,' not spinner')
                   n.momentReactionIndex[j] = a.shape[1]
                   a = np.hstack((a,np.zeros((a.shape[0],1))))
                   #Add to node moment equations for spinners
                   a[6*nms+3*nns+n.mindex[j]][n.momentReactionIndex[j]] = 1
                   
       r=np.linalg.matrix_rank(a)
       n=a.shape[1]
       
       if(r < n):
           raise ValueError("The matrix of equations is underdetermined (the rank is smaller than the number of equations). There are thus an infinite number of solutions. This occurs when the user defines a structure whose loads and moments do not have a unique solution under the approximation that all members are rigid and do not stretch. Members attached to more than two nodes, as well as moment-carrying and fixed nodes should be re-evaluated by the user to resolve the issue.",UserWarning)
       
       if(r > n):
           warnings.warn("The matrix of equations is overdetermined. It is an abnormal situation that is likely to lead to a solution that does not satisfy all conditions. The user should re-evaluate the defined structure",UserWarning)
       #print(a)
       #print(b)

       i = 0

       while(i < a.shape[0]):

           if(all(a[i] == 0)):

               if(b[i] != 0):
                   raise ValueError("Non-physical equation. This typically occurs if a load or a moment is applied, but that no member or node allows to provide a reaction to this load or moment.")

               a = np.delete(a, i, axis=0)
               b = np.delete(b, i)

           else: i = i + 1

       return a, b
    
   def PrintResults(self,s):
       print('\nResults:\n================================')
       for i,m in enumerate(self.allMembers):
           print(i,': Member ',m.name)
           if(any(m.loadSum)): print('\tInput Applied load: ',m.loadSum)
           if(any(m.momentSum)): print('\tInput Applied moment: ',m.momentSum)
       
           for n in m.nodeList:
               print('\t',n.node.index,': Node ',n.node.name)
               print('\t\tApplied load:',s[n.loadIndex:n.loadIndex+3])
               print('\t\tApplied moment: [',s[n.momentIndex[0]] if n.momentIndex[0]>-1 else 0,', ',
                     s[n.momentIndex[1]] if n.momentIndex[1]>-1 else 0,', ',
                     s[n.momentIndex[2]] if n.momentIndex[2]>-1 else 0,']')        
           
       for n in self.allNodes:
           print(n.index,': Node ',n.name)
           if(any(n.loadSum)): print('\tInput Applied load: ',n.loadSum)
           if(any(n.momentSum)): print('\tInput Applied moment: ',n.momentSum)
           print('\tApplied load reaction: [',s[n.loadReactionIndex[0]] if n.loadReactionIndex[0]>-1 else 0,', ',
                     s[n.loadReactionIndex[1]] if n.loadReactionIndex[1]>-1 else 0,', ',
                     s[n.loadReactionIndex[2]] if n.loadReactionIndex[2]>-1 else 0,']')
           print('\tApplied moment reaction: [',s[n.momentReactionIndex[0]] if n.momentReactionIndex[0]>-1 else 0,', ',
                     s[n.momentReactionIndex[1]] if n.momentReactionIndex[1]>-1 else 0,', ',
                     s[n.momentReactionIndex[2]] if n.momentReactionIndex[2]>-1 else 0,']')    

