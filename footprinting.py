import itertools
import xml.etree.ElementTree as ET
import timeit

#Generates a dictionary that contains every possible motif of a given lengt, each motif is assigned a None type value, to track whether that motif's parsimony value has been calculated yet
def generateTable(motifLength):
    motifs = list(itertools.product('acgt',repeat=motifLength))

    motifTable = {}
    for val in motifs:
        motifTable[''.join(val)] = None

    return motifTable

#Class used to store information of a given node in the phylogenetic tree such as it's nucelotide sequence (if any) and the name of the organism at that leaf. The table is the table of motifs that will store all this node's parsimony values. The distTable is the reverse of the table in that it uses the parsimony value as a key, with each value containing a list of motifs that have that parsimony value. This is important for the generation of new phases in the edge table algorithm. Edges contains a list of all EdgeTables for this node. It was assumed that this would always be 2, since almost all phylogenetic trees branch off in pairs.
class Node:
    def __init__(self, sequence, organism, motifTable):
        self.table = motifTable
        self.distTable = {}
        self.sequence = sequence
        self.organism = organism
        self.edges = []
        self.parent = -1
        self.children = []

#Calculate's all the parsimony values for a given node. If the node is a leaf (as determined by  it not having an empty string as its sequence), a motif's parsimony value is 0 if the motif is present in the sequence, and infinity if it is not. Non-leaf nodes have a more complicated algorithm explained in the EdgeTable class.
    def calculateParsimonyValue(self, d):
        if self.sequence != '':
            for motif in self.table:
                parsimonyVal = float('inf')
                if motif in self.sequence:
                    parsimonyVal = 0
                
                self.table[motif] = parsimonyVal
                if parsimonyVal in self.distTable:
                    self.distTable[parsimonyVal].append(motif)
                else:
                    self.distTable[parsimonyVal] = [motif]
        else:
            phase = 0

#Load motifs with a parsimony value of 0 into the queue's of both edges so that the first phase is loaded up (from there each phase generated the next phases values)
            if 0 in self.edges[0].childDistTable:
                for newMotif in self.edges[0].childDistTable[0]:
                    self.edges[0].queue.append(newMotif)
                    self.edges[0].queueDict[newMotif] = 0
                self.edges[0].childDistTable[0] = []
            
            if 0 in self.edges[1].childDistTable:
                for newMotif in self.edges[1].childDistTable[0]:
                    self.edges[1].queue.append(newMotif)
                    self.edges[1].queueDict[newMotif] = 0
                self.edges[1].childDistTable[0] = []

            while len(self.edges[0].queue) > 0 or len(self.edges[1].queue) > 0 or phase <= d:
                self.edges[0].newPhase(phase, d)
                self.edges[1].newPhase(phase, d)
                
                phase += 1

#Adds a child to the current node, also generating an EdgeTable for it
    def addChild(self, key, childDistTable, childTable, motifLength):
        self.children.append(key)
        self.edges.append(EdgeTable(generateTable(motifLength) ,childDistTable, self.distTable, self.table, childTable))

#This class is where the determination of parsimony values in non-leaf nodes occurs. One of the key optimizations proposed in the paper is to generate a table for each edge and search all the motifs of those tables in phases that significantly reduces the amount of computations. Each phase (p) is comprised of all the motifs of the childDistTable that have a parsimony value of p and that are adjacent to the motifs of the previous phase. The phase number is the parsimony value of every motif in the phase. However, before considering a value to be part of a phase (and using it to expand to the queue of the next phase), a check is done to make sure that it is in bound (through the isInBounds method). The queue is what stores all the motifs of a given phase. queueDict is a dictionary version of the queue, while it might seem more inefficient to have 2 data structures to represent the queue, they both serve different important purposes for optimizing the runtime. A queue is typically best suited for a list because of their ordered nature, but for every new motif to be added to the queue, a check has to be done to make sure that said item isn't already in this queue. This lookup operation is much faster using a dictionary, and so using a dictionary pays off in dividends (the runtime from one of the tests went from over 10 minutes to 50 seconds). 
class EdgeTable:
    def __init__(self, table, childDistTable, parentDistTable, parentTable, childTable):
        self.table = table
        self.visited = {}
        self.childDistTable = childDistTable
        self.parentDistTable = parentDistTable
        self.childTable = childTable
        self.parentTable = parentTable
        self.neighborTable = -1
        self.neighborQueue = -1
        self.queue = []
        self.queueDict = {}

    def newPhase(self, phase, d):
        for i in range(len(self.queue)):
            motif = self.queue.pop(0)
            self.queueDict.pop(motif)
            self.visited[motif] = 0
            if self.isInBounds(phase, motif, d):
                self.table[motif] = phase

#isDone is used to track whether the other edge has calculated a value for this motif, if so, the two values can be added together and set as the final parsimony value for the motif in the parentTable
                isDone = False
                if not self.neighborTable[motif] is None:
                    self.parentTable[motif] = phase + self.neighborTable[motif]
                    isDone = True

#if the motif isDone for the parent node, then we also need to add it to the parent's distTable
                if isDone:
                    if self.parentTable[motif] in self.parentDistTable:
                        self.parentDistTable[self.parentTable[motif]].append(motif)
                    else:
                        self.parentDistTable[self.parentTable[motif]] = [motif]

                bases = ['a','c','t','g']

#Generate every motif that is a Hamming distance of 1 away from the current motif, check whether it has been visited, is in the queue, or has a value of None in the childTable (indicating that it has been pruned from the tree). If so, add it to the queue for the next phase. 
                for i in range(len(motif)):
                    baseIndex = bases.index(motif[i])
                    
                    for j in range(1,4):
                        newIndex = (baseIndex + j)%4
                        newChar = bases[newIndex]
                        newMotif = motif[:i] + newChar + motif[i+1:]
                        if not newMotif in self.visited and not newMotif in self.queueDict and not self.childTable[newMotif] is None:
                            self.queue.append(newMotif)
                            self.queueDict[newMotif] = 0

#If there are motifs in the childDistTable that have a parsimony value equal to the value of the next phase, it needs to be added to the phase, assuming it meets the above described conditions                
        if phase+1 in self.childDistTable:
            for newMotif in self.childDistTable[phase+1]:
                if not newMotif in self.visited and not newMotif in self.queueDict and not self.childTable[newMotif] is None:
                    self.queue.append(newMotif)
                    self.queueDict[newMotif] = 0
            self.childDistTable[phase+1] = []

#Checks whether a given motif is in bounds using the sibling bounding method proposed in the paper. By processing the edge tables of a node in parallel, doing each edge's phases at the same time, it allows for us to determined whether a given motif's parsimony value will surpass d, the cap for the value that is set by the user. In the end, for each motif in a node, it's parsimony value is calculated as the sum of its values in the edge tables (essentially this aggregates the value by taking into account both child nodes). The parsimony value of the motif might not yet have been calculated in the neighbor edges table, if that is the case, we know that the value will be at least phase + 1 (because we check to make sure it is not in the current phase). Then we check to make sure that d >= phase + neighborVal, if so then we know there combined parsimony value is in the bounds of the d value
    def isInBounds(self, phase, motif, d):
        neighborVal = self.neighborTable[motif] 

        if motif in self.neighborQueue:
            neighborVal = phase

        if neighborVal is None:
            neighborVal = phase + 1
        
        if d >= phase + neighborVal:
            return True
        else:
            return False

def main(xmlFile, motifLength, d):
#Open and parse xml file that contains phylogenetic tree
    xmlTree = ET.parse(xmlFile)
    root = xmlTree.getroot()

#Initialize the tree structure, each key in the tree contains the Node class
    tree = {0:Node('','',generateTable(motifLength))}

#treeQueue is used in the construction of the tree from the xml file, it contains all the children of each node that is detected by the xml parser. It is initialized with the children that are in the root node. Only 'clade' nodes are added to the tree. Queue is structured as a nested list, containing the child node and the key that is used to store it in the tree structure. The key is used to allow it's children to add it to their parent list and to let them add themselves to this node's child list. 
    treeQueue = []
    for child in root:
        if child.tag == 'clade':
            treeQueue.append([child, 0])

    calcQueue = []
    i = 1
#Build the tree, adding each node along with its parents and children, the value i is used to generate a different key for each node in the tree. 
    while len(treeQueue) > 0:
        node = treeQueue.pop(0)
        tree[i] = Node('','',generateTable(motifLength))
        tree[i].parent = node[1]
        tree[node[1]].addChild(i, tree[i].distTable, tree[i].table, motifLength)

#Setup the edges so that they have access to their neighbors queueDict and parsimony table 
        if tree[i].parent != -1:
            parentIndex = tree[i].parent
            if len(tree[parentIndex].edges) > 1:
                tree[parentIndex].edges[0].neighborTable = tree[parentIndex].edges[1].table
                tree[parentIndex].edges[1].neighborTable = tree[parentIndex].edges[0].table

                tree[parentIndex].edges[0].neighborQueue = tree[parentIndex].edges[1].queueDict
                tree[parentIndex].edges[1].neighborQueue = tree[parentIndex].edges[0].queueDict

#check if the node is a parent based on whether it contains 'clade' children nodes, if not, then it must be a leaf, which is where the algorithm starts, so it is added to the calcQueue
        isParent = False
        for child in node[0]:
            if child.tag == 'clade':
                treeQueue.append([child,i])
                isParent = True
            if child.tag == 'sequence':
                tree[i].sequence = child.text
            if child.tag == 'taxonomy':
                tree[i].organism = child.text
        
        if not isParent:
            calcQueue.append(i)

        i += 1

    while len(calcQueue) > 0:
        node = calcQueue.pop(0)
        tree[node].calculateParsimonyValue(d)
        if tree[node].parent != -1 and not tree[node].parent in calcQueue:
            calcQueue.append(tree[node].parent)
    
        calcQueue.sort(reverse=True)

    print(tree[0].distTable)

#main('test.xml',4, 2)

print(timeit.timeit('main(\'fullTree.xml\', 8, 9)', 'from __main__ import main',number=1))