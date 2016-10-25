def displayTree(tree):
	'''
	displays FOGSAA tree dictionary
	'''
	#recursive function for displaying branches
	def printBranches(node, level):

		if node[5]==[]:
			return ''
		else:
			nodelist = ''
			indent = '  ' * level
			for c in node[5]:
				nodelist += indent + '[' + str(c[1]) + ',' + str(c[2]) + ']\n' + printBranches(c,level+1)
			return nodelist
	#create base
	print '[0,0]'
	currentnode = (0,0)
	return printBranches(tree[currentnode],1)

def recreateBestBranch(tree,bestpathnode):
	'''
	recreates the best path chosen by the alogrithm as a list of tuples containing pointer positions P1,P2
	'''
	bestalignment = [(bestpathnode[1],bestpathnode[2])]
	currentnode = tree[(bestpathnode[1],bestpathnode[2])]
	while currentnode[1] != 0 or currentnode[2] != 0:
		if currentnode[3] == 1:
			#if child is of type 1
			bestalignment.append((currentnode[1]-1,currentnode[2]-1))
			currentnode = tree[(currentnode[1]-1,currentnode[2]-1)]

		elif currentnode[3] == 2:
			#type 2
			bestalignment.append((currentnode[1]-1,currentnode[2]))
			currentnode = tree[(currentnode[1]-1,currentnode[2])]
			
		elif currentnode[3] == 3:
			#type 3
			bestalignment.append((currentnode[1],currentnode[2]-1))
			currentnode = tree[(currentnode[1],currentnode[2]-1)]
	return bestalignment


def formatSeq(dnasequence):
	"""
	takes in a DNA sequence and converts it to all capital letters, removes spaces, 
	and throws an exception when characters other than ATCG are present.
	Input: dnasequence is a string containing the sequence of interest. 
	Output: copy of dnasequence with all caps and no spaces.
	"""
	
	#remove all spaces and convert to uppercase
	editedsequence = dnasequence.replace(' ','')
	editedsequence = editedsequence.upper()

	#check for unwanted characters
	for char in editedsequence:
		if char not in ['A','T','C','G']:
			raise ValueError('Sequence contains non-nucleotide characters. Please remove any characters other than A,T,C,G.')
	return editedsequence

def findF(x1,x2,M,Ms,G):
	#Finds the Fitness Score of a FOGSAA tree node; returns Fmin and Fmax values
	if x2 < x1:
		Fmin = x2*Ms+G*(x1-x2)
		Fmax = x2*M+G*(x1-x2)
	else:
		Fmin = x1*Ms+G*(x2-x1)
		Fmax = x1*M+G*(x2-x1)
	return Fmin,Fmax

def pluckNode(q):
	#get the highest priority node from the prioity queue, q
	#returns the new node
	newnode = q.get()
	return newnode[1]

def FOGSAA(S1,S2, M,Ms,G):
	#Implements FOGSAA (Chakraborty & Bandyopadhyay, Scientific Reports 2013) on two sequences
	# to find an optimal global alignment.
	#Input: S1 and S2 are two strings containing the sequences to be compared. M, Ms, and G are ints indicating
	# match score, mismatch score, and gap score respectively
	import copy
	try:
		import Queue as Q  # ver. < 3.0
	except ImportError:
		import queue as Q
    #create copies of S1 and S2 with a space at the beginning (the root node)
	S1_copy = ' ' + S1
	S2_copy = ' ' + S2
	tree = dict()
	treeoptimal = dict()

	#create priority queue
	q = Q.PriorityQueue()
	
	#create the root node
	#pointers, begin at gap
	P1 = 0
	P2 = 0
	#(-Tmin,Tmax) = (0,0), position 0,0, type 0 (unassigned), PrS=0, list of dict values of children
	Fmin,Fmax = findF(len(S2_copy)-1,len(S1_copy)-1,M,Ms,G)
	tree[(0,0)] = [(-Fmin,Fmax),0,0,0,0,[]]
	print 'root node is' + str(tree[(0,0)])
	#set the optimal score to the Tmin of the root node
	optimal = -tree[(0,0)][0][0]
	m = len(S1_copy) - 1
	n = len(S2_copy) - 1
	#counter = 14
	#check that sequences are both nonempty
	if m != 0 and n != 0:
		#traverse a branch by choosing the best child at each node
		
		while P1 <= m-1 or P2 <= n-1: 
			#choose best child from remaining children of node according to Tmax
			restartloop = False
			possiblechildren =[]
			#expansion type 1. advance both pointers by 1
			#TODO: eliminate children that are not worthwhile
			try:
				#if it's a match
				if S1_copy[P1+1] == S2_copy[P2+1]: prs = M
				#if it's a mismatch
				elif S1_copy[P1+1] != S2_copy[P2+1]: prs = Ms
				#compute future score
				x1 = m-(P1+1)
				x2 = n-(P2+1)

				Fmin,Fmax = findF(x1,x2,M,Ms,G)
				#get present score (sum of all m/ms/g scores that have been encountered so far start from root)
				#present score of current node added to present score of previous node

				PrS = prs+tree[(P1,P2)][4]
				#get fitness score
				Tmin = PrS+Fmin

				Tmax = PrS+Fmax

				#add to possiblechildren list
				possiblechildren.append([(-Tmin,Tmax),P1+1,P2+1,1,PrS,[]])

			except IndexError:
				print 'not possible to advance both'
				continue
			#expansion type 3. keep P1 fixed (gap in S1)
			x1 = m-P1
			x2 = n-(P2+1)
		
			PrS = G+tree[(P1,P2)][4]
	
			Fmin,Fmax = findF(x1,x2,M,Ms,G)

			Tmin = PrS+Fmin
			Tmax = PrS+Fmax
			possiblechildren.append([(-Tmin,Tmax),P1,P2+1,3,PrS,[]])


			#expansion type 2. keep P2 fixed (gap in S2)
			x1 = m-(P1+1)
			x2 = n-P2

			PrS = G+tree[(P1,P2)][4]
			Fmin,Fmax = findF(x1,x2,M,Ms,G)

			Tmin = PrS+Fmin
			Tmax = PrS+Fmax
			possiblechildren.append([(-Tmin,Tmax),P1+1,P2,2,PrS,[]])

			#find the best child according to Tmax

			sortedpossiblechildren = sorted(possiblechildren, key=lambda child: (child[0][1], -child[0][0]))
			print sortedpossiblechildren
			fittestchild = sortedpossiblechildren.pop()
			print 'the fittest child is' + str(fittestchild)
			#add the rest to the priority queue 
			for c in sortedpossiblechildren:
				q.put((-c[0][1],c))
			#cell chosen based on Tmax
			#collisions resolved by Tmin
			#default for PriorityQueue is to return highest priority item
			#if more than one item has the same priority, sorted according to Tmin 
			#(due to chosen ordering in dict values and negative sign on Tmin)
	
			#check to see if this node has been traversed before in an equal or better way
			try:
				if fittestchild[4] <= tree[(fittestchild[1],fittestchild[2])][4] and q.empty() == False:
					print fittestchild[4]
					print tree[(fittestchild[1],fittestchild[2])][4]
					print 'aborting this child'
					#abort this child and start a new child from the priority queue
					
					newnode = pluckNode(q)
					P1 = newnode[1]
					P2 = newnode[2]
					try:
						while newnode[4] <= tree[(P1,P2)] and q.empty() == False:
							print 'aborting this new node : (' + str(P1) + ',' + str(P2) + ')'
							newnode = pluckNode(q)
							P1 = newnode[1]
							P2 = newnode[2]
					except KeyError:
						pass
					if q.empty() == False:
						
						tree[(newnode[1],newnode[2])] = newnode
						

						print 'new node is (' + str(P1) + ',' + str(P2) + ')'
						print newnode
						#add child to appropriate parent
						if newnode[3] == 1:
							#if child is of type 1
							tree[(newnode[1]-1,newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
						elif newnode[3] == 2:
							#type 2
							tree[(newnode[1]-1,newnode[2])][5].append(tree[(newnode[1],newnode[2])])
						elif newnode[3] == 3:
							#type 3
							tree[(newnode[1],newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
						restartloop = True
				else: 
					if q.empty() == False:
					#replace this node with the better node
						tree[(fittestchild[1],fittestchild[2])] = fittestchild
					#add child to child list of parent
						tree[(P1,P2)][5].append(tree[(fittestchild[1],fittestchild[2])])
						#choose fittest child
						P1 = fittestchild[1]
						P2 = fittestchild[2]
			except KeyError:
				#create new node
				tree[(fittestchild[1],fittestchild[2])] = fittestchild

				#add child to child list of parent
				tree[(P1,P2)][5].append(tree[(fittestchild[1],fittestchild[2])])
				P1 = fittestchild[1]
				P2 = fittestchild[2]
			if restartloop: continue
			
			
			
			
			print 'traversing to:'
			print (P1,P2)
			if P1 > m-1 or P2 > n-1 and q.empty() == False:
				if tree[(P1,P2)][0][1] >= optimal:
					print 'set new optimal'
					optimal = tree[(P1,P2)][0][1]
					bestpathnode = tree[(P1,P2)]
				print 'plucked a new node'
				newnode = pluckNode(q)
				P1 = newnode[1]
				P2 = newnode[2]
				print P1,P2
				try:
					while newnode[4] <= tree[(P1,P2)] and q.empty() == False:
						print 'aborting this new node : (' + str(P1) + ',' + str(P2) + ')'
						newnode = pluckNode(q)
						P1 = newnode[1]
						P2 = newnode[2]

				except KeyError:
					pass
				if q.empty() == False:
					tree[(P1,P2)] = newnode
					#add child to appropriate parent
					if newnode[3] == 1:
						#if child is of type 1
						tree[(newnode[1]-1,newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
					elif newnode[3] == 2:
						#type 2
						tree[(newnode[1]-1,newnode[2])][5].append(tree[(newnode[1],newnode[2])])
					elif newnode[3] == 3:
						#type 3
						tree[(newnode[1],newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
			if tree[(P1,P2)][0][1] <= optimal and q.empty() == False:
				newnode = pluckNode(q)
				P1 = newnode[1]
				P2 = newnode[2]
				while newnode[0][1] <= optimal and q.empty()==False:
					newnode = pluckNode(q)
					P1 = newnode[1]
					P2 = newnode[2]
					if newnode[0][1] > optimal:
						tree[(P1,P2)] = newnode
						#add child to appropriate parent
						if newnode[3] == 1:
						#if child is of type 1
							tree[(newnode[1]-1,newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
						elif newnode[3] == 2:
						#type 2
							tree[(newnode[1]-1,newnode[2])][5].append(tree[(newnode[1],newnode[2])])
						elif newnode[3] == 3:
						#type 3
							tree[(newnode[1],newnode[2]-1)][5].append(tree[(newnode[1],newnode[2])])
			if q.empty() == True:
				return tree,bestpathnode
			#counter -=1
			#if counter == 0:
			#	return tree,bestpathnode
			
		#set the lower bound of Tmax to be score of this branch since we must obtain
		#more than this to find a better branch
		# if prune == False and fittestchild[0][1] >= optimal: 
		# 	optimal = fittestchild[0][1]
		# 	bestbranch = currentbranch

	else:
		print 'one or both sequences is empty'
		return None	
	return tree,bestpathnode
			


#Needleman-Wunsch (NW)