def formatSeq(dnasequence):
	#takes in a DNA sequence and converts it to all capital letters, removes spaces, 
	#and throws an exception when characters other than ATCG are present.
	#Input: dnasequence is a string containing the sequence of interest. 
	#Output: copy of dnasequence with all caps and no spaces.

	
	#remove all spaces and convert to uppercase
	editedsequence = dnasequence.replace(' ','')
	editedsequence = editedsequence.upper()

	#check for unwanted characters
	for char in editedsequence:
		if char not in ['A','T','C','G']:
			raise ValueError('Sequence contains non-nucleotide characters. Please remove any characters other than A,T,C,G.')
	return editedsequence

def findF(x1,x2,M,Ms,G):
	if x2 < x1:
		Fmin = x2*Ms+G*(x1-x2)
		Fmax = x2*M+G*(x1-x2)
	else:
		Fmin = x1*Ms+G*(x2-x1)
		Fmax = x1*M+G*(x2-x1)
	return Fmin,Fmax

def FOGSAA(S1,S2):
	#Implements FOGSAA (Chakraborty & Bandyopadhyay, Scientific Reports 2013) on two sequences
	# to find an optimal global alignment.
	#Input: S1 and S2 are two strings containing the sequences to be compared.

	try:
		import Queue as Q  # ver. < 3.0
	except ImportError:
		import queue as Q
     #create copies of S1 and S2 with a space at the beginning
	S1_copy = ' ' + S1
	S2_copy = ' ' + S2

	#create node dictionary
	nodes = dict()
	#create priority queue
	q = Q.PriorityQueue()
	#match score
	M = 1
	#mismatch score
	Ms = -1
	#gap penalty
	G = -2
	#pointers, begin at gap
	P1 = 0
	P2 = 0
	#start root node with (-Tmin,Tmax) = (0,0), position 0,0, type 0 (unassigned), PrS=0
	x1 = len(S2_copy)-1
	x2 = len(S1_copy)-1
	Fmin,Fmax = findF(x1,x2,M,Ms,G)
	nodes[(0,0)] = ((-Fmin,Fmax),(0,0),0,0)
	lowerBound = -10
	if len(S1_copy) != 1 and len(S2_copy) != 1:
		#terminate when Tmax of the top node of the priority queue indicates less than 30% similarity
		#between sequences.
		while lowerBound < -2:
			#if we aren't starting from the root node, then choose the best node from 
			#the priority queue as long as we can potentially get a higher Tmax from it

			if P1 != 0 or P2 != 0:
				newnode = q.get()
				newnode = newnode[1]
				print newnode
				print newnode[0][1]
				print lowerBound

				if newnode[0][1] >=lowerBound:
					P1 = newnode[1][0]
					P2 = newnode[1][1]
					print P1,P2
				else:
					print "not worth pursuing this branch"
			prune = False
			#traverse a branch by choosing the best child at each node
			while P1 <= len(S1_copy)-2 or P2 <= len(S2_copy)-2:
				#check three possible children of node.
				#advance both by 1
				fittestchild = []
				
				try:
					if S1_copy[P1+1] == S2_copy[P2+1]:
						prs = M
					elif S1_copy[P1+1] != S2_copy[P2+1]:
						prs = Ms
					x1 = len(S2_copy)-1-(P2+1)
					x2 = len(S1_copy)-1-(P1+1)
					Fmin,Fmax = findF(x1,x2,M,Ms,G)
					PrS = prs+nodes[(P1,P2)][3]
					Tmin = PrS+Fmin
					Tmax = PrS+Fmax
					#check to see if there is an already expanded node with the same or better Tmax. if 
					#so, don't use.
					try:
						if nodes[(P1+1,P2+1)][0][1] >= PrS:
							#don't add this to the node dict but save for comparision to other children
							fittestchild = ((-Tmin,Tmax),(P1+1,P2+1),1,PrS)
							nobetter = True
						else:
							nodes[(P1+1,P2+1)] = ((-Tmin,Tmax),(P1+1,P2+1),1,PrS)
							fittestchild = nodes[(P1+1,P2+1)]
							nobetter = False
					except KeyError:

						nodes[(P1+1,P2+1)] = ((-Tmin,Tmax),(P1+1,P2+1),1,PrS)
						fittestchild = nodes[(P1+1,P2+1)]
						nobetter = False
				except IndexError:
					print 'not possible to advance both'
					continue
				#keep P2 fixed (gap in S2)
				x1 = len(S2_copy)-1-(P2)
				x2 = len(S1_copy)-1-(P1+1)
				PrS = G+nodes[(P1,P2)][3]
				Fmin,Fmax = findF(x1,x2,M,Ms,G)
				Tmin = PrS+Fmin
				Tmax = PrS+Fmax
				#check to see if there is an already expanded node with the same or better Tmax. if 
					#so, don't use.
				try:
					if nodes[(P1+1,P2)][0][1] >= PrS:
						#don't add this to the node dict but save for comparison to other children
						if fittestchild[0][1] < Tmax:
							fittestchild = ((-Tmin,Tmax),(P1+1,P2),2,PrS)
							nobetter = True
					else:
						nodes[(P1+1,P2)] = ((-Tmin,Tmax),(P1+1,P2),2,PrS)
						if fittestchild[0][1] < Tmax:
							fittestchild = nodes[(P1+1,P2)]
							nobetter = False
				except KeyError:
					nodes[(P1+1,P2)] = ((-Tmin,Tmax),(P1+1,P2),2,PrS)
					if fittestchild[0][1] < Tmax:
						fittestchild = nodes[(P1+1,P2)]
						nobetter = False

				
				
				
				#keep P1 fixed (gap in S1)
				x1 = len(S2_copy)-1-(P2+1)
				x2 = len(S1_copy)-1-(P1)
				PrS = G+nodes[(P1,P2)][3]
				Fmin,Fmax = findF(x1,x2,M,Ms,G)
				Tmin = PrS+Fmin
				Tmax = PrS+Fmax
				#check to see if there is an already expanded node with the same or better Tmax. if 
					#so, don't use.
				try:
					if nodes[(P1,P2+1)][0][1] >= PrS:
						#don't add this to the node dict but save for comparision to other children
						if fittestchild[0][1] < Tmax:
							fittestchild = ((-Tmin,Tmax),(P1,P2+1),3,PrS)
							nobetter = True
					else:
						nodes[(P1,P2+1)] = ((-Tmin,Tmax),(P1,P2+1),3,PrS)
						if fittestchild[0][1] < Tmax:
							fittestchild = nodes[(P1,P2+1)]
							nobetter = False
				except KeyError:
					nodes[(P1,P2+1)] = ((-Tmin,Tmax),(P1,P2+1),3,PrS)
					if fittestchild[0][1] < Tmax:
						fittestchild = nodes[(P1,P2+1)]
						nobetter = False
				
				if nobetter:
					print fittestchild
					print nodes[(fittestchild[1][0],fittestchild[1][1])]
					print 'terminating branch at' + str((fittestchild[1][0],fittestchild[1][1]))
					prune = True
					break
					
				#save less fit to priority queue
				#default for PriorityQueue is to return highest priority item
				#if more than one item has the same priority, sorted according to Tmin (due to chosen ordering in dict values and negative sign on Tmin)
				for n in [nodes[(P1+1,P2+1)],nodes[(P1+1,P2)],nodes[(P1,P2+1)]]:
					if n != fittestchild:
						q.put((-n[0][1],n))
				
				#choose fittest child
				P1 = fittestchild[1][0]
				P2 = fittestchild[1][1]
				print 'traversing to:'
				print (P1,P2)
				print fittestchild
				print nodes[(fittestchild[1][0],fittestchild[1][1])]
				#DELETE LATER; FOR TESTING
				#return fittestchild
			
			#set the lower bound of Tmax to be score of this branch since we must obtain
			#more than this to find a better branch
			if prune == False: 
				lowerBound = fittestchild[0][1]
	return 'lower bound is' + str(lowerBound)			
			







#Needleman-Wunsch (NW)