def buildlat(G, H):
	if H.order()!= 2 or not H.is_normal(G):
		raise ValueError("The subgroup must be normal of order 2")
	GM = GLattice([])
	IL = GM.induced_lattice(G)
	m = IL.action_matrix(H[1])
	EV = m.eigenspaces_right()
	if EV[0][0] != 1:
		raise ValueError("Hmmmm")
	else:
		B = EV[0][1].basis()
	SL1 = IL.sublattice(B)
	SL2 = SL1.zero_sum_sublattice()
	L = IL.quotient_lattice(SL2)
	Lder = IL.quotient_lattice(SL1) 	
	return [L, Lder]

def sagesmallgroup(n, m):
	Gp = libgap.SmallGroup(n, m)
	return PermutationGroup(Gp.IsomorphismPermGroup().Image().GeneratorsOfGroup().sage())

def computor(o):
	n = gap.NrSmallGroups(o).sage()
	for i in range(n):
		print(i+1)
		G = sagesmallgroup(o, i+1)
		SBG = [h for h in G.normal_subgroups() if h.order()==2]
		for h in SBG :
			print(h)
			[L, Lder] = buildlat(G, h)
			print("lat : ")
			print(L.Tate_Cohomology(1))
			print("lat_der : ")
			print(Lder.Tate_Cohomology(1))
			print("------------")
		print("=============")

def isinrad(g, e):
	for i in range(g.order()):
		if g**i == e:
			return g
	return e

def radbuilder(G, e):
	elts = []
	for g in G:
		elts += [isinrad(g, e)]
	return G.subgroup(elts)

def computor2(o):
	n = gap.NrSmallGroups(o).sage()
	for i in range(n):
		print(i+1)
		G = sagesmallgroup(o, i+1)
		gord = G.order()
		SBG = [h for h in G.normal_subgroups() if h.order()==2]
		for h in SBG :
			print(h)
			[L, Lder] = buildlat(G, h)
			radord = radbuilder(G, h[1]).order()
			rat = gord / radord
			if is_odd(rat):
				print("odd")
			else:
				print("even")
			print(L.Tate_Cohomology(1))
			print("------------")
		print("=============")


def cycliccontains(G, e):
	subgs = [h for h in G.subgroups() if (h.is_cyclic() and is_power_of_two(h.order()))]
	size = 2
	winner = G.subgroup([e])
	for h in subgs:
		if e in h and h.order()>size:
			winner = h
			size = h.order()
	return [winner, size]

def computor3(o):
	n = gap.NrSmallGroups(o).sage()
	for i in range(n):
		G = sagesmallgroup(o, i+1)
		gord = G.order()
		SBG = [h for h in G.normal_subgroups() if h.order()==2]
		for h in SBG :
			[L, Lder] = buildlat(G, h)
			[cycle, radord] = cycliccontains(G, h[1])
			rat = gord / radord
			if is_odd(rat) and (L.Tate_Cohomology(1) == [2]):
				return False
			elif is_even(rat) and (L.Tate_Cohomology(1) == []):
				return False
	return True

def computor4(o):
	n = gap.NrSmallGroups(o).sage()
	for i in range(n):
		print(i+1)
		G = sagesmallgroup(o, i+1)
		gord = G.order()
		SBG = [h for h in G.normal_subgroups() if h.order()==2]
		for h in SBG :
			print(h)
			[L, Lder] = buildlat(G, h)
			print(L.Tate_Cohomology(1))
			print("------------")
			print("syltest")
			syl = G.sylow_subgroup(2)
			index = G.order()/syl.order()
			print(syl.is_cyclic())
			if syl.is_cyclic():
				print(len(L.subgroup_lattice(syl).Tate_Cohomology(1)) == index-1)
			else:
				print(len(L.subgroup_lattice(syl).Tate_Cohomology(1)) == index)
			#normsyl = G.normalizer(syl)
			#print("VVVVVVVVVVVV")
			#print(normsyl.is_abelian())
			#print(L.subgroup_lattice(normsyl).Tate_Cohomology(1))
		print("=============")

def buildlat2(G, H):
	if H.order()!= 2 or not H.is_normal(G):
		raise ValueError("The subgroup must be normal of order 2")
	GM = GLattice([])
	IL = GM.induced_lattice(G)
	m = IL.action_matrix(H[1])
	EV = m.eigenspaces_right()
	if EV[0][0] != 1:
		raise ValueError("Hmmmm")
	else:
		B = EV[0][1].basis()
	SL1 = IL.sublattice(B)
	SL2 = SL1.zero_sum_sublattice()
	L = IL.quotient_lattice(SL2)
	Lder = IL.quotient_lattice(SL1) 	
	return [L, Lder, SL1, SL2]


def another_test(o):
	print("---------")
	print(o)
	print("---------")
	n = gap.NrSmallGroups(o).sage()
	for i in range(n):
		print(i+1)
		G = sagesmallgroup(o, i+1)
		SYL = G.sylow_subgroup(2)
		index = G.order()/SYL.order()
		norms = [h.order() for h in G.normal_subgroups()]
		if 2 in norms:
			print((index in norms))
		print("=================")


def countord2(group):
	counter = 0
	for g in group:
		if g.order() ==2 :
			counter+=1
	return counter

def first_coboundary_space(l):

    lat = l.isomorphic_ambient_lattice()
    r = lat.rank()
    group = lat.group()
    grouplist = list(group)
    coboundary_basis = []
    if r != 0:
        actionmat = [matrix(lat._action_morphism.Image(g).sage()) for g in grouplist]
        B = lat.basis()
        if not grouplist[0] == group.identity():
            grouplist = [group.identity()] + grouplist.remove(group.identity())
        for b in B:
            l = []
            for i in actionmat:
                l += i*b-b
            coboundary_basis.append(matrix(group.order(), l))
    return coboundary_basis


def returnindex(elt, li):
	for i in range(len(li)):
		if elt in li[i]:
			return i

def checkgroup(group, h):
	cos = group.cosets(h)
	counter = 0 
	for g in group:
		for c in cos:
			elt = g*c[0]
			i = returnindex(elt, cos)
			if elt == cos[i][1] :
				counter += 1
		if mod(counter, 2) :
			return "odd"
		else :
			counter = 0
	return "even" 

def checkgroupall(group):
	subs = [h for h in group.normal_subgroups() if h.order()==2]
	for h in subs:
		print(checkgroup(group, h))
		print("---------------")