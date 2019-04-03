def saha_E(temp, elpress, ionstage):
	#Defining constants that vary with T
	keVT = keV * temp
	kergT = kerg * temp
	eldens = elpress / kergT
	
	u = partfunc_E(temp)
	u = np.append(u, 2) 	# Adding new element to array
	sahaconst = (2. * np.pi * elmass * kergT / (h**2))**1.5 * 2. / eldens
	nstage = np.zeros(5)
	nstage[0] = 1. 	# Setting initial value to 1
	for r in range(4):
		nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-chiion[r] / keVT)
	ntotal = np.sum(nstage)
	nstagerel = nstage / ntotal
	return nstagerel[ionstage - 1]

