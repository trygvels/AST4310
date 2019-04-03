def partfunc_E(temp):
	u = np.zeros(4)	# declare a 4 zero-element array
	for r in range(4):
		for s in range(chiion[r]):
			u[r] = u[r] + np.exp(-s / keV / temp)
	return u # returns all the values of u array

