def boltz_E(temp, r, s):
	u = partfunc_E(temp)
	relnrs = 1. / u[r - 1] * np.exp(-(s - 1) / (keV* temp))
	return relnrs
