from auto3dgm.mesh.meshfactory import MeshFactory
from numpy import all, amin, any, argmax, array, isclose, ndarray, where, empty
from scipy.spatial.distance import cdist

import random

class Subsample(object):
	"""Class stub, but subsample is for producing subsampled meshes """

	@staticmethod
	def far_point_subsample(mesh, n, seed=empty([0,0])):
		v = mesh.vertices
		if n > v.shape[0] or n < seed.shape[0]:
			raise ValueError('n larger than number of vertices or smaller than number of seed points')
		if isinstance(seed, ndarray) and seed.size:
			if v.shape[1] == 3 and v.ndim == 2:
				# are s in v (or close enough?)
				if all([any(all(isclose(x, v), 1)) for x in seed]):
					# get ind for seed points
					seedint = [where(all(isclose(x, v), axis=1))[0][0] for x in seed]
			else:
				raise ValueError('seed improperly formed, expecting n x 3 array')
		else:
			random.seed()
			seedint = [random.randint(0, v.shape[0]-1)]
		subint = seedint
		for i in range(len(subint),n):
			subint.append(argmax(amin(cdist(v[subint], v), axis=0)))

		return MeshFactory.mesh_from_data(v[subint])