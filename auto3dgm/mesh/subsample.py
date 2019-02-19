"""
TODO: Please fix this comment
"Subsample (Class with Instance Attributes and Methods, and Static Methods)"
"Class summary: Subsamples one or more meshes"
"Constructor input: Optional parameters including…Number of points Subsample methodSet of meshes"
"Instance attributes: Number of points Subsample method Set of meshes"
"Instance methods: "
"prepare_analysis "
"Receives one or more meshes, and optionally number of points and subsampling method"
"Sets parameters and meshes as object attributes"
"Calls and returns result of export_analysis"
"export_analysis"
"Returns a data object capable of being input into JobRun. "
"This data object could be something like a dict similar to "
"{‘data’: {‘a’: mesh1, ‘b’: mesh2, ‘c’: mesh3, etc.}, "
" ‘params’: {‘point_number’: 200, ‘subsample_method’: ‘GPR’}, ‘func’: function_reference }"
"process_result"
"Receives a JubRun output data object, and returns the dict of Mesh objects included with that"
"Static methods:"
"Methods specific to each subsampling algorithm, each method receives a mesh "
"and number of points and returns a subsampled mesh"
"""

import vtk 
from numpy import all, amin, any, argmax, array, isclose, ndarray, where, empty
import random
from scipy.spatial.distance import cdist
from auto3dgm.mesh.mesh import Mesh
from auto3dgm.mesh.meshfactory import MeshFactory
import numpy as np

class Subsample:
    def __init__(self,numberofpoints,subsamplemethod,setofmeshes):
        self.pts=numberofpoints
        self.subsmpl=subsamplemethod
        self.meshset=setofmeshes

    def prepare_analysis(self):
        pass
         
    ## class method
    def export_analysis(job_run_object):
        pass

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
        # list of integers that subsampled
        return MeshFactory.mesh_from_data(v[subint])

    # far_point_subsample('mesh') TODO: What is this
     
    # TODO: All of this causes errors because of Matlab syntax
    # @staticmethod
    # def gpl_subsample(mesh, n, seed=empty([0,0])):
    #     v = mesh.vertices
    #     if n > v.shape[0] or n < seed.shape[0]:
    #         raise ValueError('n larger than number of vertices or smaller than number of seed points')
    #     f=mesh.faces
    #     nV=np.size(v,1)
    #     Center = np.mean(v, axis=1)
    #     v = v - np.matlib.repmat(Center, 1, nV)
    #     Area = compute_surface_area(np.transpose(v),np.transpose(f))
    #     v = v * sqrt(1/Area)
    #     curvature = findPointNormals(np.transpose(v),10)
    #     Lambda = curvature/np.sum(curvature)
    #     I = np.array([F[1,:], F[2,:],F[3,:]])
    #     J = np.array([F[2,:], F[3,:],F[1,:]])
    #     E = np.array([[I], [J]])

    #     EdgeIdxI = E[1,:]
    #     EdgeIdxJ = E[2,:]
    #     bandwidth = np.mean(np.sqrt(np.sum((V(:,EdgeIdxI)-V(:,EdgeIdxJ)).^2)))/5

    # @staticmethod
    # def compute_surface_area(mesh):
    #     v = mesh.vertices
    #     f = mesh.faces
    #     L1 = np.sqrt(np.sum((v[F[:,2],:]-V[F[:,3],:]).^2,2));
    #     L2 = np.sqrt(np.sum((v[F[:,1],:]-V[F[:,3],:]).^2,2));
    #     L3 = np.sqrt(np.sum((v[F[:,1],:]-V[F[:,2],:]).^2,2));
    #     S=(L1+L2+L3)/2;
    #     TriArea=np.sqrt(np.absolute(S.*(S-L1).*(S-L2).*(S-L3)));
    #     Area=np.sum(TriArea);

