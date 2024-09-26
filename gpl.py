#### python translation from https://github.com/trgao10/PuenteAlignment/blob/master/code/gplmk.m
import numpy as np
import using_kdtree
from auto3dgm.mesh.mesh import Mesh
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import KDTree

class GPL():
     def surfacearea(mesh):  
         v=mesh.vertices
         f=mesh.faces
         L1 = np.sqrt(np.sum((v[f[:,2],:]-v[f[:,3],:]).^2,2))
         L2 = np.sqrt(np.sum((v[f[:,1],:]-v[f[:,3],:]).^2,2))
         L3 = np.sqrt(np.sum((v[f[:,1],:]-v[f[:,2],:]).^2,2))
         S=(L1+L2+L3)/2
         TriArea=np.sqrt(np.absolute(S.*(S-L1).*(S-L2).*(S-L3)))# semiperimeter formula for the area of a triangle whose sides have lengths L1,L2,L3
         Area=np.sum(TriArea) ########## summing the area 
         return TriArea, Area

     def findpointnormals(points,numnbrs=None, viewPoint = None,dirLargest=None):
#FINDPOINTNORMALS Estimates the normals of a sparse set of n 3d points by
# using a set of the closest neighbours to approximate a plane.
#
#   Required Inputs:
#   points- nx3 set of 3d points (x,y,z)
#
#   Optional Inputs: (will give default values on empty array [])
#   numNeighbours- number of neighbouring points to use in plane fitting
#       (default 9)
#   viewPoint- location all normals will point towards (default [0,0,0])
#   dirLargest- use only the largest component of the normal in determining
#       its direction wrt the viewPoint (generally provides a more stable
#       estimation of planes near the viewPoint, default true)
#
#   Outputs:
#   normals- nx3 set of normals (nx,ny,nz)
#   curvature- nx1 set giving the curvature
#
#   References-
#   The implementation closely follows the method given at
#   http://pointclouds.org/documentation/tutorials/normal_estimation.php
#   This code was used in generating the results for the journal paper
#   Multi-modal sensor calibration using a gradient orientation measure 
#   http://www.zjtaylor.com/welcome/download_pdf?pdf=JFR2013.pdf
#
#   This code was written by Zachary Taylor
#   zacharyjeremytaylor@gmail.com
#   http://www.zjtaylor.com

## check inputs
         if v.shape[0]!=3:
	    raise ValueError('unacceptable input')
         if numnbrs is None:
            numnbrs=9
         if isinstance(numnbrs,int) and numnbrs>0:
            print ("Acceptable input")
            if numnbrs >100:
              print("Too many neighbouring points will be used in plane, expect long run times, large ram usage and poor results near edges") 
         else:
	    raise ValueError('unacceptable input') 



         if viewPoint is None:
            viewPoint=[[0,0,0]]
         if np.size(a,0)==1 and np.size(a,1)==3:
            print ("Acceptable input")
         else
	    raise ValueError('unacceptable input') 
         if dirLargest is None:
            dirLargest=True



         #ensure inputs of correct type
         points = np.float32(points)
         viewPoint = np.float32(viewPoint)
         #create kdtree

         kdtreeobj = using_kdtree.KDTreeSearcher.kdtree(mesh)
         neigh = KNeighborsClassifier(n_neighbors=numnbrs+1) 
         neigh.fit(kdtreeobj) 
         #get nearest neighbours
         n=neigh.kneighbors(points, return_distance=False)       
         # remove self  
         n=np.delete(n,1,axis = 0) 
         #find difference in position from neighbouring points
         p = np.matlib.repmat(points[:,1:3], numnbrs,1)-points[np.reshape(n, (np.size(n), 1)),1:3]
         p = reshape(p,  np.size(points,0),numnbrs,3)
         #calculate values for covariance matrix
         C=np.zeros((np.size(points,0), 6))
         C[:,0] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,1]),axis=1)
         C[:,1] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,2]),axis=1)
         C[:,2] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,3]),axis=1)
         C[:,3] = np.sum(np.array(p[:,:,2])*np.array(p[:,:,2]),axis=1)
         C[:,4] = np.sum(np.array(p[:,:,2])*np.array(p[:,:,3]),axis=1)
         C[:,5] = np.sum(np.array(p[:,:,3])*np.array(p[:,:,3]),axis=1)
         C = C/ numnbrs
         #normals and curvature calculation
         normals = np.zeros(np.size(points))
         curvature = np.zeros(np.size(points,1),1)
         for i in range(np.size(points,1)+1):
         #form covariance matrix
            Cmat = [[C[i,1], C[i,2], C[i,3]],[C[i,2], C[i,4], C[i,5]],[C[i,3], C[i,5], C[i,6]]] 
         #get eigen values and vectors
            v, d = LA.eig(Cmat)
            d = np.diag(d)
            lambda1=np.minimum(d)
            k=np.unravel_index(np.argmin(d, axis=None), d.shape)
         #store normals
            normals[i,:] = np.transpose(v[:,k])
         #store curvature
            curvature[i] = lambda1 / np.sum(d)
         #ensure normals point towards viewPoint
            points = points - np.matlib.repmat(viewPoint,np.size(points,1),1)

     ####this is the main function ########
     def gplmk(mesh , N, seed=empty([0,0])):
# Subsample N points from V, using Gaussian Process landmarking
# Arguments:
#   seed - 3 x M matrix with points to be respected in V, i.e. points that
#   belong to V and that should be the first M points in the resulting X.
#   E.g. when you had a previously subsampled set and you want to just
#   increase the number of sampled points
	 V = mesh.vertices
	 F = mesh.faces
         nV = np.size(V,2)
         Center = np.mean(V,axis=1)
         V = V - np.matlib.repmat(Center,1,nV)
         Area = surfacearea(mesh) #### uses surfacearea  
         V = V*np.sqrt(1/Area)
         curvature=findPointNormals(np.transpose(V),10) 
         Lambda = curvature/np.sum(curvature)
         I=np.transpose(np.array([[F[1,:],F[2,:],F[3,:]]]))
         J=np.transpose(np.array([[F[2,:],F[3,:],F[1,:]]]))
         E=np.array([[I],[J]])
         E = [I  J]
         E = np.transpose(E.sort(axis=0))
         E = np.sort(E, axis=0)   
         E1 = E.ravel().view(np.dtype((np.void, E.dtype.itemsize*E.shape[1])))
         E = np.unique(E1, return_index=True)
         E = E[np.sort(unique_idx)]
         EdgeIdxI = E[1,:]
         EdgeIdxJ = E[2,:]
         bandwidth = np.mean(np.sqrt(np.sum(np.square(V[:,EdgeIdxI]-V[:,EdgeIdxJ]))))/5
         BNN = np.minimum(500,nV)
         nbrs = NearestNeighbors(n_neighbors=BNN+1,algorithm='ball_tree').fit(np.transpose(V))
         q_points=  np.take(np.transpose(V),0:nV)                       ## query points
         distances, indices = nbrs.kneighbors(q_points)     
         fullphi = sparse.coo_matrix((exp(-dist.^2/bandwidth),(np.repmat(1:nV,1,BNN+1),idx)),shape=(nV,nV))
         fullMatProd = fullPhi * sparse.coo_matrix(Lambda,(1:nV,1:nV),shape=(nV,nV))* fullPhi
         KernelTrace = np.diag(fullMatProd)
         ind_seed = knnsearch( np.transpose(V), np.transpose(seed) )########################
         nbrs_seed = NearestNeighbors(n_neighbors=1,algorithm='ball_tree').fit(np.transpose(V))
         distseed,indseed=nbrs.kneighbors(np.transpose(seed)) 
         if linalg.norm(seed - V[:,ind_seed], 'fro')> 1e-10:
            raise ValueError("Some seed point did not belong to the set of points to subsample from")

         n_seed = ind_seed.size
         ind=np.column_stack((np.reshape(ind_seed, (1, n_seed))  np.zeros((1, N - n_seed))))
         invKn = np.zeros(N,N)
         invKn[1:ind[1:n_seed],1:ind[1:n_seed]] = inv[fullMatProd[ind[1:n_seed],ind[1:n_seed]]] 


         for j in range(n_seed + 1 : N+1):
                 if j == 2:
                    invKn[1:(j-1),1:(j-1)] = 1/fullMatProd[ind[1],ind[1]]
                    ptuq = KernelTrace -np.sum(np.transpose(fullMatProd[:,ind[1:j]])*(invKn[1:(j-1),1:(j-1)]*np.transpose(fullMatProd[ind[1:(j-1)],:]),1)), axis=0)  ##########
                 else
                    p = fullMatProd[ind(1:[j-2]],ind[j-1]]
                    mu = 1/(fullMatProd[ind[j-1],ind[j-1]]-np.transpose(p)*invKn[1:(j-2),1:(j-2)]*p)
                    invKn[1:(j-2),1:(j-1)] = invKn[1:(j-2),1:(j-2)]*(np.ones(j-2)+mu*(p*np.transpose(p))*invKn[1:(j-2),1:(j-2)],-mu*p)
                    invKn[j-1,1:(j-1)] =np.column_stack((np.transpose(invKn[1:(j-2),j-1]),mu)) 
                    productEntity = invKn[1:(j-1),1:(j-1)]*fullMatProd[ind[1:(j-1)],:]
                    ptuq = KernelTrace - np.sum(np.transpose(fullMatProd[:,ind[0:j-1]])*productEntity,axis=0)
                 maxUQIdx=np.argmax(ptuq, axis=0)
                 ind[j] = maxUQIdx




