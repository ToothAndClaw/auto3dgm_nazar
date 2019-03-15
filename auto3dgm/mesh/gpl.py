################33 need to check this once more for accuracy


import numpy as np
import using_kdtree
from sklearn.neighbors import NearestNeighbors

class GPL():



     def surfacearea(mesh):
         v=mesh.vertices
         f=mesh.faces
         L1 = np.sqrt(np.sum((v[f[:,2],:]-v[f[:,3],:]).^2,2));
         L2 = np.sqrt(np.sum((v[f[:,1],:]-v[f[:,3],:]).^2,2));
         L3 = np.sqrt(np.sum((v[f[:,1],:]-v[f[:,2],:]).^2,2));
         S=(L1+L2+L3)/2;
         TriArea=np.sqrt(np.absolute(S.*(S-L1).*(S-L2).*(S-L3)));
         Area=np.sum(TriArea);

         return TriArea, Area






     def findpointnormals(points,numnbrs=9, viewPoint = [[0,0,0]]):
         if v.shape[0]!=3:
	    raise ValueError('unacceptable input')
         points = np.float32(points);
         viewPoint = np.float32(viewPoint);
         kdtreeobj = using_kdtree.KDTreeSearcher.kdtree(mesh)
         neigh = KNeighborsClassifier(n_neighbors=numnbrs+1) 
         neigh.fit(kdtreeobj) 
         n=neigh.kneighbors(points, return_distance=False)
         n=np.delete(n,1,axis = 0) 
         p = np.matlib.repmat(points[:,1:3], numnbrs,1)-points[np.reshape(n, (np.size(n), 1)),1:3]
         p = reshape(p,  np.size(points,0),numnbrs,3);
         C=np.zeros((np.size(points,0), 6))
         C[:,0] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,1]),axis=1);
         C[:,1] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,2]),axis=1);
         C[:,2] = np.sum(np.array(p[:,:,1])*np.array(p[:,:,3]),axis=1);
         C[:,3] = np.sum(np.array(p[:,:,2])*np.array(p[:,:,2]),axis=1);
         C[:,4] = np.sum(np.array(p[:,:,2])*np.array(p[:,:,3]),axis=1);
         C[:,5] = np.sum(np.array(p[:,:,3])*np.array(p[:,:,3]),axis=1);

         C = C/ numnbrs;

         normals = np.zeros(np.size(points));
         curvature = np.zeros(np.size(points,1),1);
         for i in range(np.size(points,1)+1):
    
            Cmat = [[C[i,1], C[i,2], C[i,3]],[C[i,2], C[i,4], C[i,5]],[C[i,3], C[i,5], C[i,6]]]
    

            v, d = LA.eig(Cmat)
            d = np.diag(d);
            lambda=np.minimum(d)
            k=np.unravel_index(np.argmin(d, axis=None), d.shape)

            normals[i,:] = np.transpose(v[:,k]);
            curvature(i) = lambda / np.sum(d);

            points = points - np.matlib.repmat(viewPoint,np.size(points,1),1);






     def gplmk( V, F, N, seed):
         nV = np.size(V,2)
         Center = np.mean(V,axis=1)
         V = V - np.matlib.repmat(Center,1,nV)
         Area = surfacearea(np.transpose(V),np.transpose(F));   
         V = V*np.sqrt(1/Area)
         curvature=findPointNormals(np.transpose(V),10); 
         Lambda = curvature/np.sum(curvature);
         I=np.transpose(np.array([[F[1,:],F[2,:],F[3,:]]]))
         J=np.transpose(np.array([[F[2,:],F[3,:],F[1,:]]]))
         E=np.array([[I],[J]])
         E = [I ; J];
         E = np.transpose(E.sort(axis=0))
         E = np.sort(E, axis=0)  ; 
         E1 = E.ravel().view(np.dtype((np.void, E.dtype.itemsize*E.shape[1])))
         E = np.unique(E1, return_index=True)
         E = E[np.sort(unique_idx)]
         EdgeIdxI = E[1,:];
         EdgeIdxJ = E[2,:];
         bandwidth = np.mean(np.sqrt(np.sum(np.square(V[:,EdgeIdxI]-V[:,EdgeIdxJ]))))/5;
         BNN = np.minimum(500,nV)
         atria = nn_prepare(np.transpose(V)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         fullPhi = sparse(np.repmat(1:nV,1,BNN+1),idx,exp(-dist.^2/bandwidth),nV,nV);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         fullMatProd = fullPhi * sparse(1:nV,1:nV,Lambda,nV,nV) * fullPhi; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         KernelTrace = np.diag(fullMatProd);


         ind_seed = knnsearch( np.transpose(V), np.transpose(seed) );%%%%%%%%%%%%%%%%%%%%%%%%
         if linalg.norm(seed - V[:,ind_seed], 'fro')> 1e-10:
            raise ValueError("Some seed point did not belong to the set of points to subsample from")

         n_seed = ind_seed.size
         ind=np.column_stack((np.reshape(ind_seed, (1, n_seed))  np.zeros((1, N - n_seed))))
         invKn = np.zeros(N,N);
         invKn[1:ind[1:n_seed],1:ind[1:n_seed]] = inv[fullMatProd[ind[1:n_seed],ind[1:n_seed]]];  %%%% 





         for j in range(n_seed + 1 : N+1):

                 if j == 2:
                    invKn[1:(j-1),1:(j-1)] = 1/fullMatProd[ind[1],ind[1]];
                    ptuq = KernelTrace - sum(fullMatProd(:,ind(1:(j-1)))'...
                    .*(invKn(1:(j-1),1:(j-1))*fullMatProd(ind(1:(j-1)),:)),1)'############### reqires correction
                 else
                    p = fullMatProd[ind(1:[j-2]],ind[j-1]];
                    mu = 1/(fullMatProd[ind[j-1],ind[j-1]]-np.transpose(p)*invKn[1:(j-2),1:(j-2)]*p);
                    invKn[1:(j-2),1:(j-1)] = invKn[1:(j-2),1:(j-2)]*(np.ones(j-2)+mu*(p*np.transpose(p))*invKn[1:(j-2),1:(j-2)],-mu*p);
                    invKn[j-1,1:(j-1)] =np.column_stack((np.transpose(invKn[1:(j-2),j-1]),mu)) 
                    productEntity = invKn[1:(j-1),1:(j-1)]*fullMatProd[ind[1:(j-1)],:]
                    ptuq = KernelTrace - sum(fullMatProd(:,ind(1:(j-1)))'...
                    .*productEntity,1)'############### reqires correction

                 maxUQIdx=np.argmax(ptuq, axis=0)
                 ind[j] = maxUQIdx;


