from numpy import linalg
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix, KDTree
from scipy.optimize import linear_sum_assignment as Hungary
from scipy.sparse import csr_matrix, identity, find, lil_matrix, coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
from scipy import sparse as sp
import numpy as np
from pprint import pprint
from auto3dgm_nazar import jobrun
from auto3dgm_nazar.jobrun import jobrun
from auto3dgm_nazar.jobrun import job
from auto3dgm_nazar.jobrun.jobrun import JobRun
from auto3dgm_nazar.jobrun.job import Job
import os
file_path = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MOSEKLM_LICENSE_FILE', os.path.join(file_path, '../lib/mosek.lic'))
import sys
import mosek
import pdb

class Correspondence:
    #params: self,
    #meshes: either a mesh object or a list (should be list) or a dictionary of mesh objects
    #meshes: a list of meshes 
    #globalize: flag for performing globalized pairwise alignment (default:yes)
    #mirror: flag for allowing mirrored meshes (default: no)
    #reference_index= index of the mesh that other meshes are aligned against
    
    def __init__(self, meshes, globalize=1, mirror=0, initial_alignment=None):
        #assumes that meshes is an ordered list of meshes
        self.globalizeparam=globalize
        self.mirror=mirror
        self.meshes=meshes
        self.initial_alignment=initial_alignment
        if self.initial_alignment:
            self.initial_pa_alignment=self.localize()
        n = len(meshes)
        n_pts = meshes[0].vertices.shape[0]

        job_data = self.generate_job_data()
        job_params = self.generate_params()
        job_func = self.generate_func()

        new_job = Job(data=job_data, params=job_params, func=job_func)
        new_jobrun = JobRun(job=new_job)
        output = new_jobrun.execute_jobs()
        #dependent on changing line 52 from jobrun (job_dict['output']['results'].... -> job_dict['output'] = results_dict)
        #also depends on results_dict being in d->float, p->permutation r-> rotaiton mapping
        self.results_dict = output['output']
        
        self.d_ret = np.empty([n, n])
        #self.p_ret = np.empty([n, n, n_pts], dtype=int)
        self.p_ret = {}
        #what is permutation
        self.r_ret = np.array([[None for x in range(n)] for y in range(n)])


        for key in self.results_dict.keys():
            #key will be a tuple
            self.d_ret[key[0]][key[1]] = self.results_dict[key]['d']

            x = self.results_dict[key]['p']
            if key[0] not in self.p_ret:
                self.p_ret[key[0]] = {}
            self.p_ret[key[0]][key[1]] = self.results_dict[key]['p']

            self.r_ret[key[0]][key[1]] = self.results_dict[key]['r']
        self.pairwise_alignment = {'d': self.d_ret, 'r': self.r_ret, 'p': self.p_ret}    

        print('+++++')
        print(self.d_ret)
        #print(self.p_ret)
        print(self.r_ret)

        if self.globalizeparam:
            self.mst_matrix = Correspondence.find_mst(self.d_ret)
            self.globalized_alignment = self.globalize(self.pairwise_alignment, self.mst_matrix)

    def localize(self):
        if not self.initial_alignment:
            return
        p = {}
        r = np.array([[None for x in range(len(self.initial_alignment['r']))] for y in range(len(self.initial_alignment['r']))])
        for i in range(0, len(self.initial_alignment['r'])):
            for j in range(0, len(self.initial_alignment['r'])):            
                r[i][j] = self.initial_alignment['r'][i].T @ self.initial_alignment['r'][j]
                if i not in p:
                    p[i] = {}
                p_temp = self.initial_alignment['p'][j] @ self.initial_alignment['p'][i].T
                p[i][j] = self.permutation_sparse_to_flat(p_temp)
        return {'r': r, 'p': p}


    def generate_job_data(self):
        ret = {}
        for indexf, first in enumerate(self.meshes):
            for indexs, second in enumerate(self.meshes):
                if not Correspondence.has_pair(indexf, indexs, ret):
                    if self.initial_alignment is None:
                        initial_alignment = None
                    else:
                        initial_alignment = {
                            'r': self.initial_pa_alignment['r'][indexf][indexs],
                            'p': self.initial_pa_alignment['p'][indexf][indexs]
                        }
                    val = Correspondence.dict_val_gen(first, second, initial_alignment=initial_alignment)
                    toopl = (indexf, indexs)
                    ret[toopl] = val
                
        return ret

    def generate_params(self):
            return {'mirror': self.mirror}

    def generate_func(self):
        return self.align


    @staticmethod
    def dict_val_gen(first, second, initial_alignment=None):
        return {'mesh1': first, 'mesh2': second, 'initial_alignment': initial_alignment}

    @staticmethod
    def has_pair(key1, key2, dictionary):
        tooplf = (key1, key2)
        toopls = (key2, key1)
        if (tooplf in dictionary.keys()) or (toopls in dictionary.keys()):
            #for testing, this has been set to False, Typically, should be set to True
            return False
        return False

    @staticmethod
    def find_mst(distance_matrix):
        return minimum_spanning_tree(csr_matrix(distance_matrix)).toarray()

    @staticmethod
    def graphshortestpaths(graph, index, reference):
        '''
        returns a [dist, pat] object with dist (double of distance between index and reference) and path (1xN vector representing path from reference->index)
        '''
        [dist_matrix, predecessors] = shortest_path(graph, return_predecessors=True)
        dist = dist_matrix[reference, index]
        pat = Correspondence.getpath(predecessors, index, reference)
        return [dist, pat]

    @staticmethod
    def globalize(pa, tree, base=0, type='mst'):
        '''
        takes a pairwise alignment (pa) formatted as a 3-entry array [D, R, P] and a tree (NP-matrix) and returns
        the global alignment obtained by propagating the tree as a list
        '''
        n = len(tree)
        [r, c] = np.nonzero(tree)
        mm = min(r[0], c[0])
        MM = max(r[0], c[0])
        pa_r = pa['r']
        pa_p = pa['p']
        N = pa_p[mm][MM].shape[0]

        retR = [None] * n
        retP = [None] * n

        if type is 'mst':
            for li in range(0, n):
                [dist, pat] = Correspondence.graphshortestpaths(tree+tree.transpose(), li, base)
                P = identity(N, dtype=int)
                R = np.identity(3)
                for lj in range(1, len(pat)):
                    #P = P @ Correspondence.permutation_flat_to_sparse(pa_p[pat[lj]][pat[lj-1]])
                    R = pa_r[pat[lj], pat[lj-1]] @ R
                retR[li] = R
                retP[li] = P
        return {'r': retR, 'p': retP}               

    @staticmethod
    def getpath(predecessors, index, reference):
        ancestors = predecessors[reference]
        tracer=index
        ret=[]
        while(tracer != reference):
            ret.append(tracer)
            tracer = ancestors[tracer]
        ret.append(reference)
        #ret.reverse()
        return ret

    @staticmethod
    def permutation_flat_to_sparse(p):
        data = np.ones(p.shape[0], dtype=int)
        cols = np.array(range(0,p.shape[0]), dtype=int)
        return csr_matrix((data, (p, cols)), shape=(p.shape[0], p.shape[0]))

    @staticmethod
    def permutation_sparse_to_flat(p):
        return find(p)[0]

    @staticmethod
    # An auxiliary method for computing the initial pairwise-alignment
    # Computes the principal components of two meshes and all possible rotations of the 3-axes)
    # params: mesh1, mesh2 meshes that have vertices that are 3 x n matrices
    #        mirror: a flag for whether or not mirror images of the shapes should be considered
    def principal_component_alignment(mesh1, mesh2, mirror):
        mirror = 1
        X = mesh1.vertices.T
        Y = mesh2.vertices.T
        #print(X)
        #print(Y)
        UX, DX, VX = linalg.svd(X, full_matrices=False)
        UY, DY, VY = linalg.svd(Y, full_matrices=False)
        #print(UX)
        #print(UY)
        #print(UX.T)
        #print(UY.T)
        P=[]
        R=[]

        P.append(np.array([1, 1, 1]))
        P.append(np.array([1, -1, -1]))
        P.append(np.array([-1, -1, 1]))
        P.append(np.array([-1, 1, -1]))
        if (mirror == 1):
            P.append(np.array([-1, 1, 1]))
            P.append(np.array([1, -1, 1]))
            P.append(np.array([1, 1, -1]))
            P.append(np.array([-1, -1, -1]))
        for i in P:
            R.append(np.dot(np.dot(UX, np.diag(i)), UY.T))
        return R

    @staticmethod
    # Computes the best possible initial alignment for meshes 1 and 2
    # Mesh 1 is used as the reference
    # params: mesh1, mesh2 meshes that have vertices that are 3 x n matrices
    #        mirror: a flag for whether or not mirror images of the shapes should be considered
    def best_pairwise_PCA_alignment(mesh1, mesh2, mirror):
        R = Correspondence.principal_component_alignment(mesh1, mesh2, mirror)
        M_0 = np.ones([len(mesh1.vertices), len(mesh1.vertices)])
        permutations = []
        min_cost = np.ones(len(R)) * np.inf
        for rot, i in zip(R, range(len(R))):
            cost = distance_matrix(mesh1.vertices, np.dot(rot, mesh2.vertices.T).T)**2
            # The hungarian algorithm:
            P, d = Correspondence.linassign(M_0, cost)
            min_cost[i] = d
            print('eight rotation error is', np.sqrt(d))
            permutations.append(P)

        best_rot_ind = np.argmin(min_cost)
        best_permutation = permutations[best_rot_ind]
        best_rot = R[best_rot_ind]
        d = min_cost[best_rot_ind]

        return {'p': best_permutation, 'r': best_rot, 'd': d}

    @staticmethod
    # Returns the meshed aligned by the initial PCA component pairing.
    # Everything is aligned against the mesh specified by the reference_index
    def initial_rotation(self):
        mesh1 = self.mesh_list[self.reference_index]
        Rotations = []
        for i in self.mesh_list:
            Rotations.append(self.best_pairwise_PCA_alignment(mesh1, i, self))
        return Rotations

    @staticmethod
    def permutation_from_rotation(mesh1, mesh2, R):
        cost = distance_matrix(mesh1.vertices, np.dot(R, mesh2.vertices.T).T)**2
        # The hungarian algorithm:
        M_0 = np.ones([len(mesh1.vertices), len(mesh1.vertices)])
        P, d = Correspondence.linassign(M_0, cost)
        return P, d

    @staticmethod
    # Computes the Local Generalized Procrustes Distance between meshes.
    # NOTE: Params R_0, M_0 needs specification.
    def locgpd(mesh1, mesh2, R_0=None, M_0=None, max_iter=1000, mirror=False):
        # best_permutation and best_rot come from PCA
        if M_0 is None:
            M_0 = np.ones([len(mesh1.vertices), len(mesh1.vertices)])

        if R_0 is None:
            best_permutation, best_rot, best_d = Correspondence.best_pairwise_PCA_alignment(mesh1, mesh2, mirror)
        else:
            best_rot = R_0
            best_permutation, best_d = Correspondence.permutation_from_rotation(mesh1, mesh2, R_0)

        
        print(mesh1.name)
        print(mesh2.name)
        #print(best_rot)
        #print(best_permutation)

        V1 = mesh1.vertices.T
        V2 = mesh2.vertices.T
        #newV2 = np.dot(best_rot, V2)
        newV2 = V2

        i = 0
        while True:
            newV2= np.squeeze(np.asarray(newV2 @ best_permutation))
            err1 = V1 - np.dot(best_rot, newV2)
            print("after Hungary", np.linalg.norm(err1))
            
            # Do Kabsch
            cur_rot = Correspondence.Kabsch(V1.T, newV2.T)
            #newV2 = np.dot(cur_rot, newV2)
            err2 = V1 - np.dot(cur_rot, newV2)
            cur_d = np.linalg.norm(err2)
            print("after Kabsch", np.linalg.norm(err2))
            
            # Do Hungary
            #cur_cost = distance_matrix(V1.T, np.dot(cur_rot, newV2).T)**2
            gamma = 1.5 * Correspondence.ltwoinf(V1 - newV2)
            M, MD2 = Correspondence.jrangesearch(V1.T, np.dot(cur_rot, newV2).T, gamma)
            #pdb.set_trace()
            #print(M)
            #cur_trash, cur_permutation, garbage = lap.lapjv(cur_cost)
            cur_permutation, trash = Correspondence.linassign(M, MD2)
            
            #if i > max_iter or cur_trash - best_trash > 0 or ((abs(cur_trash - best_trash)<1e-5*best_trash)):
            if cur_d < 0.00000001 or i > max_iter or abs(best_d-cur_d)<0.00000001:
                break
            else:
                if i % 1 == 0:
                    #print('start')
                    print("Current error is: ", cur_d)
                    # print("current iteration is:", i)
                    # print("current error is: ", np.sqrt(cur_trash))
            # update
            best_permutation = cur_permutation
            best_rot = cur_rot
            best_d = cur_d
            i += 1

        d = np.sqrt(best_d)
        Rotate = best_rot
        #print(Rotate)
        #print(d)
        Permutate = best_permutation
        return {'d': d, 'r': Rotate, 'p': Permutate, 'g': gamma}

    @staticmethod
    def ltwoinf(X):
        """l2-inf norm of x, i.e.the maximum 12 norm of the columns of x"""
        d = np.sqrt(max(np.square(X).sum(axis=0)))
        return d

    @staticmethod
    def Kabsch(A, B):
        assert len(A) == len(B)

        N = A.shape[0]  # total points

        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)

        # centre the points
        AA = A - np.tile(centroid_A, (N, 1))
        BB = B - np.tile(centroid_B, (N, 1))

        # dot is matrix multiplication for array
        H = np.dot(AA.T, BB)
        U, S, Vt = linalg.svd(H)
        R = np.dot(U, Vt)
        return R

    @staticmethod
    def align(mesh1, mesh2, mirror, initial_alignment=None):
        n = len(mesh1.vertices)
        if not initial_alignment:
            initial_alignment = Correspondence.best_pairwise_PCA_alignment(mesh1=mesh1, mesh2=mesh2, mirror=mirror)
        return Correspondence.locgpd(mesh1=mesh1, mesh2=mesh2, R_0=initial_alignment['r'], M_0=np.ones((n, n)), mirror=mirror)

    @staticmethod
    def jrangesearch(X,Y,epsilon):
        I = X.shape[0]
        J = X.shape[0]
        MD2 = np.zeros((I, J))
        #MD2 = [[0 for i in range(I)] for j in range (J)]
        M = np.ones((I,J))
        for i in range(I):
            for j in range(J):
                A = X[i]
                B = Y[j]
                d2 = np.linalg.norm(A-B)**2
                MD2[i][j]=d2
                if (d2>epsilon):
                    M[i,j]=0
        #MD2 = np.asarray(MD2)
        #M = lil_matrix(M)
        return(M,MD2)
    
    '''
    @staticmethod
    def ltwoinf(X):
        Y = np.sqrt(max(np.sum(np.multiply(X,X),axis=0)))
        return Y
    '''

    @staticmethod
    def linassign(A, D):
        N = A.shape[0]
        if np.array_equal(A, np.ones((N, N))):
            rowId, colId = Hungary(D)
            P = sp.coo_matrix((np.ones(N),(rowId,colId)),shape=(N, N))
            P = P.T
            d = D[rowId, colId].sum()
        else:
            #use mosek to speed up
            A = lil_matrix(A)
            tmpD = np.reshape(D.T, (N*N))
            tmpA = np.reshape(A.T, (N*N))
            ivars = np.nonzero(tmpA)
            ivars = ivars[1]
            n_vars = len(ivars)
            #pdb.set_trace()
            
            # mosek enviroment
            inf = 0.0
            
            # Define a stream printer to grab output from MOSEK
            def streamprinter(text):
                sys.stdout.write(text)
                sys.stdout.flush()
            
            numvar = n_vars
            numcon = 2 * N
            
            # Make mosek environment
            with mosek.Env() as env:
                # Create a task object
                with env.Task(0, 0) as task:
                    # Attach a log stream printer to the task
                    #task.set_Stream(mosek.streamtype.log, streamprinter)
                    
                    #build equality constraints for all variables
                    Aeq1 = np.kron(np.eye(N), np.ones(N))
                    Aeq2 = np.kron(np.ones(N), np.eye(N))
                    Aeq = np.concatenate((Aeq1, Aeq2), axis = 0)
                    
                    #reduce the number of variables using ivars
                    Aeq_red = Aeq[:, ivars]
                    Asparse = Aeq_red
                    #Asparse = lil_matrix(Aeq_red)
                    dissimilarity_for_lp_red = tmpD[ivars]
                    c = dissimilarity_for_lp_red
                    
                    #positivity constraints and rhs of equality constraints
                    task.appendcons(numcon)
                    
                    # Append 'numvar' variables.
                    # The variables will initially be fixed at zero (x=0).
                    task.appendvars(numvar)
                    
                    # Input A
                    Arows, Acols, Avals = find(Asparse)
                    task.putaijlist(Arows, Acols, Avals)
                    
                    # Input c
                    ccols = list(range(numvar))
                    task.putclist(ccols, c)
                    
                    # Input constraints
                    task.putconboundlistconst(range(numcon), mosek.boundkey.fx, 1, 1)
                    
                    # Input bounds
                    task.putvarboundlistconst(range(numvar), mosek.boundkey.lo, 0, +inf)
                    
                    # Input the objective sense (minimize/maximize)
                    task.putobjsense(mosek.objsense.minimize)
                    
                    # Input optimizer
                    task.putintparam(mosek.iparam.optimizer, mosek.optimizertype.primal_simplex)
                    
                    # Solve the problem
                    task.optimize()
                    
                    # Print a summary containing information
                    # about the solution for debugging purposes
                    #task.solutionsummary(mosek.streamtype.msg)
                    
                    # Get status information about the solution
                    #solsta = task.getsolsta(mosek.soltype.bas)
                    
                    # Output the xx
                    xx = [0.] * numvar
                    task.getxx(mosek.soltype.bas, # Request the basic solution.
                               xx)
                    
                    d = task.getprimalobj(mosek.soltype.bas)
                    P = sp.coo_matrix((xx,(ivars, np.zeros(n_vars))),shape=(N*N, 1))
                    P = np.reshape(P, (N, N))
                    P = P.T
        return P, d
