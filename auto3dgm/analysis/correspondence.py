from numpy import linalg
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix, KDTree
from scipy.optimize import linear_sum_assignment as Hungary
from scipy.sparse import csr_matrix, identity
from scipy.sparse.csgraph import minimum_spanning_tree, shortest_path
import numpy as np
from auto3dgm import jobrun
from auto3dgm.jobrun import jobrun
from auto3dgm.jobrun import job
from auto3dgm.jobrun.jobrun import JobRun
from auto3dgm.jobrun.job import Job


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
        n = len(meshes)

        job_data = self.generate_job_data()
        job_params = self.generate_params()
        job_func = self.generate_func()

        new_job = Job(data=job_data, params=job_params, func=job_func)
        new_jobrun = JobRun(job=new_job)
        output = new_jobrun.execute_jobs()
        #dependent on changing line 52 from jobrun (job_dict['output']['results'].... -> job_dict['output'] = results_dict)
        #also depends on results_dict being in d->float, p->permutation r-> rotaiton mapping
        results_dict = output['output']
        
        self.d_ret = np.array([[None for x in range(n)] for y in range(n)])
        self.p_ret = np.array([[None for x in range(n)] for y in range(n)])
        #what is permutation
        self.r_ret = np.array([[None for x in range(n)] for y in range(n)])


        for key in results_dict.keys():
            #key will be a tuple
            self.d_ret[key[0]][key[1]] = results_dict[key]['d']
            self.p_ret[key[0]][key[1]] = results_dict[key]['p']
            self.r_ret[key[0]][key[1]] = results_dict[key]['r']

        if globalize:
            self.mst_matrix = Correspondence.find_mst(self.d_ret)


    def generate_job_data(self):
        ret = {}
        for first in self.meshes:
            for second in self.meshes:
                if not Correspondence.has_pair(first, second, ret):
                    val = Correspondence.dict_val_gen(self.meshes[first], self.meshes[second], first, second)
                    toopl = (self.meshes[first], self.meshes[second])
                    ret[toopl] = val
        return ret


    @staticmethod
    def dict_val_gen(firstindex, secondindex, first, second):
        return {firstindex: first, secondindex: second}

    @staticmethod
    def has_pair(key1, key2, dictionary):
        str1 = key1 + key2
        str2 = key2 + key1
        if (str1 in dictionary.keys()) or (str2 in dictionary.keys()):
            return True
        return False

    def generate_params(self):
        return {'mirror': self.mirror, 'inital_alignment': self.initial_alignment}

    def generate_func(self):
        return self.initial_alignment

    @staticmethod
    def find_mst(distance_matrix):
        X = csr_matrix([t for t in distance_matrix])
        output = minimum_spanning_tree(X)
        return output.toarray()

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
    def globalize(pa, tree, base, type='mst'):
        '''
        takes a pairwise alignment (pa) formatted as a 3-entry array [D, R, P] and a tree (NP-matrix) and returns
        the global alignment obtained by propagating the tree as a list
        '''
        n = len(tree)
        [r, c] = np.nonzero(tree)
        mm = min(r[0], c[0])
        MM = max(r[0], c[0])
        N = len(pa[2][mm, MM][0])

        retR = [None] * n
        retP = [None] * n
        '''
        BaseDistMatrix = pa[0]
        temp = np.diag(np.diag(BaseDistMatrix))
        BaseDistMatrix = BaseDistMatrix - temp
        BaseDistMatrix += BaseDistMatrix.transpose()

        tD = BaseDistMatrix
        tD += np.diag(np.matrix(np.ones(len(tD)) * np.inf))
        epsilon = np.mean(min(tD, [], 2))

        adjMat = np.exp(-np.square(tD) / np.square(epsilon))

        RCell = [[None] * len(pa[1]) for n in len(pa[1][0])]
        
        for j in range(1, len(RCell)):
            for k in range(j+1, len(RCell[0])):
                RCell[j][k] = pa[1][j][k]
                RCell[k][j] = RCell[j][k].transpose()
        '''
        if type is 'mst':
            for li in range(1, n):
                [dist, pat] = Correspondence.graphshortestpaths(tree+tree.transpose(), li, base)
                P = identity(N)
                R = np.identity(3)
                for lj in range(2, len(pat)):
                    if pat[lj - 1] > pat[lj]:
                        P = np.matmul(P, pa[2][pat[lj], pat[lj-1]])
                        R = np.matmul(pa[1][pat[lj], pat[lj-1]], R)
                    else:
                        P = np.matmul(P, pa[2][pat[lj-1], pat[lj]])
                        R = np.matmul(pa[1][pat[lj-1], pat[lj]], R)
                retR[li] = R
                retP[li] = P
        return [retR, retP]                

    @staticmethod
    def getpath(predecessors, index, reference):
        ancestors = predecessors[reference]
        tracer=index
        ret=[]
        while(tracer != reference):
            ret.append(tracer)
            tracer = ancestors[tracer]
        ret.append(reference)
        ret.reverse()
        return ret

    @staticmethod
    # An auxiliary method for computing the initial pairwise-alignment
    # Computes the principal components of two meshes and all possible rotations of the 3-axes)
    # params: mesh1, mesh2 meshes that have vertices that are 3 x n matrices
    #        mirror: a flag for whether or not mirror images of the shapes should be considered

    def principal_component_alignment(mesh1, mesh2, mirror):
        X = mesh1.vertices.T
        Y = mesh2.vertices.T
        UX, DX, VX = linalg.svd(X, full_matrices=False)
        UY, DY, VY = linalg.svd(Y, full_matrices=False)
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
            R.append(np.dot(UX * i, UY.T))

        return R

    @staticmethod
    # Computes the best possible initial alignment for meshes 1 and 2
    # Mesh 1 is used as the reference
    # params: mesh1, mesh2 meshes that have vertices that are 3 x n matrices
    #        mirror: a flag for whether or not mirror images of the shapes should be considered
    def best_pairwise_PCA_alignment(mesh1, mesh2, self):
        R = self.principal_component_alignment(mesh1, mesh2, self.mirror)
        permutations = []
        for rot, i in zip(R, range(len(R))):
            min_cost = np.ones(len(R)) * np.inf
            cost = distance_matrix(mesh1.vertices, np.dot(rot, mesh2.vertices.T).T)
            # The hungarian algorithm:
            V1_ind, V2_ind = linear_sum_assignment(cost)
            min_cost[i] = np.sqrt(np.sum(cost[V1_ind, V2_ind]))
            permutations.append(V2_ind)

        best_rot_ind = np.argmin(min_cost)
        best_permutation = permutations[best_rot_ind]
        best_rot = R[best_rot_ind]

        return best_permutation, best_rot



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
    def locgpd(self,mesh1, mesh2, R_0, M_0, max_iter):

        best_permutation, best_rot = self.best_pairwise_PCA_alignment(mesh1, mesh2,self)

        if R_0 != 0:
            best_rot = R_0

        V1_sub = mesh1.vertices.T
        V2_sub = mesh2.vertices.T

        newV2_sub = np.dot(best_rot.T, V2_sub)

        i = 0
        while True:
            newV2_sub = newV2_sub[:, best_permutation]
            #  Do Kabsch
            cur_rot = self.Kabsch(newV2_sub.T, V1_sub.T)
            newV2_sub = np.dot(cur_rot.T, newV2_sub)
            print("after Kab cost = ", np.linalg.norm(V1_sub - newV2_sub))
            # Do Hungary
            cur_cost = distance_matrix(V1_sub.T, newV2_sub.T)
            cur_V1_ind, cur_permutation = Hungary(cur_cost)
            # print(cur_permutation)
            print("after Hungary cost = ", np.sqrt(np.sum(cur_cost[cur_V1_ind, cur_permutation] ** 2)))
            if i > max_iter:
                break
            else:
                if i % 100 == 0:
                    # print("Current error is: ", np.sum((cur_permutation - best_permutation) != 0))
                    print("current iteration is:", i)
            # update
            best_permutation = cur_permutation
            best_rot = cur_rot
            i += 1

        d = np.sum((cur_permutation - best_permutation))
        Rotate = best_rot
        Permutate = best_permutation
        gamma = 1.5 * self.ltwoinf(V1_sub - newV2_sub)

        return d, Rotate, Permutate, gamma

    def ltwoinf(self,X):
        """l2-inf norm of x, i.e.the maximum 12 norm of the columns of x"""
        d = np.sqrt(max(np.square(X).sum(axis=0)))
        return d


    def Kabsch(self,A, B):
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

        R = np.dot(Vt.T, U.T)
        # # special reflection case
        # if linalg.det(R) < 0:
        #     print("Reflection detected")
        #     Vt[-1,:] *= -1
        #     R = Vt.T * U.T
        return R
