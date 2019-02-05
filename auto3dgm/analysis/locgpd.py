import linalg
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix, KDTree
from scipy.optimize import linear_sum_assignment as Hungary
import numpy as np


class Correspondence:
    # params: self,
    # meshes: either a mesh object or a list or a dictionary of mesh objects
    # meshes: a list of meshes
    # globalize: flag for performing globalized pairwise alignment (default:yes)
    # mirror: flag for allowing mirrored meshes (default: no)
    # reference_index= index of the mesh that other meshes are aligned against
    def __init__(self, meshes, globalize=1, mirror=0, reference_index=0):
        self.globalize = globalize
        self.mirror = mirror
        if isinstance(meshes, (list, dict)):
            if isinstance(meshes, (list,)):
                self.mesh_list = meshes
            if isinstance(meshes, (dict,)):
                self.mesh_list = meshes.keys()
            if (reference_index > len(meshes) - 1):
                msg = 'There are fewer meshes than the given reference_index'
                raise OSError(msg)
            else:
                self.reference_index = reference_index

    @staticmethod
    # An auxiliary method for computing the initial pairwise-alignment
    # Computes the principal components of two meshes and all possible rotations of the 3-axes)
    # params: mesh1, mesh2 meshes that have vertices that are 3 x n matrices
    #        mirror: a flag for whether or not mirror images of the shapes should be considered

    def principal_component_alignment(mesh1, mesh2, mirror):
        X = mesh1.vertices
        Y = mesh2.vertices
        UX, DX, VX = linalg.svd(X, full_matrices=False)
        UY, DY, VY = linalg.svd(Y, full_matrices=False)
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
            cost = distance_matrix(mesh1.T, np.dot(rot, mesh2).T)
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

    
    def locgpd(mesh1, mesh2, R_0=0, M_0, max_iter):

        best_permutation, best_rot = self.best_pairwise_PCA_alignment(mesh1, mesh2)

        if R_0 != 0:
            best_rot = R_0

        V1_sub = mesh1.vertices
        V2_sub = mesh2.vertices

        newV2_sub = np.dot(best_rot.T, V2_sub)

        i = 0
        while True:
            newV2_sub = newV2_sub[:, best_permutation]
            #  Do Kabsch
            cur_rot = Kabsch(newV2_sub.T, V1_sub.T)
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
        gamma = 1.5 * ltwoinf(X - Rotate.dot(Y).dot(Permutate))

        return d, Rotate, Permutate, gamma

    def ltwoinf(X):
        """l2-inf norm of x, i.e.the maximum 12 norm of the columns of x"""
        d = np.sqrt(max(np.square(X).sum(axis=0)))
        return d


    def Kabsch(A, B):
        assert len(A) == len(B)

        N = A.shape[0]  # total points

        centroid_A = mean(A, axis=0)
        centroid_B = mean(B, axis=0)

        # centre the points
        AA = A - np.tile(centroid_A, (N, 1))
        BB = B - np.tile(centroid_B, (N, 1))

        # dot is matrix multiplication for array
        H = np.dot(AA.T, BB)

        U, S, Vt = LA.svd(H)

        R = np.dot(Vt.T, U.T)
        # # special reflection case
        # if linalg.det(R) < 0:
        #     print("Reflection detected")
        #     Vt[-1,:] *= -1
        #     R = Vt.T * U.T
        return R
