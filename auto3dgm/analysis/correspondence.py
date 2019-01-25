class Correspondence:
    #params: self,
    #meshes: either a mesh object or a list or a dictionary of mesh objects
    #meshes: a list of meshes 
    #globalize: flag for performing globalized pairwise alignment (default:yes)
    #mirror: flag for allowing mirrored meshes (default: no)
    #reference_index= index of the mesh that other meshes are aligned against
    def __init__(self, meshes, globalize=1, mirror=0, reference_index=0):
        self.globalize=globalize
        self.mirror=mirror
        if isinstance(meshes,(list,dict)):
            if isinstance(meshes,(list,)):
                self.mesh_list=meshes
            if isinstance(meshes,(dict,)):
                self.mesh_list=meshes.keys()
            if(reference_index>len(meshes)-1):
                msg = 'There are fewer meshes than the given reference_index'
                raise OSError(msg)
            else:
                 self.reference_index=reference_index


