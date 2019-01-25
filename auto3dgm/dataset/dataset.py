class Dataset:
    #params: self,
    #meshes: either a mesh object or a list or a dictionary of mesh objects
    #analyses: an analysis object (list, dictionary of such) associated with the mesh objects
    def __init__(self, meshes, analyses={}):
        if isinstance(meshes,dict):
            mesh_set=meshes
            if isinstance(analyses,dict):
                analysis_set=analyses
        elif isinstance(meshes, (list,)):
            ID=0
            self.mesh_set={}
            self.analysis_set={}
            for i in range(len(meshes)):
                ID=ID+1
                self.mesh_set[ID]=meshes[i]
                if (isinstance(analyses, (list,)) and len(meshes)==len(analyses)): 
                    self.analysis_set[ID]=analyses[i]
        elif hasattr(meshes,'vertices'):
            self.mesh_set={1:meshes}
            self.analysis_set={1:analyses} 
        else:
            msg = 'Input not a list nor a valid mesh object'
            raise OSError(msg)
