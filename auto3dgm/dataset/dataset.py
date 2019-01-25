class Dataset:
    #params: self,
    #meshes: either a mesh object or a list of mesh objects
    #analyses: either an analysis object or a list of analysis objects
    def __init__(self, meshes, analyses={}):
        if isinstance(meshes, (list,)):
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
