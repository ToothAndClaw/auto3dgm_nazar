
print("Begin Test")

factory = meshfactory.MeshFactory()
x = factory.mesh_from_file("test.ply")
print("This is mesh"), x
print("This is mesh.faces"), x.faces
print("This is mesh.vetices"), x.vertices