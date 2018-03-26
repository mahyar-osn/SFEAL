
from sfeal import core as sf
from sfeal import visualise as fig
reload(fig)
import os


"""
Directory setup
"""
root = "/hpc/mosa004/Lung/Data"  # specify where the root directory for lung meshes are
study = 'Human_Lung_Atlas'  # specify the study
volume = 'TLC'  # specify the imaged volume
path = os.path.join(volume, 'Lung/SurfaceFEMesh')
subjects = ['P2BRP032-H684', 'P2BRP001-H653', 'P2BRP143-H683', 'P2BRP241-H11296']  # list the subjects

"""
Mesh Generation & Alignment
"""
myMesh = sf.MESH()

for sub in range(len(subjects)):
    full_path = os.path.join(root, study, subjects[sub], path)  # setting up the path to the subjects fitted mesh
    output_path = myMesh.convert_cm_mesh(full_path, lung='LR')  # converting cm mesh into matrices
    myMesh.generate_mesh(output_path, lung='L', save=True)  # generate mesh from the matrices

meshes = []
for sub in range(len(subjects)):
    full_path = os.path.join(root, study, subjects[sub], path, 'morphic_original')
    for mesh_file in os.listdir(full_path):
        if mesh_file.endswith("Lung_fitted.mesh"):
            meshes.append(os.path.join(full_path, mesh_file))

reference_mesh = meshes[0]  # randon mesh as reference for alignment
d, Z, tform, mesh = [myMesh.align_mesh(reference_mesh, meshes[i], scaling=True, reflection='best')
                     for i in range(len(meshes))]  # mesh alignment (a General Procrsutes Alignment)

"""
Statistical Shape Modeling (PCA training | Score calculation | New mesh projection)
"""
mySF = sf.SSM()  # Create a SFEAL SSM object

meshes = []
for sub in range(len(subjects)):
    full_path = os.path.join(root, study, subjects[sub], path, 'morphic_aligned')
    for mesh_file in os.listdir(full_path):
        if mesh_file.endswith("Lung_fitted.mesh"):
            meshes.append(os.path.join(full_path, mesh_file))

for m in range(len(meshes)):
    mySF.add_mesh(meshes[m])  # Add the meshes

pmesh, _ = mySF.pca_train()  # Train the PCA and store in pmesh

score, ratio = mySF.calculate_score(meshes, meshes[0], save=False)  # Calculate scores

weights = [1.5]
mySF.export_to_cm(pmesh, weights, name='default', lung='L', show_mesh=True)  # Export a mesh to a CMISS mesh

new_mesh = "hpc/mosa004/Lung/Human_Aging/test/AGING041/Insp/Lung/SurfaceFEMesh/morphic_aligned/Left_fitted.mesh"
score, ratio = mySF.project_new_mesh(meshes, new_mesh)  # Project a new mesh on to the pca to calculate its score

"""
Data visualization
"""
myFig = fig.FIGURE()

myFig.animation(pmesh, 1, -1, 1, t=10, fissure=False)  # plotting lung pca animation

myFig.spectrum(pmesh, 1, -1, 1, fissure=False)  # plotting shape variation as color spectrum

myFig.show_lung(meshes[0], fissure=False)  # plotting a specific lung
