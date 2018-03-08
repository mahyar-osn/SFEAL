import numpy
import morphic

class SSM(object):

    def __init__(self):
        self.X = []
        self.input_mesh = None
        self.groups = None
        self.mesh = None
        self.pcamesh = None
        self.pmesh = None
        self.score_0 = None
        self.z_score = []
        self.ratio = {}

    def add_mesh(self, mesh, index = 0):
        mesh = morphic.Mesh(str(mesh))
        if self.input_mesh == None:
            self.input_mesh = mesh
        if isinstance(mesh, str):
            mesh = morphic.Mesh(mesh)
        x = []
        if self.groups is None:
            for node in mesh.nodes:
                if not isinstance(node, morphic.mesher.DepNode):
                    x.extend(node.values.flatten().tolist())

        else:
            for node in mesh.nodes:
                if node.in_group(self.groups):
                    x.extend(node.values.flatten().tolist())

        self.X.append(x)

    def pca_setup(self, num_modes = 2):
        from sklearn import decomposition
        self.X = numpy.array(self.X)
        self.num_modes = num_modes
        self.pca = decomposition.PCA(n_components=num_modes)
        self.pca.fit(self.X)
        self.mean = self.pca.mean_
        self.components = self.pca.components_.T
        self.variance = self.pca.explained_variance_
        self.generate_mesh()
        return (self.mesh, self.X)

    def generate_mesh(self):
        self.mesh = morphic.Mesh()
        weights = numpy.zeros(self.num_modes + 1)
        weights[0] = 1.0
        self.mesh.add_stdnode('weights', weights)
        variance = numpy.zeros(self.num_modes + 1)
        variance[0] = 1.0
        variance[1:] = numpy.sqrt(self.variance)
        self.mesh.add_stdnode('variance', variance)
        idx = 0
        if self.groups is None:
            for node in self.input_mesh.nodes:
                nsize = node.values.size
                x = self.get_pca_node_values(node, idx)
                self.mesh.add_pcanode(node.id, x, 'weights', 'variance', group='pca')
                idx += nsize

        else:
            for node in self.input_mesh.nodes:
                nsize = node.values.size
                if node.in_group(self.groups):
                    x = self.get_pca_node_values(node, idx)
                    self.mesh.add_pcanode(node.id, x, 'weights', 'variance', group='pca')
                    idx += nsize
                else:
                    if isinstance(node, morphic.mesher.StdNode):
                        self.mesh.add_stdnode(node.id, node.values)
                    elif isinstance(node, morphic.mesher.DepNode):
                        self.mesh.add_depnode(node.id, node.element, node.node, shape=node.shape, scale=node.scale)
                    if isinstance(node, morphic.mesher.PCANode):
                        raise Exception('Not implemented')

        for element in self.input_mesh.elements:
            self.mesh.add_element(element.id, element.basis, element.node_ids)

        self.mesh.generate()

    def get_pca_node_values(self, node, idx):
        nsize = node.values.size
        if len(node.shape) == 1:
            pca_node_shape = (node.shape[0], 1, self.num_modes)
            x = numpy.zeros((node.shape[0], 1, self.num_modes + 1))
            x[:, 0, 0] = self.mean[idx:idx + nsize].reshape(node.shape)
            x[:, :, 1:] = self.components[idx:idx + nsize, :].reshape(pca_node_shape)
            return x
        if len(node.shape) == 2:
            pca_node_shape = (node.shape[0], node.shape[1], self.num_modes)
            x = numpy.zeros((node.shape[0], node.shape[1], self.num_modes + 1))
            x[:, :, 0] = self.mean[idx:idx + nsize].reshape(node.shape)
            x[:, :, 1:] = self.components[idx:idx + nsize, :].reshape(pca_node_shape)
            return x
        print 'Cannot reshape this node when generating pca mesh'

    def calculate_score(self, mesh_file_names, mesh_file, save = False):
        subject_name = mesh_file
        print '\n\t=========================================\n'
        print '\t   Please wait... \n'
        totalSubjects = 0
        subjectStore = []
        x = []
        for i in range(len(mesh_file_names)):
            single_mesh = mesh_file_names[i]
            mesh = morphic.Mesh(str(single_mesh))
            totalSubjects += 1
            subjectStore.append(mesh_file_names[i])
            nodes = mesh.get_nodes()
            size = len(nodes)
            for node in mesh.nodes:
                x.extend(node.values)
                X1 = numpy.asarray(x)

        X = X1.reshape((totalSubjects, size * 12))
        print '\t   Total number of subjects in pca = %d' % totalSubjects
        num_modes = totalSubjects - 1
        from sklearn import decomposition
        pca = decomposition.PCA(n_components=num_modes)
        pca.fit(X)
        pca_mean = pca.mean_
        pca_mean = pca_mean.reshape((1, size * 12))
        pca_components = pca.components_.T
        pca_variance = pca.explained_variance_
        pca_explained_variance = pca.explained_variance_ratio_
        self.ratio = {'MODE_{}'.format(m + 1):'{:.2f}'.format(float(pca_explained_variance[m])) for m in range(len(pca_explained_variance))}
        datasets = dict()
        for i in range(len(subjectStore)):
            datasets.update({subjectStore[i]: i})

        count = len(pca_variance)
        modeCount = []
        for i in range(len(pca_variance)):
            modeCount.append(i + 1)

        print '\t   Total modes of variation = %d' % count
        print '\t   Projecting Subject: %s' % subject_name
        mode_scores = []
        mah_dist = []
        for j in range(len(datasets)):
            mean = numpy.average(X[j], axis=0)
            subject = X[j] - pca_mean
            score = numpy.dot(subject, pca_components)
            mode_scores.append(score[0][0:count])

        self.SD = numpy.std(mode_scores, axis=0)
        self.mean = numpy.average(mode_scores, axis=0)
        number = datasets[subject_name]
        subject_0 = X[number] - pca_mean
        self.score_0 = numpy.dot(subject_0, pca_components)
        self.score_0 = self.score_0[0][0:count]
        self.score_z = self.convert_scores(self.score_0, self.SD, self.mean)
        return (self.score_z, self.ratio)

    def convert_scores(self, scores, SD, mean):
        for i in range(len(scores)):
            self.z_score.append((scores[i] - mean[i]) / SD[i])

        return self.z_score