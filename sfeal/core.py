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
        self.score_1 = []
        self.z_score = {}
        self.ratio = {}
        self.fname = None
        self.mean = None
        self.score_z = None
        self.SD = None
        self.nodes = None
        self.new_data = []
        self.weights = []
        self.lung = None

    def add_mesh(self, mesh, index=0):
        mesh = morphic.Mesh(str(mesh))
        if self.input_mesh is None:
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

    def pca_train(self, num_modes=2):
        from sklearn import decomposition
        self.X = numpy.array(self.X)
        self.num_modes = num_modes
        self.pca = decomposition.PCA(n_components=num_modes)
        self.pca.fit(self.X)
        self.mean = self.pca.mean_
        self.components = self.pca.components_.T
        self.variance = self.pca.explained_variance_
        self.generate_mesh()
        return self.mesh, self.X

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
        print ('Cannot reshape this node when generating pca mesh')

    def calculate_score(self, mesh_file_names, mesh_file, save=False):
        if not self.new_data:
            pass
        else:
            self.new_data = []
        subject_name = mesh_file
        print ('\n\t=========================================\n')
        print ('\t   Please wait... \n')
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
        print ('\t   Total number of subjects in pca = %d') % totalSubjects
        num_modes = totalSubjects - 1
        from sklearn import decomposition
        pca = decomposition.PCA(n_components=num_modes)
        pca.fit(X)
        pca_mean = pca.mean_
        pca_mean = pca_mean.reshape((1, size * 12))
        pca_components = pca.components_.T
        pca_variance = pca.explained_variance_
        pca_explained_variance = pca.explained_variance_ratio_

        print 'PCA COMPONENETS'
        print pca_components

        self.ratio = {}
        self.ratio = {'MODE_{} RATIO'.format(m + 1): '{:.2f}'.format(float(pca_explained_variance[m])) for m in
                      range(len(pca_explained_variance))}
        dataset = dict()
        for i in range(len(subjectStore)):
            dataset.update({subjectStore[i]: i})

        count = len(pca_variance)
        mode_count = []
        for i in range(len(pca_variance)):
            mode_count.append(i + 1)

        print ('\t   Total modes of variation = %d') % count
        print ('\t   Projecting Subject: %s') % subject_name
        mode_scores = []
        for j in range(len(dataset)):
            subject = X[j] - pca_mean
            score = numpy.dot(subject, pca_components)
            mode_scores.append(score[0][0:count])

        self.SD = numpy.std(mode_scores, axis=0)
        self.mean = numpy.average(mode_scores, axis=0)
        number = dataset[subject_name]
        subject_0 = X[number] - pca_mean
        self.score_0 = numpy.dot(subject_0, pca_components)
        self.score_0 = self.score_0[0][0:count]
        self.score_1 = []
        self.score_1 = self.convert_scores(self.score_0, self.SD, self.mean)
        self.score_z = {}
        self.score_z = {'MODE_{} SCORE'.format(m + 1): '{:.2f}'.format(float(self.score_1[m])) for m in
                        range(len(self.score_1))}

        if save:
            # save scores to a csv file

            print ('\n\t   ______________\n')
            print ('\t   SORRY!')
            print ('\t   SAVE OPTION NOT IMPLEMENTED YET')
            print ('\n\t   ______________\n')

            """
            TODO

            import os
            import csv
            if self.fname == None:
                print '\n\t   ______________\n'
                self.fname = raw_input('\t   PLEASE NAME YOUR FILE: ')
            else:
                output_dir = 'output/'
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                weights_file = output_dir+'%s.csv' % self.fname
                if os.path.isfile(weights_file):
                    print '\n\t   ______________\n'
                    print '\t   File already exists'
                    print '\t   Appending...\n'
                else:
                    weights_file = output_dir+'%s.csv' % self.fname
                    text_file = open(weights_file, "a")
                    text_file.write("Mesh")
                    text_file.close()
                    for i in range(len(mode_count)):
                        f = open(weights_file)
                        data = [item for item in csv.reader(f)]
                        f.close()
                        new_column = mode_count
                        self.new_data = []
                        for j, item in enumerate(data):
                            try:
                                item.append(new_column[i])
                            except IndexError, e:
                                item.append("placeholder")
                            self.new_data.append(item)
                        f = open(weights_file, 'w')
                        csv.writer(f).writerows(self.new_data)
                        f.close()
                if not self.new_data:
                    pass
                else:
                    self.new_data = []
                # saving = True
                # if saving:
                text_file = open(weights_file, "a")
                text_file.write("%s" % subject_name)
                text_file.close()

                # if os.path.isfile(weights_file):
                #     self.new_data = []
                #     new = []
                #     for i in range(len(mode_count)):
                #         f = open(weights_file)
                #         data = [item for item in csv.reader(f)]
                #         f.close()
                #
                #         data = [data[-1]]
                #         for j, item in enumerate(data):
                #             try:
                #                 item.append("%.2f" % self.score_z[i])
                #             except IndexError, e:
                #                 item.append("placeholder")
                #             self.new_data.append(item)
                #         f = open(weights_file, 'a')
                #         csv.writer(f).writerows(self.new_data)
                #         f.close()
                print ("\t   Scores saved to %s \n") % weights_file
                print ('\t   ______________')
                
            """
        print ('\n\t=========================================\n')
        return self.score_z, self.ratio

    def convert_scores(self, scores, SD, mean):
        self.score_1 = []
        for i in range(len(scores)):
            self.score_1.append((scores[i] - mean[i]) / SD[i])

        return self.score_1

    def export_to_cm(self, pmesh, weights, name='default', lung='l', show_mesh=False):
        if not self.weights:
            pass
        else:
            self.weights = []

        if lung == 'R':
            self.lung = 'Right'
        elif lung == 'L':
            self.lung = 'Left'
        else:
            raise Exception("'lung' argument can ONLY be L OR R!")

        self.weights = weights
        self.pmesh = pmesh
        self.pmesh.nodes['weights'].values[1:] = 0  # reset weights to zero
        self.pmesh.nodes['weights'].values[0] = 1  # adding average
        for numMode in range(len(self.weights)):
            self.pmesh.nodes['weights'].values[numMode + 1] = self.weights[numMode]
            self.pmesh.update_pca_nodes()

        # saving
        import os
        import pandas as pd
        import subprocess
        from useful_files import nodes

        if self.nodes is not None:
            self.nodes = None

        self.nodes = nodes.Nodes()

        output_dir = 'output/export_to_cm/%s' % name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        input_folder = '../useful_files'

        path_to_export_mesh = output_dir
        temp_file = '%s_reconstructed_temp.csv' % self.lung
        output_file = '%s_reconstructed' % self.lung
        save_temp_file = output_dir + '/%s' % temp_file
        save_output_file = output_dir + '/%s' % output_file
        ipnode_file = '%s_reconstructed' % self.lung
        path_to_ipnode_file = '%s/%s' % (output_dir, ipnode_file)

        path_to_com_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'ipnode2exnode.com')
        path_to_cmgui_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'cmgui.com')
        ip2ex_perl = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'ipnode2exnode.pl')
        ip2ex_cm = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'ipnode2exnode')
        cmgui_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'cmgui')
        param_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', '3d_fitting')
        versions_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'versions')
        base_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'BiCubic_Surface_Unit')

        if self.lung == 'Right':
            node_file = 'nodes_%s.csv' % self.lung
            input_file = os.path.join(os.path.dirname(__file__), input_folder, node_file)
            elem_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'templateRight')

            nodes = self.nodes.set_nodes(lung='right')
            for node_number in nodes:
                node = self.pmesh.nodes[node_number]
                nodeValues = node.values
                with open(save_temp_file, 'a') as f:
                    numpy.savetxt(f, nodeValues)

            a = pd.read_csv(input_file)
            b = pd.read_csv(save_temp_file, delimiter=' ')
            result = pd.concat([a, b], axis=1)
            result.to_csv('%s.csv' % save_output_file, sep=' ', index=False)
            os.remove(save_temp_file)
            py2ip_right_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'py2ip_right.pl')

            subprocess.call(["perl", py2ip_right_file, "%s.csv" % save_output_file, "%s.ipnode" % path_to_ipnode_file])

            with open(path_to_com_file, 'wb') as comfile:
                comfile.write(" set echo on;\n")
                comfile.write(" fem def param;r;{0}".format("%s;\n" % param_file))
                comfile.write(" fem def coor;r;{0}".format("%s;\n" % versions_file))
                comfile.write(" fem def base;r;{0}".format("%s;\n" % base_file))
                comfile.write(" fem def node;r;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem def elem;r;{0}".format("%s;\n" % elem_file))
                comfile.write(" fem export node;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem export elem;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem def node;w;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem quit;\n")

            if show_mesh:
                with open(path_to_cmgui_file, 'wb') as comfile:
                    comfile.write(" gfx read node {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx read elem {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx cre egroup fissure;\n")
                    comfile.write(" gfx mod egroup fissure add 51..62;\n")
                    comfile.write(
                        " gfx mod g_e {0} general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e {0} lines coordinate coordinates select_on material green selected_material default_selected;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e fissure general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n")
                    comfile.write(" gfx mod g_e fissure surfaces material tissue;\n")
                    comfile.write(" gfx edit scene;\n")
                    comfile.write(" gfx cre win;\n")
                show_cmgui = 'show'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui],
                                shell=False)

            else:
                show_cmgui = 'no'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui], shell=True)

            print ("\n\t=========================================\n")
            print ("\t   ALL MESH FILES EXPORTED TO:")
            print ("\n\t\t   %s ") % path_to_export_mesh
            print ("\n\t=========================================\n")

            os.remove(save_output_file + '.csv')

        elif self.lung == 'Left':
            node_file = 'nodes_%s.csv' % self.lung
            input_file = os.path.join(os.path.dirname(__file__), input_folder, node_file)
            elem_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'templateLeft')

            nodes = self.nodes.set_nodes(lung='left')
            for node_number in nodes:
                node = self.pmesh.nodes[node_number]
                nodeValues = node.values
                with open(save_temp_file, 'a') as f:
                    numpy.savetxt(f, nodeValues)

            a = pd.read_csv(input_file)
            b = pd.read_csv(save_temp_file, delimiter=' ')
            result = pd.concat([a, b], axis=1)
            result.to_csv('%s.csv' % save_output_file, sep=' ', index=False)
            os.remove(save_temp_file)

            py2ip_left_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'py2ip_left.pl')
            subprocess.call(["perl", py2ip_left_file, "%s.csv" % save_output_file, "%s.ipnode" % path_to_ipnode_file])

            with open(path_to_com_file, 'wb') as comfile:
                comfile.write(" set echo on;\n")
                comfile.write(" fem def param;r;{0}".format("%s;\n" % param_file))
                comfile.write(" fem def coor;r;{0}".format("%s;\n" % versions_file))
                comfile.write(" fem def base;r;{0}".format("%s;\n" % base_file))
                comfile.write(" fem def node;r;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem def elem;r;{0}".format("%s;\n" % elem_file))
                comfile.write(" fem export node;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem export elem;{0} as {1};\n".format("%s" % path_to_ipnode_file, self.lung))
                comfile.write(" fem def node;w;{0}".format("%s;\n" % path_to_ipnode_file))
                comfile.write(" fem quit;\n")

            if show_mesh:
                with open(path_to_cmgui_file, 'wb') as comfile:
                    comfile.write(" gfx read node {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx read elem {0}".format("'%s';\n" % path_to_ipnode_file))
                    comfile.write(" gfx cre egroup fissure;\n")
                    comfile.write(" gfx mod egroup fissure add 111..118;\n")
                    comfile.write(
                        " gfx mod g_e {0} general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e {0} lines coordinate coordinates select_on material green selected_material default_selected;\n".format(
                            "'%s'" % self.lung))
                    comfile.write(
                        " gfx mod g_e fissure general clear circle_discretization 6 default_coordinate coordinates; element_discretization '12*12*12' native_discretization none;\n")
                    comfile.write(" gfx mod g_e fissure surfaces material tissue;\n")
                    comfile.write(" gfx edit scene;\n")
                    comfile.write(" gfx cre win;\n")
                show_cmgui = 'show'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui],
                                shell=False)

            else:
                show_cmgui = 'no'
                subprocess.call(["perl", ip2ex_perl, "%s" % ip2ex_cm, "%s" % cmgui_file, "%s" % show_cmgui], shell=True)

            print ("\n\t=========================================\n")
            print ("\t   ALL MESH FILES EXPORTED TO:")
            print ("\n\t\t   %s ") % path_to_export_mesh
            print ("\n\t=========================================\n")

            os.remove(save_output_file + '.csv')

        return None

    def project_new_mesh(self, mesh_file_names, mesh_file):
        if not self.new_data:
            pass
        else:
            self.new_data = []
        subject_name = mesh_file
        print ('\n\t=========================================\n')
        print ('\t   Please wait... \n')
        totalSubjects = 0
        subjectStore = []
        x = []
        y = []
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
        print ('\t   Total number of subjects in pca = %d') % totalSubjects
        num_modes = totalSubjects - 1
        from sklearn import decomposition
        pca = decomposition.PCA(n_components=num_modes)
        pca.fit(X)
        pca_mean = pca.mean_
        pca_mean = pca_mean.reshape((1, size * 12))
        pca_components = pca.components_.T
        pca_variance = pca.explained_variance_
        pca_explained_variance = pca.explained_variance_ratio_

        self.ratio = {}
        self.ratio = {'MODE_{} RATIO'.format(m + 1): '{:.2f}'.format(float(pca_explained_variance[m])) for m in
                      range(len(pca_explained_variance))}
        dataset = dict()
        for i in range(len(subjectStore)):
            dataset.update({subjectStore[i]: i})

        count = len(pca_variance)
        mode_count = []
        for i in range(len(pca_variance)):
            mode_count.append(i + 1)

        print ('\t   Total modes of variation = %d') % count
        print ('\t   Projecting Subject: %s') % subject_name

        mode_scores = []
        for j in range(len(dataset)):
            subject = X[j] - pca_mean
            score = numpy.dot(subject, pca_components)
            mode_scores.append(score[0][0:count])

        if self.SD is not None:
            self.SD = None
        if self.mean is not None:
            self.mean = None

        self.SD = numpy.std(mode_scores, axis=0)
        self.mean = numpy.average(mode_scores, axis=0)

        project_mesh_path = mesh_file
        project_mesh = morphic.Mesh(project_mesh_path)
        for node in project_mesh.nodes:
            y.extend(node.values)
            Y = numpy.asarray(y)

        Y = Y.reshape((size * 12))
        subject_0 = Y - pca_mean

        if self.score_0 is not None:
            self.score_0 = None
        if self.score_z is not None:
            self.score_z = None

        self.score_0 = numpy.dot(subject_0, pca_components)
        self.score_0 = self.score_0[0][0:count]
        self.score_1 = []
        self.score_1 = self.convert_scores(self.score_0, self.SD, self.mean)
        self.score_z = {}
        self.score_z = {'MODE_{} SCORE'.format(m + 1): '{:.2f}'.format(float(self.score_1[m])) for m in
                        range(len(self.score_1))}

        print ('\n\t=========================================\n')
        return self.score_z, self.ratio


class MESH(object):

    def __init__(self):
        self.lung = None
        self.count = 0
        self.elements = None
        self.mesh = None
        self.output = None
        self.file_path = None

    def generate_mesh(self, file_path, lung='L', save=True):
        """
        generate_mesh creates morphic meshes from finite element meshes
        built in cmiss. The input mesh which has already been converted
        into a matrix format is transformed into a morphic mesh using the
        node and element class within the useful_files module. These
        classes can be modified accordingly for other mesh topologies.

        Inputs:
        ------------
        file_path
            path where the .ipnode file is stored.

        file_name
            name of the .ipnode file (do not include the .ipnode itself).

        lung
            Left = L | Right = R | Both = LR

        save
            A boolean variable to save the mesh.

        Outputs:
        ------------
        None

        :param file_path
        :param file_name
        :param lung
        :param save
        :return: None
        """
        import csv
        import os, sys
        from useful_files import elements

        if lung == 'L' or lung == 'l':
            self.lung = 'Left'
        elif lung == 'R' or lung == 'r':
            self.lung = 'Right'
        elif lung == 'LR' or lung == 'lr' or lung == 'RL' or lung == 'rl':
            self.lung = 'Lung'

        if self.mesh is not None:
            self.mesh = None

        self.mesh = morphic.Mesh()
        data = {}

        if self.elements is not None:
            self.elements = None

        self.elements = elements.Elements()

        print ('\n\t=========================================\n')
        print ('\t   GENERATING MESH... \n')
        print ('\t   PLEASE WAIT... \n')

        if self.file_path is not None:
            self.file_path = None

        self.file_path = file_path

        for filenum in os.listdir(self.file_path):
            filenum_path = os.path.join(self.file_path, filenum)
            if filenum_path == os.path.join(self.file_path, self.lung + '_fitted.ip2py'):
                if os.path.isfile(filenum_path):
                    self.count += 1
                    with open(filenum_path, 'r') as csvfile:
                        data[filenum] = csv.reader(csvfile, delimiter=' ', quotechar='|')
                        for rowx in data[filenum]:
                            rowy = data[filenum].next()
                            rowz = data[filenum].next()
                            node = [[float(rowx[1]), float(rowx[2]), float(rowx[3]), float(rowx[4])],
                                    [float(rowy[1]), float(rowy[2]), float(rowy[3]), float(rowy[4])],
                                    [float(rowz[1]), float(rowz[2]), float(rowz[3]), float(rowz[4])]]
                            nd = self.mesh.add_stdnode(str(rowx[0]), node)

                        if self.lung == 'Left':
                            elements = self.elements.set_elements(lung='left')
                            for ii, elem in enumerate(elements):
                                self.mesh.add_element(ii + 1, ['H3', 'H3'], elem)
                        elif self.lung == 'Right':
                            elements = self.elements.set_elements(lung='right')
                            for ii, elem in enumerate(elements):
                                self.mesh.add_element(ii + 1, ['H3', 'H3'], elem)
                        else:
                            elements = self.elements.set_elements(lung='lr')
                            for ii, elem in enumerate(elements):
                                self.mesh.add_element(ii + 1, ['H3', 'H3'], elem)

                        self.mesh.generate()

                        if save:
                            meshOutput = os.path.normpath(filenum_path + os.sep + os.pardir)
                            self.mesh.save(meshOutput + '/' + self.lung + '_fitted.mesh')

                            print ('\t   MESH SAVED IN \n')
                            print ('\t   %s DIRECTORY \n') % meshOutput
        print ('\n\t=========================================\n')

    def align_mesh(self, reference_mesh, mesh, scaling=True, reflection='best'):
        """
        align_mesh is a method that perfomes a Procrustes analysis which
        determines a linear transformation (translation, reflection, orthogonal rotation
        and scaling) of the nodes in mesh to best conform them to the nodes in reference_mesh,
        using the sum of squared errors as the goodness of fit criterion.

        Inputs:
        ------------
        reference_mesh, mesh
            meshes (as morphic meshes) of target and input coordinates. they must have equal
            numbers of  nodes (rows), but mesh may have fewer dimensions
            (columns) than reference_mesh.

        scaling
            if False, the scaling component of the transformation is forced
            to 1

        reflection
            if 'best' (default), the transformation solution may or may not
            include a reflection component, depending on which fits the data
            best. setting reflection to True or False forces a solution with
            reflection or no reflection respectively.

        Outputs
        ------------
        d
            the residual sum of squared errors, normalized according to a
            measure of the scale of reference_mesh, ((reference_mesh - reference_mesh.mean(0))**2).sum()

        Z
            the matrix of transformed Y-values

        tform
            a dict specifying the rotation, translation and scaling that
            maps X --> Y

        self.mesh
            Aligned mesh

        :param reference_mesh
        :param mesh
        :param scaling
        :param reflection
        :return: d, Z, tform, self.mesh
        """

        import os

        print ('\n\t=========================================\n')
        print ('\t   ALIGNING MESH... \n')
        print ('\t   PLEASE WAIT... \n')

        if self.mesh is not None:
            self.mesh = None

        r = morphic.Mesh(reference_mesh)
        self.mesh = morphic.Mesh(mesh)

        X = r.get_nodes()
        Y = self.mesh.get_nodes()

        n, m = X.shape
        ny, my = Y.shape

        muX = X.mean(0)
        muY = Y.mean(0)

        X0 = X - muX
        Y0 = Y - muY

        ssX = (X0 ** 2.).sum()
        ssY = (Y0 ** 2.).sum()

        # centred Frobenius norm
        normX = numpy.sqrt(ssX)
        normY = numpy.sqrt(ssY)

        # scale to equal (unit) norm
        X0 /= normX
        Y0 /= normY

        if my < m:
            Y0 = numpy.concatenate((Y0, numpy.zeros(n, m - my)), 0)

        # optimum rotation matrix of Y
        A = numpy.dot(X0.T, Y0)
        U, s, Vt = numpy.linalg.svd(A, full_matrices=False)
        V = Vt.T
        T = numpy.dot(V, U.T)

        if reflection is not 'best':

            # does the current solution use a reflection?
            have_reflection = numpy.linalg.det(T) < 0

            # if that's not what was specified, force another reflection
            if reflection != have_reflection:
                V[:, -1] *= -1
                s[-1] *= -1
                T = numpy.dot(V, U.T)

        traceTA = s.sum()

        if scaling:

            # optimum scaling of Y
            b = traceTA * normX / normY

            # standarised distance between X and b*Y*T + c
            d = 1 - traceTA ** 2

            # transformed coords
            Z = normX * traceTA * numpy.dot(Y0, T) + muX

        else:
            b = 1
            d = 1 + ssY / ssX - 2 * traceTA * normY / normX
            Z = normY * numpy.dot(Y0, T) + muX

        # translation matrix
        if my < m:
            T = T[:my, :]
        c = muX - b * numpy.dot(muY, T)

        # transformation values
        tform = {'rotation': T, 'scale': b, 'translation': c}

        for num, object in enumerate(self.mesh.nodes):
            node = self.mesh.nodes[object.id].values[:, 0]
            Zlist = Z.tolist()
            self.mesh.nodes[object.id].values[:, 0] = Zlist[num]

        if self.output is not None:
            self.output = None

        self.output = 'morphic_aligned'

        meshOutput = os.path.normpath(mesh + os.sep + os.pardir)
        meshOutput = os.path.normpath(meshOutput + os.sep + os.pardir)

        meshOutput = os.path.join(meshOutput, self.output)

        if not os.path.exists(meshOutput):
            os.makedirs(meshOutput)

        meshName = mesh.split('/')

        self.mesh.save(os.path.join(meshOutput, str(meshName[-1])))

        print ('\t   ALIGNED MESH SAVED IN \n')
        print ('\t   %s DIRECTORY \n') % meshOutput

        print ('\n\t=========================================\n')

        return d, Z, tform, self.mesh

    def convert_cm_mesh(self, file_path, lung='L'):
        """

        :param file_path:
        :param file_name:
        :param lung:
        :return:
        """
        import os

        if lung == 'L' or lung == 'l':
            output_path, _ = self.process_cm_mesh('Left', file_path)
        elif lung == 'R' or lung == 'r':
            output_path, _ = self.process_cm_mesh('Right', file_path)
        else:
            output_path, output_lung_left = self.process_cm_mesh('Left', file_path)
            output_path, output_lung_right = self.process_cm_mesh('Right', file_path)
            final_file_path = os.path.normpath(output_lung_left + os.sep + os.pardir)
            final_file = os.path.join(final_file_path, 'Lung_fitted.ip2py')

            with open(output_lung_right) as f:
                with open(final_file, "w") as f1:
                    for line in f:
                        f1.write(line)
            with open(output_lung_left) as f:
                with open(final_file, "a") as f1:
                    for line in f:
                        f1.write(line)
        return output_path

    def process_cm_mesh(self, lung, file_path):
        """

        :param lung:
        :return:
        """
        import os
        import subprocess

        if self.lung is not None:
            self.lung = None
        self.lung = lung

        input_folder = '../useful_files'
        output_folder = 'morphic_original'
        output_path = os.path.join(file_path, output_folder)
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        perl_file = os.path.join(os.path.dirname(__file__), input_folder, 'perl_com', 'ip2py_%s.pl' % self.lung)
        try:
            _inp = file_path + '/' + self.lung + '_fitted.ipnode'
            # input_lung = os.path.join(file_path, self.lung + '_fitted.ipnode')
        except Exception:
            _inp = file_path + '/fitted' + self.lung + '.ipnode'

        input_lung = _inp
        output_lung = os.path.join(output_path, self.lung + '_fitted.ip2py')
        subprocess.call(["perl", perl_file, input_lung, output_lung])

        return output_path, output_lung
