"""

Classes to prepare plots for PEC. This code is arranged in classes to allow further improvements in a tidier way.

"""
import matplotlib

matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pec_utilities


class MorphologyVisualiser:

    def __init__(self, morph, normal, hkl, point, new_facet_edges, new_facet_ordered, crystal, area):

        ax = self.generate_BFDH_plot(morph, 'blue')
        self._add_plane(normal=normal, hkl=hkl, point=point, ax=ax)
        self._plot_intercept_points(intercept_points=new_facet_edges, ax=ax)
        facet = np.vstack(new_facet_ordered)
        Axes3D.plot(ax,
                    facet[:, 0],
                    facet[:, 1],
                    facet[:, 2])
        plt.title(f'Crystal: {crystal.identifier} \n Facet {hkl} area : {round(area, 3)}')

    @staticmethod
    def _plot_intercept_points(intercept_points, ax):
        """
        Plots the points where the plane intercepts the BFDH

        :param intercept_points: list intercept points across BFDH
        :param ax: obj axes to plot on
        :return: plots on axes
        """
        Axes3D.scatter3D(ax, xs=[point[0] for point in intercept_points],
                         ys=[point[1] for point in intercept_points],
                         zs=[point[2] for point in intercept_points], color='red')

    @staticmethod
    def _wireframe(face_for_wireframe, ax):
        """
        Plot wireframe of morphology in matplotlib

        :param face_for_wireframe: obj containing the edges for a given facet
        :param ax: obj Axes for matplotlib
        :return: adds line to given ax
        """
        for n, edge in enumerate(face_for_wireframe.edges):
            Axes3D.plot(ax,
                        [coord[0] for coord in edge],
                        [coord[1] for coord in edge],
                        [coord[2] for coord in edge],
                        c='black',
                        linewidth=1.5)

    @staticmethod
    def _skin(face_for_skin, colour, ax):
        """
        Mesh the facets for matplotlib

        :param face_for_skin: obj containing the edges for a given facet
        :param colour: str facet colour of BDFH
        :param ax: obj Axes for matplotlib
        :return: list vertices [edge1, edge2, centroid]
        """
        vert_list = []
        for edge in face_for_skin.edges:
            vertices = [(edge[0], edge[1], face_for_skin.centre_of_geometry)]
            Axes3D.add_collection3d(ax, Poly3DCollection(vertices, color=colour, linewidth=0, alpha=0.2))
            vert_list.append(vertices)
        return vert_list

    @staticmethod
    def _add_labels(face_for_labels, ax):
        """
        Add labels to facets

        :param face_for_labels: obj containing the edges for a given facet
        :param ax: obj Axes for matplotlib
        :return: adds labels to given facets on ax
        """
        ax.text(face_for_labels.centre_of_geometry[0],
                face_for_labels.centre_of_geometry[1],
                face_for_labels.centre_of_geometry[2],
                '''{}'''.format(face_for_labels.miller_indices.hkl),
                color='black',
                size=10)

        ax.scatter(face_for_labels.centre_of_geometry[0],
                   face_for_labels.centre_of_geometry[1],
                   face_for_labels.centre_of_geometry[2],
                   s=5,
                   color='black')

    @staticmethod
    def generate_morphology_plot(morphology, colour="teal"):
        """
        Main driving function for plotting a morphology

        see https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/morphology.html#morphology

        :param morphology: obj containing BFDH facet information
        :param colour: str desired colour for facets
        :return: obj matplotlib object for plotting on axes
        """
        # Set up Axes
        fig = plt.figure(figsize=(10., 10.))  # 3D graph instance
        ax = fig.add_subplot(111, projection="3d", proj_type="ortho")  # 3D Axes
        ax.grid(False)
        plt.axis('on')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        #  Plot the morphology in matplotlib
        for facet in morphology.facets:
            MorphologyVisualiser._wireframe(facet, ax)  # Plot the edges
            MorphologyVisualiser._skin(facet, colour, ax)  # Mesh the faces
            MorphologyVisualiser._add_labels(facet, ax)  # Add labels to the faces

        # Orient and position the morphology in a way that is easy to view in a 2D image
        obb = morphology.oriented_bounding_box  # First we generate an oriented bounding box
        #  We use the box to scale the image to a good size
        extent = max([pec_utilities.distance([0, 0, 0], corner) for corner in obb.corners]) / 2
        ax.set_xlim3d(-1 * extent, extent)
        ax.set_ylim3d(-1 * extent, extent)
        ax.set_zlim3d(-1 * extent, extent)

        plt.show()
        return ax

    @staticmethod
    def _add_plane(point, hkl, normal, ax):
        """
        Plots a plane on the BDFH morphology that is equivalent to the one requested

        :type normal: numpy.array
        :param point: numpy.array centre of the plane
        :param hkl: numpy.array desired hkl of facet
        :param normal: numpy.array normal vector of the plane
        :param ax: obj graph axes to be plotted upon
        :return: plots to axes
        """
        d = -point.dot(normal)
        x = np.linspace(-25, 25, 10)
        y = np.linspace(-25, 25, 10)
        z = np.linspace(-25, 25, 10)

        X, Y, Z = None, None, None
        if hkl[2] != 0:
            X, Y = np.meshgrid(x, y)
            Z = (-normal[0] * X - normal[1] * Y - d) * 1. / normal[2]
        elif hkl[1] != 0:
            X, Z = np.meshgrid(x, z)
            Y = (-normal[0] * X - normal[2] * Z - d) * 1. / normal[1]
        elif hkl[0] != 0:
            Z, Y = np.meshgrid(z, y)
            X = (-normal[2] * Z - normal[1] * Y - d) * 1. / normal[0]
        else:
            pass
        Axes3D.plot_surface(ax, X, Y, Z, alpha=0.2)


class AspectRatioVisualiser:

    def __init__(self, morphologies, surface_energy_penalties):
        ax = self.generate_AR_plot(morphologies, surface_energy_penalties)

    @staticmethod
    def generate_AR_plot(morphologies, surface_energy_penalties, colour_map="hot", levels=10):
        fig, ax = plt.subplots()
        ax.grid(False)
        plt.axis('on')
        ax.set_xlabel('L / T')
        ax.set_ylabel('W / T')

        cmap = plt.get_cmap(colour_map)

        aspect_ratios = [pec_utilities.calculate_aspect_ratio(morph) for morph in morphologies]

        x = [item[0] for item in aspect_ratios]
        y = [item[1] for item in aspect_ratios]
        z = surface_energy_penalties

        limits = [1, 14.5]

        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[0], limits[1])
        # ax.set_title("Title")
        ax.set_xlabel("L / T")
        ax.set_ylabel("W / T")
        ax.set_aspect('equal')

        tcf = ax.tricontourf(x, y, z, levels=levels, cmap=cmap)
        cbar = fig.colorbar(tcf)
        cbar.ax.tick_params(labelsize=8)
        cbar.ax.set_ylabel("surface energy penalty (kJ mol$^{-1}$)", rotation=90)
        ax.tricontour(x, y, z, levels=levels, colors='k', linewidths=0.6)

        ax.plot(limits, limits, "k", lw=0.5)

        plt.show(block=True)

        return ax


class ParticleEnergyVisualiser:
    """ """

    def __init__(self, nanoparticles_1, nanoparticles_2):
        ax = self.generate_energy_plot(nanoparticles_1, nanoparticles_2)

    @staticmethod
    def generate_energy_plot(nanoparticles_1, nanoparticles_2, smallest_size=False):
        fig, ax = plt.subplots()
        ax.grid(False)
        plt.axis('on')

        ax.set_ylabel('Particle Energy (kJ/mol)')

        # we divide the particle size by 10 to express results in nm
        if smallest_size:
            ax.set_xlabel('Smallest size (nm)')
            size_1 = [particle.smallest_dimension / 10 for particle in nanoparticles_1]
            energy_1 = [particle.particle_energy for particle in nanoparticles_1]
            size_2 = [particle.smallest_dimension / 10 for particle in nanoparticles_2]
            energy_2 = [particle.particle_energy for particle in nanoparticles_2]
        else:
            size_1 = [particle.equivalent_diameter / 10 for particle in nanoparticles_1]
            energy_1 = [particle.particle_energy for particle in nanoparticles_1]
            size_2 = [particle.equivalent_diameter / 10 for particle in nanoparticles_2]
            energy_2 = [particle.particle_energy for particle in nanoparticles_2]
            ax.set_xlabel('PED (nm)')

        ax.axhline(y=nanoparticles_1[-1].bulk_energy, c="r", ls=":")
        ax.axhline(y=nanoparticles_2[-1].bulk_energy, c="k", ls=":")

        ax.plot(size_1, energy_1, lw=1.5, c="red", label="form 1")
        ax.plot(size_2, energy_2, lw=1.5, c="black", label="form 2")

        plt.legend(loc="best", frameon=False)

        plt.show(block=True)

        return ax


class CrossCheckVisualiser:

    def __init__(self, results):
        ax = self.generate_cross_check_plot(results)

    @staticmethod
    def generate_cross_check_plot(sizes, fractions):
        fig, ax = plt.subplots()

        # set the thickness of plot frame spines
        ax.spines['bottom'].set_linewidth(0.5)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['top'].set_linewidth(0.5)
        ax.spines['right'].set_linewidth(0.5)

        plt.ylim(0, 100)

        ax.set_box_aspect(1)  # if the axes are set to be of equal size this is superfluous

        ax.set_xlabel("PED (nm)", name='Arial')
        ax.set_ylabel("Fraction of stability switches (%)", name='Arial')
        ax.tick_params(axis='both', which='major', labelsize=8)

        ax.bar(sizes, fractions, color="#f3b61f", width=5, edgecolor="k", linewidth=0.5)

        plt.show(block=True)

        return ax
