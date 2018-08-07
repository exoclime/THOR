
import numpy as np
import h5py
import pathlib
import re

from math import sqrt
from utilities import HSV_to_RGB, H_to_RGB, spherical


def normalize(data):
    # rescale
    min_data = np.min(data)
    max_data = np.max(data)

    delta = max_data - min_data

    data[delta != 0.0] = (data[delta != 0] - min_data)/(max_data - min_data)
    return data


class simdataset:
    def __init__(self, folder):
        self.dataloader = dataloader(folder)

    def select_planet(self, planet):
        self.dataloader.select_planet(planet)

        self.currently_loaded_data = None

        self.num_samples = self.get_num_datasets()
        self.num_points, self.num_levels = self.get_dim()
        self.color_data_buffer = np.zeros(
            (self.num_samples, self.num_levels + 1, self.num_points), dtype=np.float32)

    def get_planets(self):
        return self.dataloader.get_planets()

    def get_num_datasets(self):
        return self.dataloader.get_num_datasets()

    def get_dim(self):
        pts, lvl = self.dataloader.get_dim()
        return pts, lvl + 1

    def get_grid(self):
        return self.dataloader.get_grid()

    def get_data_types(self):
        datatypes = {}
        datatypes["Pressure"] = {'type': 'scalar'}
        datatypes["Rho"] = {'type': 'scalar'}
        datatypes["Momentum"] = {'type': 'vector'}
        datatypes["Momentum_norm"] = {'type': 'scalar'}
        datatypes["Momentum_norm_vert"] = {'type': 'scalar'}
        datatypes["Momentum_norm_horiz"] = {'type': 'scalar'}
        datatypes["Momentum_horiz"] = {'type': 'vector'}
        datatypes["Momentum_vert"] = {'type': 'vector'}

        return datatypes

    def get_relative_radii(self):
        radius = self.dataloader.get_radius()
        altitudes = self.dataloader.get_altitudes()

        lonlat = self.dataloader.get_grid()

        alt_min = altitudes[0]
        alt_max = altitudes[-1]

        r_max = radius+alt_max
        rad = []
        for a in altitudes:
            rad.append((radius+a)/r_max)

        return rad

    def select_vector_data_type(self, dataname):
        self.field_dataname = dataname
        radius = self.dataloader.get_radius()
        altitudes = self.dataloader.get_altitudes()

        lonlat = self.dataloader.get_grid()

        radii = self.get_relative_radii()

        field = np.zeros((self.num_samples,
                          self.num_levels,
                          self.num_points,
                          2, 3), dtype=np.float32)

        if dataname == "Momentum":
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels+1), dtype=np.float32)

            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels -
                         1] = self.dataloader.get_data(i, "Mh")
                vector_v[i, :, :self.num_levels] = self.dataloader.get_data(
                    i, "Wh")

            # interpolate vertical value
            # blah

            # create field vector

            scale = 0.01
            for l in range(self.num_levels-1):
                for i in range(self.num_points):
                    normal = spherical(1.0,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    r = radii[l]
                    vertex = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    field[:, l, i, 0, :] = vertex
                    field[:, l, i, 1, :] = vertex + scale * \
                        (vector_h[:, i, l, :] +
                         normal)  # *vector_v[:, i, l])
                    # field[:, l, i, 1, :] = vertex + normal
                    # field[:, l, i, 0, :] = (0.0, 0.0, 0.0)
                    # field[:, l, i, 1, :] = normal

        elif dataname == "Momentum_horiz":
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels-1, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels), dtype=np.float32)

            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels-1] = self.dataloader.get_data(
                    i, "Mh")
                vector_v[i, :, :] = self.dataloader.get_data(
                    i, "Wh")

            # interpolate vertical value
            # blah

            # create field vector

            scale = 0.01
            for l in range(self.num_levels-1):
                for i in range(self.num_points):
                    normal = spherical(1.0,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    r = radii[l]
                    vertex = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    field[:, l, i, 0, :] = vertex
                    field[:, l, i, 1, :] = vertex + scale * \
                        vector_h[:, i, l, :]
                    # field[:, l, i, 1, :] = vertex + normal
                    # field[:, l, i, 0, :] = (0.0, 0.0, 0.0)
                    # field[:, l, i, 1, :] = normal

        elif dataname == "Momentum_vert":
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels-1, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels), dtype=np.float32)

            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels -
                         1] = self.dataloader.get_data(i, "Mh")
                vector_v[i, :, :] = self.dataloader.get_data(
                    i, "Wh")
                print("vectorvnan", np.isnan(vector_v))

            # interpolate vertical value
            # blah
            vector_v = normalize(vector_v)
            # create field vector
            scale = 0.01
            dr = 0.01
            for l in range(self.num_levels-1):
                for i in range(self.num_points):
                    normal = spherical(1.0,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    r = radii[l]
                    vertex = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    field[:, l, i, 0, :] = (1.0+dr)*vertex
                    field[:, l, i, 1, :] = (1.0+dr)*vertex + (scale *
                                                              vector_v[:, i, l]*normal.reshape((3, 1))).T
                    # print(vector_v)
                    # vertex + \
                    # scale * vector_v[:, i, l]*normal.reshape((1, 3))
                    # field[:, l, i, 1, :] = vertex + normal
                    # field[:, l, i, 0, :] = (0.0, 0.0, 0.0)
                    # field[:, l, i, 1, :] = normal
        else:
            print("Unrecognised data name for vector field: ", dataname)

        self.field_data = field

    def select_scalar_data_type(self, dataname):
        print("update scalar display: ", dataname)

        radius = self.dataloader.get_radius()
        altitudes = self.dataloader.get_altitudes()

        lonlat = self.dataloader.get_grid()

        radii = self.get_relative_radii()
        field = np.zeros((self.num_samples,
                          self.num_levels,
                          self.num_points,
                          2, 3), dtype=np.float32)
        self.dataname = dataname
        self.data_color = np.zeros(
            (self.num_samples, self.num_levels, self.num_points, 3), dtype=np.float32)
        self.data_scalar = np.zeros(
            (self.num_samples, self.num_levels, self.num_points), dtype=np.float32)
        # load data set
        if dataname == "Pressure":
            data = np.zeros((self.num_samples, self.num_points,
                             self.num_levels), dtype=np.float32)

            print(self.data_color.shape)
            for i in range(self.num_samples):
                data[i, :, :self.num_levels -
                     1] = self.dataloader.get_data(i, dataname)

            # rescale
            min_data = np.min(data)
            max_data = np.max(data)

            zeros = (max_data - min_data) == 0
            data = (data - min_data)/(max_data - min_data)
            data[zeros] = 0.0

            # self.data_color = np.array(
            #    np.swapaxes(H_to_RGB(data), 1, 2), copy=True, dtype=np.float32)
            self.data_scalar = np.array(np.swapaxes(data, 1, 2))
            for i in range(self.num_samples):
                for j in range(self.num_levels):
                    for k in range(self.num_points):
                        #self.data_color[i, j, k] = (data[i, k, j], 0.0, 0.0)
                        self.data_color[i, j, k] = HSV_to_RGB(
                            360.0*data[i, k, j], 1.0, 1.0)

        elif dataname == "Rho":
            data = np.zeros((self.num_samples, self.num_points,
                             self.num_levels), dtype=np.float32)

            for i in range(self.num_samples):
                data[i, :, :self.num_levels -
                     1] = self.dataloader.get_data(i, dataname)

            # rescale
            min_data = np.min(data)
            max_data = np.max(data)

            data = (data - min_data)/(max_data - min_data)
            self.data_scalar = np.array(np.swapaxes(data, 1, 2))
            self.data_color = np.array(np.swapaxes(
                H_to_RGB(data), 1, 2), copy=True, dtype=np.float32)

        elif dataname == "Momentum_norm":
            print("computing momentum norm")
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels+1), dtype=np.float32)
            field = np.zeros((self.num_samples,
                              self.num_levels,
                              self.num_points, 3), dtype=np.float32)

            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels -
                         1] = self.dataloader.get_data(i, "Mh")
                vector_v[i, :, :self.num_levels] = self.dataloader.get_data(
                    i, "Wh")

            # interpolate vertical value
            # blah

            # create field vector

            for l in range(self.num_levels-1):
                for i in range(self.num_points):
                    normal = spherical(1.0,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    r = radii[l]
                    vertex = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    field[:, l, i, :] = vector_h[:, i, l, :] + \
                        (vector_v[:, i, l]*normal.reshape((3, 1))).T
                    # field[:, l, i, 0, :] = (0.0, 0.0, 0.0)
                    # field[:, l, i, 1, :] = normal
            data = np.zeros((self.num_samples, self.num_levels, self.num_points,
                             3), dtype=np.float32)

            if np.any(np.isnan(data)):
                print("NaN")

            # compute norm
            data = np.sqrt(np.sum(np.power(field, 2.0), axis=3))
            self.data_scalar = data
            # data = data[:, :, :, 0]
            # rescale
            for i in range(self.num_samples):
                for j in range(self.num_levels):
                    d = data[i, j, :]
                    min_data = np.min(d)
                    max_data = np.max(d)

                    data[i, j, :] = (d - min_data)/(max_data - min_data)

                    if np.any(max_data - min_data == 0):
                        print("zeros")

            # min_data = np.min(data)
            # max_data = np.max(data)

            # print(min_data, max_data)
            # zeros = (max_data - min_data) == 0
            # non_zeros = (max_data - min_data) != 0
            # data[non_zeros] = (data[non_zeros] - min_data[non_zeros]) / \
            #     (max_data[non_zeros] - min_data[non_zeros])
            # data[zeros] = 0.0

            # self.data_color = data
            self.data_color[:, :, :, :] = np.array(
                H_to_RGB(360.0*data), copy=True, dtype=np.float32)
        elif dataname == "Momentum_norm_horiz":
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels+1), dtype=np.float32)
            field = np.zeros((self.num_samples,
                              self.num_levels,
                              self.num_points, 3), dtype=np.float32)

            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels -
                         1] = self.dataloader.get_data(i, "Mh")
                vector_v[i, :, :self.num_levels] = self.dataloader.get_data(
                    i, "Wh")

            # interpolate vertical value
            # blah

            # create field vector

            for l in range(self.num_levels-1):
                for i in range(self.num_points):
                    normal = spherical(1.0,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    r = radii[l]
                    vertex = spherical(r,
                                       lonlat[i, 0],
                                       lonlat[i, 1])
                    field[:, l, i, :] = vector_h[:, i, l, :]

                    # field[:, l, i, 0, :] = (0.0, 0.0, 0.0)
                    # field[:, l, i, 1, :] = normal
            data = np.zeros((self.num_samples, self.num_levels, self.num_points,
                             3), dtype=np.float32)

            if np.any(np.isnan(data)):
                print("NaN")

            # compute norm
            data = np.sqrt(np.sum(np.power(field, 2.0), axis=3))
            # data = data[:, :, :, 0]
            # rescale
            for i in range(self.num_samples):
                for j in range(self.num_levels):
                    d = data[i, j, :]
                    min_data = np.min(d)
                    max_data = np.max(d)
                    if np.any(max_data - min_data == 0):
                        print("zeros")
                    else:
                        data[i, j, :] = (d - min_data)/(max_data - min_data)

            # min_data = np.min(data)
            # max_data = np.max(data)

            # print(min_data, max_data)
            # zeros = (max_data - min_data) == 0
            # non_zeros = (max_data - min_data) != 0
            # data[non_zeros] = (data[non_zeros] - min_data[non_zeros]) / \
            #     (max_data[non_zeros] - min_data[non_zeros])
            # data[zeros] = 0.0
            self.data_scalar = data
            # self.data_color = data
            self.data_color[:, :, :, :] = np.array(
                H_to_RGB(360.0*data), copy=True, dtype=np.float32)

        elif dataname == "Momentum_norm_vert":
            vector_h = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels, 3), dtype=np.float32)
            vector_v = np.zeros((self.num_samples, self.num_points,
                                 self.num_levels), dtype=np.float32)
            print(self.num_samples)
            # get the data
            for i in range(self.num_samples):
                vector_h[i, :, :self.num_levels -
                         1] = self.dataloader.get_data(i, "Mh")
                vector_v[i, :, :self.num_levels] = self.dataloader.get_data(
                    i, "Wh")
            data = normalize(vector_v)
            self.data_scalar = np.array(np.swapaxes(data, 1, 2))
            self.data_color = np.array(np.swapaxes(
                H_to_RGB(360.0*data[:, :, :]), 1, 2), copy=True, dtype=np.float32)

        else:
            print("Unrecognised data name for scalar field: ", dataname)

    def get_color_data(self, dataset_idx):
        # print(self.data_color[dataset_idx, :, :, :])
        # return self.data_color[dataset_idx, :, :, :]
        return self.data_color[dataset_idx, :, :, :], self.data_scalar[dataset_idx, :, :]

    def get_field_data(self, idx):
        return self.field_data[idx, :, :, :, :]


class dataloader:
    def __init__(self, folder):
        self.folder = pathlib.Path(folder)

        # grid_dataset = folder / "esp_output_grid_Earth.h5"

        # #print("grid def")
        # grid = h5py.File(grid_dataset)
        # self.lonlat = grid['lonlat']
        # self.num_points = int(len(self.lonlat)/2)
        # # colors = np.zeros((num_points, 3), dtype=np.float32)

        # # grid color indexing
        # # indexing function through rhombis
        # # level
        # print("num points", self.num_points)
        # self.g = int(pow((self.num_points - 2)/10, 1/4)) - 2
        # print("level: ", self.g)
        # num_rhombi = 10
        # # nfaces
        # num_subrhombi = int(pow(4.0, self.g - 4))
        # nfaces = num_subrhombi
        # num_points_side_region = int(pow(2, 4))
        # nl_reg = num_points_side_region
        # nl2 = int(pow(num_points_side_region, 2))
        # kxl = int(sqrt(num_subrhombi))
        # print("nl_reg: ", nl_reg)
        # print("kxl: ", kxl)

        # def idx(fc, kx, ky, i, j):
        #     return nl2*(fc*nfaces + ky*kxl + kx) + j*nl_reg + i

        # read and sort all files
        datasets = folder.glob('esp_output_*_*.h5')

        dataset_re = re.compile('esp_output_(.*)_(\d+)')

        self.planets = {}

        # get dataset files
        for f in datasets:
            match = dataset_re.match(f.stem)
            if match is not None:
                basename = match.group(1)
                number = match.group(2)
                if basename not in self.planets:
                    self.planets[basename] = {'datasets': [(number, f)]}
                else:
                    self.planets[basename]['datasets'].append((number, f))

        # sort datasets files, get grid and planet file
        for planet, values in self.planets.items():
            # sort them numericaly
            d = values['datasets']
            values['datasets'] = [a[1]
                                  for a in sorted(d, key=lambda k:int(k[0]))]

            # open grid file
            grid = folder / ('esp_output_grid_{}.h5'.format(planet))
            self.planets[planet]['grid'] = grid
            # open planet file
            planet_def = folder / ('esp_output_planet_{}.h5'.format(planet))
            self.planets[planet]['def'] = planet_def

    def get_radius(self):
        return self.planet_def['A'][0]

    def get_altitudes(self):
        return self.grid['Altitude']

    def get_altitudesh(self):
        return self.grid['Altitudeh']

    def select_planet(self, planet):
        self.planet = planet
        self.grid = h5py.File(self.planets[self.planet]['grid'])
        self.planet_def = h5py.File(self.planets[self.planet]['def'])

        self.datasets = self.planets[self.planet]['datasets']
        self.num_samples = len(self.datasets)
        self.lonlat = self.grid['lonlat']
        self.num_points = len(self.lonlat)//2
        self.num_levels = int(self.grid['nv'][0])

    def get_planets(self):
        return list(self.planets.keys())

    def get_num_datasets(self):
        return self.num_samples

    def get_dim(self):
        return self.num_points, self.num_levels

    def get_grid(self):
        return np.array(self.lonlat).reshape(self.num_points, 2)

    def get_data(self, dataset_idx, data_name):
        print("loading: ", self.datasets[dataset_idx])
        data = None
        if (data_name == "Pressure"
                or data_name == "Rho"):
            data = np.array(h5py.File(self.datasets[dataset_idx])[data_name], dtype=np.float32).reshape(
                (self.num_points, self.num_levels))
        elif data_name == "Mh":
            data = np.array(h5py.File(self.datasets[dataset_idx])[data_name], dtype=np.float32).reshape(
                (self.num_points, self.num_levels, 3))
        elif data_name == "Wh":
            data = np.array(h5py.File(self.datasets[dataset_idx])[data_name], dtype=np.float32).reshape(
                (self.num_points, self.num_levels+1))
        return data

    def get_vector_data(self, dataset_idx, data_name):
        pass

        # grd = None
#         for planet, files in planets.items():
#             print("Planet: ", planet)
#             print("Number of datafiles:", len(files['datasets']))
#             num_samples = len(files['datasets'])
#             print("Planet def")
#             planet_def = h5py.File(files['def'])
#             for k, v in planet_def.items():
#                 print(k, v[0])

#             print("grid def")
#             grid = h5py.File(files['grid'])
#     for k, v in grid.items():
#         print(k, v)

#     num_levels = int(grid['nv'][0])
#     num_datas = 0
#     if len(files['datasets']) > 0:
#         print("datafiles def")
#         datafile0 = h5py.File(files['datasets'][0])
#         num_datas = len(datafile0['Pressure'])

#         for k, v in datafile0.items():
#             print(k, v)

#     if grd is None:
#         grd = grid['lonlat']
#         # neighbours = grid['pntloc']
#     num_points = len(grd)//2
#     ground_moment = np.zeros(
#         (len(files['datasets']),   num_points, num_levels, 3), dtype=np.float32)
#     print(ground_moment.shape)
#     for i in range(len(files['datasets'])):
#         d = h5py.File(files['datasets'][i])
#         ground_moment[i, :, :, :] = np.array(d['Mh']).reshape(num_points,
#                                                               num_levels,
#                                                               3)
#     ground_moment = np.transpose(ground_moment, (0, 2, 1, 3))
#     pressure = np.zeros(
#         (len(files['datasets']), num_datas//(num_levels)), dtype=np.float32)

#     for i in range(len(files['datasets'])):
#         d = h5py.File(files['datasets'][i])
#         pressure[i, :] = np.array(d['Pressure']).reshape(
#             num_points, num_levels)[:, 0]
#     print(pressure.shape)
# #    for n in neighbours:
# #        print(n)
# #    print(neighbours)

# colors2 = np.zeros((num_samples, num_levels, num_points, 3), dtype=np.float32)
# # colors[:, :, 0] = pressure/np.max(pressure)
# k = 0.2
# colors2[:, :, :, 0] = k + (1.0-k)*ground_moment[:, :, :, 0] / \
#     np.max(ground_moment[:, :, :, 0])
# colors2[:, :, :, 1] = k + (1.0-k)*ground_moment[:, :, :, 1] / \
#     np.max(ground_moment[:, :, :, 1])
# colors2[:, :, :, 2] = k + (1.0-k)*ground_moment[:, :, :, 2] / \
#     np.max(ground_moment[:, :, :, 2])
# # for i in range(n):
# #     colors[i, :, 0] = i/n

# # print(colors)

# print("color data size: ", colors2.nbytes/(1024*1024), "MiB")

# # moments = np.zeros((num_samples, len(grid['lonlat'])//2, 3), dtype=np.float32)
# # # colors[:, :, 0] = pressure/np.max(pressure)
# # k = 0.2
# # moments[:, :, :] = ground_moment[:, :, :]
