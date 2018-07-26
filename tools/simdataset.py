
import numpy as np
import h5py
import pathlib
import re

from math import sqrt


class simdataset:
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

    def get_scalar_data(self, dataset_idx, data_name):
        data = np.array(h5py.File(self.datasets[dataset_idx])[data_name], dtype=np.float32).reshape(
            (self.num_points, self.num_levels))

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
