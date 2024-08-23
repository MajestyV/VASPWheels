import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
from src import GetStructure, Supercell
# from mda.mdio import write, read

if __name__ == '__main__':
    POSCAR = 'D:/PhD_research/Jingfang Pei/Solution-processed IC/Simulation/relaxed_structure.vasp'  # Lingjiang

    GS = GetStructure()  # 实例化GetStructure类

    structure_info = GS.GetPOSCAR(POSCAR)  # 读取POSCAR文件

    # print(structure_info['atom'])  # 输出原子种类信息

    atom_pos_Si = structure_info['atom_pos_Si']  # 提取原子坐标信息
    atom_pos_O = structure_info['atom_pos_O']  # 提取原子坐标信息

    supercell = (3,3,3)

    atom_pos_Si = Supercell(atom_pos_Si, structure_info['lattice_vector'], supercell)
    atom_pos_O = Supercell(atom_pos_O, structure_info['lattice_vector'], supercell)

    pair_dist = distance.pdist(atom_pos_O)
    # pair_dist = distance.cdist(atom_pos_Si,atom_pos_O)  # 计算原子间的距离
    # dists = distance.squareform(dists)  # 将距离向量转换为距离矩阵
    # print(dists.reshape(-1))
    pair_dist = pair_dist.reshape(-1)

    # 读取分子坐标文件
    # traj = read('trajectory.xtc')
    # 提取每个分子的坐标
    # atom_coords = traj.xyz[0]
    # 计算所有原子之间的距离
    # dists = pairwise_distances(atom_pos_O)
    # print(dists)
    # 计算RDF
    # rdf = np.histogram(dists, bins=200, range=[0, 10])[0]

    r_cutoff = 10
    dr = 0.1

    RDF, bin = np.histogram(pair_dist, bins=int(r_cutoff / dr), range=(0, r_cutoff))
    r = (bin[:-1] + bin[1:]) / 2.0  # 计算PCF的横坐标

    np.savetxt('D:/PhD_research/Jingfang Pei/Solution-processed IC/Simulation/RDF_O.txt', np.array([r, RDF]).T)

    plt.plot(r, RDF)

    plt.show(block=True)

    # rdf = np.histogram2d(np.argwhere(dists < 10), bins=50, range=[[0, trj.n_atoms]])[0]