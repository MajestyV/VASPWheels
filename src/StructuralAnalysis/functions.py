import itertools as it  # 引入itertools模块, 用于生成排列和组合
import numpy as np
import matplotlib.pyplot as plt

# 杂记 ------------------------------------------------------------------------------------------------------------------
# 1. Radial Distribution Function (RDF) 和 Pair Correlation Function (PCF) 的区别：
# 定义其实没有任何区别，但是在此代码中，RDF和PCF是通过不同的算法计算的，RDF是计算原子间距离，而PCF是通过原子对排列组合计算键长实现的
# ----------------------------------------------------------------------------------------------------------------------

def RDF(atom_pos: np.ndarray, r_cutoff: float, dr: float):
    '''
    Calculate the Radial Distribution Functions (RDFs) of the target systems from the atomic position list
    :param atom_pos: the coordinates of the atomic positions
    :param r_cutoff: the cutoff radius
    :param dr: the resolution of the RDF
    :return: the Radial Distribution Functions (RDFs)
    '''
    num_atoms = len(atom_pos)  # 体系内的原子个数
    r = np.arange(0, int(r_cutoff/dr))*dr  # 生成RDF的横坐标
    RDF = np.zeros(int(r_cutoff/dr), dtype=float)  # 初始化RDF数组

    # 历遍计算体系中所有键的长度
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):  # 避免重复计算
            r_ij = np.linalg.norm(atom_pos[i] - atom_pos[j])
            if r_ij < r_cutoff:
                bin_index = int(r_ij/dr)  # 计算键长所在的bin
                RDF[bin_index] += 1  # 累加（Count each pair only once）

    # Normalization
    dV = 4.0/3.0*np.pi*(np.arange(1, len(RDF)+1)*dr)**3  # 计算每个bin的体积微元
    RDF /= (num_atoms*(num_atoms-1)/2.0)*dV  # 计算RDF

    return r, RDF

def PCF(atom_pos: np.ndarray, r_cutoff: float, dr: float):
    '''
    Calculate the Pair Correlation Functions (PCFs) of the target systems from the atomic position list
    :param atom_pos: the coordinates of the atomic positions
    :param r_cutoff: the cutoff radius
    :param dr: the resolution of the PCF
    :return: the Pair Correlation Functions (PCFs)
    '''

    atomic_pairs = list(it.combinations(atom_pos, 2))  # 利用itertools模块生成原子对（其实就是原子坐标列表元素的组合）
    bonds = np.array([np.linalg.norm(pair[0] - pair[1]) for pair in atomic_pairs])  # 计算原子对之间的距离（键长）

    PCF, bin = np.histogram(bonds, bins=int(r_cutoff/dr), range=(0,r_cutoff))
    r = (bin[:-1] + bin[1:])/2.0   # 计算PCF的横坐标

    # Normalization
    num_atoms = len(atom_pos)  # 体系内的原子个数
    PCF = np.array(PCF, dtype=float)  # 将PCF转换为浮点数数组
    dV = 4.0 / 3.0 * np.pi * (np.arange(1, len(PCF) + 1) * dr) ** 3  # 计算每个bin的体积微元
    PCF /= (num_atoms * (num_atoms - 1) / 2.0) * dV  # 计算RDF

    return r, PCF

if __name__ == '__main__':
    # example_datafile = 'D:/Projects/PINN/Data/Si_melting_example/traj_nvt.lammpstrj'  # 办公室电脑
    # example_datafile = 'D:/Projects/EchoTrajNetwork/Test/traj_nvt.lammpstrj'  # 宿舍电脑
    # example_datafile = '/Users/liusongwei/Desktop/traj_nvt.lammpstrj'  # Macbook
    example_datafile = 'D:/PhD_research/EchoTrajNetwork/Test/traj_nvt.lammpstrj'  # 家里电脑

    MD_data = GetTraj_LAMMPS(example_datafile)

    # data_reorganized = MD_data.Reoraganized()
    # print(data_reorganized[0,27])

    pos_traj = MD_data.ExtractTrajectoryData((2, 5))
    # print(pos_traj[0][27])

    frame = 500

    atoms = pos_traj[frame]

    # 快速可视化，检验数据
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(pos_traj[frame, :, 0], pos_traj[frame, :, 1], pos_traj[frame, :, 2])
    #plt.show(block=True)

    # RDF测试 - 1

    # radius_search, RDF = SA.RDF(atoms, 20, 500)

    # plt.plot(radius_search, RDF)

    # RDF测试 - 2

    r, RDF = RDF(atoms, 10, 0.05)

    plt.plot(r, RDF)

    # PCF测试
    r, PCF = PCF(atoms, 10, 0.05)

    print(r)

    plt.plot(r, PCF)

    # 非console环境下，显示图片的命令
    plt.show(block=True)