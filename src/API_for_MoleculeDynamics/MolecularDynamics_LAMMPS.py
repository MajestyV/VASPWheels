# This code is designed for extracting the molecular dynamics data from LAMMPS trajectory files.

import linecache
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class GetTraj_LAMMPS:
    """
    This class is designed for extracting the molecular dynamics data from LAMMPS trajectory files.
    """
    def __init__(self, data_file, degree_of_freedom=3):
        system_info = [[linecache.getline(data_file, i + 1)] for i in range(9)]  # 读取目标系统的信息，应注意lincache函数默认从1开始

        num_atom = int(system_info[3][0])  # 读取目标系统的原子个数，即为系统信息的第四行

        nline_per_step = num_atom + 9  # 每一个分块的行数

        # 统计总行数
        with open(data_file) as file:
            for num_line, _ in enumerate(file, 1):
                pass

        num_step = int(num_line/nline_per_step)  # 计算总共的采样步数
        time = np.array([float(linecache.getline(data_file, 2 + i * nline_per_step)) for i in range(num_step)])  # 时间步

        rows_to_skip = []  # 跳过一些冗余的信息行
        for i in range(num_step):
            for j in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
                rows_to_skip.append(i * nline_per_step + j)

        # 利用pandas读取molecular dynamic轨迹
        data_chunks = pd.read_csv(data_file, header=None, skiprows=rows_to_skip, sep='\s+', chunksize=num_atom)  # 分块读取
        # 数据整理及转化
        data_DataFrame = [chunk for chunk in data_chunks]  # DataFrame格式的数据
        data_list = [data_DataFrame[i].values for i in range(num_step)]  # 提取具体的分子运动轨迹数据

        # 将一些变量转换为实例变量，方便后续处理
        self.deg_free = degree_of_freedom  # 空间自由度
        self.num_atom = num_atom    # 原子个数
        self.num_step = num_step    # 总步数
        self.time = time            # 时间步
        self.data_list = np.array(data_list)  # 分子运动轨迹数据（转换为numpy数组，方便后续运算处理）
        self.num_feature = len(data_list[0][0])  # 轨迹的特征数

    def Reoraganized(self):
        # 创建一个空数组用于存放整理后的数据
        data_reorganized = np.zeros((self.num_step, self.num_atom, self.num_feature), dtype=float)

        for i in range(self.num_step):
            for j in range(self.num_atom):
                atom_index = int(self.data_list[i][j][0])  # 提取原子序号原子序号
                data_reorganized[i,atom_index-1] = self.data_list[i,j]  # 重新排列数据

        return data_reorganized

    # 此函数可以提取体系特定特征的轨迹，默认的2:5是原子的位置信息
    def ExtractTrajectoryData(self,data_range=(2,5)):
        start, end = data_range
        data_reorganized = self.Reoraganized()
        return data_reorganized[:,:,start:end]

if __name__ == '__main__':
    example_datafile = 'D:/Projects/PINN/Data/Si_melting_example/traj_nvt.lammpstrj'

    MD_data = GetTraj_LAMMPS(example_datafile)

    # data_reorganized = MD_data.Reoraganized()
    # print(data_reorganized[0,27])

    pos_traj = MD_data.ExtractTrajectoryData((2,5))
    print(pos_traj[0][27])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(pos_traj[0,:,0], pos_traj[0,:,1], pos_traj[0,:,2])
    plt.show(block=True)
