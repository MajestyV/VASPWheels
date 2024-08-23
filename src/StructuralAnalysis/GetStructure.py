# 此代码专门用于提取VASP可读的晶体结构文件POSCAR或者CONTCAR中的数据，或者根据指定构型生成POSCAR文件用于第一性原理计算

import numpy as np
import pandas as pd
import codecs, linecache

class GetStructure:
    """
    结构数据读取类
    """
    def __init__(self):
        self.name = 'GetStructure'

    # 一些后续可能会重复调用的通用函数，在一开头定义，方便后续代码的编写跟维护
    # 这个函数利用linecache模块，可以从数据文件中读出指定行的信息，并以字符串形式返回
    # 应注意，这个函数的行数从1开始，即line_index=5指要读文件中的第五行
    def GrepLineContent(self, file, line_index): return linecache.getline(file,line_index).strip()

    def GetPOSCAR_old(self, POSCAR: str) -> dict:
        '''
        该函数可以从POSCAR（CONTCAR）文件中提取晶体结构信息
        This function can extract the information of the lattice structure from the POSCAR (CONTCAR) file for V.A.S.P..
        param POSCAR: absolute address of the POSCAR file
        return: dict including structural info
        '''

        # 按行读取POSCAR文件
        file = codecs.open(POSCAR, 'rb', 'utf-8', 'ignore')
        line = file.readline()
        lindex = 0
        lattice_vector = []
        atomic_position_raw = []
        while line:
            content = line.split()  # 以空字符（空格，换行'\n'，制表符'\t'等）为分隔符对本行内容做切片 （什么都不填默认为空字符）
            if lindex == 0:
                system = line.split()
            elif lindex == 1:
                scale = float(content[0])
            elif lindex >= 2 and lindex <= 4:
                content = list(map(float, content))
                lattice_vector.append(content)
            elif lindex == 5:
                atom_species = content
            elif lindex == 6:
                num_atom = list(map(float, content))
            elif lindex == 7:
                coordinate = content[0]
            else:
                content = list(map(float, content))
                atomic_position_raw.append(content)

            line = file.readline()
            lindex += 1
        file.close()

        # 重整原子坐标信息，丢弃无意义的部分
        num_total = int(sum(num_atom))  # 将各原子数加起来得到原子总数,记得要替换为整型
        atomic_position = []
        for i in range(num_total):
            atomic_position.append(atomic_position_raw[i])

        structure_info = {'system': system,
                          'scale': scale,
                          'atom': atom_species,
                          'num_atom': num_atom,
                          'coordinate_system': coordinate,
                          'lattice_vector': lattice_vector,
                          'atom_pos': np.array(atomic_position)}  # 将原子坐标信息转换为数组，以防出错
        return structure_info

    def GetPOSCAR(self,POSCAR: str) -> dict:
        '''
        该函数可以从POSCAR（CONTCAR）文件中提取晶体结构信息
        This function can extract the information of the lattice structure from the POSCAR (CONTCAR) file for V.A.S.P..
        param POSCAR: absolute address of the POSCAR file
        return: dict including structural info
        '''
        system_name = self.GrepLineContent(POSCAR, 1)  # 系统名称

        scaling = self.GrepLineContent(POSCAR, 2)  # 缩放尺度

        lattice_vector = np.array([list(map(float,self.GrepLineContent(POSCAR, 3).split())),  # 晶格向量
                                   list(map(float,self.GrepLineContent(POSCAR, 4).split())),
                                   list(map(float,self.GrepLineContent(POSCAR, 5).split()))])

        # print(list(map(float,self.GrepLineContent(POSCAR, 3).split())))

        atom_species = self.GrepLineContent(POSCAR, 6).split()  # 原子种类
        num_atom = self.GrepLineContent(POSCAR, 7).split()  # 每个原子种类分别的原子数
        num_atom_tot = sum([int(i) for i in num_atom])  # 总原子数
        num_atom_species = len(atom_species)  # 原子种类的数目

        coordinate_system = self.GrepLineContent(POSCAR, 8)  # 坐标系

        atom_pos_DF = pd.read_csv(POSCAR, header=None, skiprows=8, sep='\s+', nrows=num_atom_tot)  # 分块读取
        atom_pos = atom_pos_DF.values  # DataFrame格式的数据转换为array格式
        for i in range(num_atom_tot):
            # print(atom_pos[i].reshape(1,-1))
            atom_pos_real = np.dot(atom_pos[i].reshape(1,-1),lattice_vector)  # 将原子分数坐标转换为实空间坐标
            # print(atom_pos_real)
            atom_pos[i] = atom_pos_real[0]
        # atom_pos = np.array([atom_pos[i]*lattice_vector for i in range(num_atom_tot)])  # 将原子分数坐标转换为实空间坐标

        structure_info = dict()  # 创建结构信息字典用于保存数据
        # 基本信息
        structure_info['system_name'] = system_name
        structure_info['scaling'] = scaling
        structure_info['lattice_vector'] = lattice_vector
        structure_info['atom_species'] = atom_species
        structure_info['num_atom'] = num_atom
        structure_info['num_atom_tot'] = num_atom_tot
        structure_info['num_atom_species'] = num_atom_species
        structure_info['coordinate_system'] = coordinate_system
        structure_info['atom_pos'] = atom_pos
        # 对于非单质而言，还需要分别提取每个原子种类的坐标信息

        for i in range(num_atom_species):
            if i == 0:
                # print(int(num_atom[i]))
                structure_info['atom_pos_' + atom_species[i]] = atom_pos[:int(num_atom[i])]
            else:
                structure_info['atom_pos_' + atom_species[i]] = atom_pos[
                                                                int(num_atom[i - 1]):int(num_atom[i - 1]) + int(
                                                                    num_atom[i])]
        return structure_info

if __name__ == '__main__':
    POSCAR = 'D:/PhD_research/Jingfang Pei/Solution-processed IC/Simulation/relaxed_structure.vasp'
    # structure_info = GetStructure().GetPOSCAR(POSCAR)

    GS = GetStructure()  # 实例化GetStructure类

    a = GS.GrepLineContent(POSCAR, 1)

    system_name = GS.GrepLineContent(POSCAR, 1)  # 系统名称

    scaling = GS.GrepLineContent(POSCAR, 2)  # 缩放尺度

    lattice_vector = np.array([GS.GrepLineContent(POSCAR, 3), # 晶格向量
                               GS.GrepLineContent(POSCAR, 4),
                               GS.GrepLineContent(POSCAR, 5)])

    atom_species = GS.GrepLineContent(POSCAR, 6).split()  # 原子种类
    num_atom = GS.GrepLineContent(POSCAR, 7).split()  # 每个原子种类分别的原子数
    num_atom_tot = sum([int(i) for i in num_atom])  # 总原子数
    num_atom_species = len(atom_species)  # 原子种类的数目

    coordinate_system = GS.GrepLineContent(POSCAR, 8)  # 坐标系

    atom_pos_DF = pd.read_csv(POSCAR, header=None, skiprows=8, sep='\s+', nrows=num_atom_tot)  # 分块读取
    atom_pos = atom_pos_DF.values  # DataFrame格式的数据转换为array格式

    structure_info = dict()  # 创建结构信息字典用于保存数据
    # 基本信息
    structure_info['system_name'] = system_name
    structure_info['scaling'] = scaling
    structure_info['lattice_vector'] = lattice_vector
    structure_info['atom_species'] = atom_species
    structure_info['num_atom'] = num_atom
    structure_info['num_atom_tot'] = num_atom_tot
    structure_info['num_atom_species'] = num_atom_species
    structure_info['coordinate_system'] = coordinate_system
    structure_info['atom_pos'] = atom_pos
    # 对于非单质而言，还需要分别提取每个原子种类的坐标信息

    for i in range(num_atom_species):
        if i == 0:
            # print(int(num_atom[i]))
            structure_info['atom_pos_'+atom_species[i]] = atom_pos[:int(num_atom[i])]
        else:
            structure_info['atom_pos_'+atom_species[i]] = atom_pos[int(num_atom[i-1]):int(num_atom[i-1])+int(num_atom[i])]



    print(structure_info['atom_pos_O'].shape)
    print(structure_info['atom_pos_Si'].shape)
    # for i in range(num_atom_species):




    # print(atom_pos)