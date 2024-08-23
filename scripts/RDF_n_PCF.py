import matplotlib.pyplot as plt
from src import RDF, PCF, GetStructure

if __name__ == '__main__':
    POSCAR = 'D:/PhD_research/Jingfang Pei/Solution-processed IC/Simulation/relaxed_structure.vasp'  # Lingjiang

    GS = GetStructure()  # 实例化GetStructure类

    structure_info = GS.GetPOSCAR(POSCAR)  # 读取POSCAR文件

    # print(structure_info['atom'])  # 输出原子种类信息

    atom_pos_O = structure_info['atom_pos_O']  # 提取原子坐标信息
    # print(atom_pos_O)

    # RDF测试 - 1

    # radius_search, RDF = SA.RDF(atoms, 20, 500)

    # plt.plot(radius_search, RDF)

    # RDF测试 - 2

    # r, RDF = RDF(atom_pos_O, 10, 0.05)

    # plt.plot(r, RDF)

    # PCF测试
    r, PCF = PCF(atom_pos_O, 10, 0.05)
    print(r)

    # print(r)

    # plt.plot(r, PCF)

    # 非console环境下，显示图片的命令
    plt.show(block=True)