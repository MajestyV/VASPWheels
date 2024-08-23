import numpy as np

def Supercell(atom_pos: np.ndarray, lattice_vector: np.ndarray, supercell_shape: tuple) -> np.ndarray:
    '''
    此函数用于生成超胞结构
    :param atom_pos:
    :param lattice_vector:
    :param supercell_shape:
    :return:
    '''

    num_atom = len(atom_pos)  # 原子个数
    num_cell = np.prod(supercell_shape)  # 超胞中的晶胞数（np.prod可以实现连乘功能）

    nx, ny, nz = supercell_shape  # 超胞的形状

    a1, a2, a3 = lattice_vector  # 提取晶格向量

    atom_pos_supercell = np.zeros((num_atom*num_cell, 3), dtype=float)  # 初始化超胞中的原子坐标数组
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for n in range(num_atom):
                    atom_pos_supercell[n + num_atom*(i*ny*nz + j*nz + k)] = atom_pos[n] + i*a1 + j*a2 + k*a3

    return atom_pos_supercell

if __name__ == '__main__':
    pass