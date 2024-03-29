import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from VaspWheels import GetKpath,GetElectronicBands,GeneralAnalyzer

GK = GetKpath.vasp()              # 调用GetKpath模块
GEB = GetElectronicBands.vasp()   # 调用GetElectronicBands模块
GA = GeneralAnalyzer.functions()  # 调用GeneralAnalyzer模块

##################################################################################################################
# 先定义一些在接下来的计算中可能用到的参数

# 晶体结构信息

crystal_info = ['HEX', [3.1473295667554400, 3.1473295667554400, 43.9122903625234997, 90, 90, 120], 'primitive']

real_lattice = [[3.1473295667554400, 0.0000000000000000, 0.0000000000000000],
                [-1.5736647833777198, 2.7256673589277995, 0.0000000000000000],
                [0.0000000000000000, 0.0000000000000000, 43.9122903625234997]]

reciprocal_lattice = GK.CalculateReciprocalLattice(real_lattice)

# 一些关键的高对称点（High symmetry point, HSP）
HSP = {'G':  [0, 0, 0],
       'M':  [1/2.0,0,0],
       'M1': [1/2.0,-1/2.0,0],
       'K':  [1/3.0, 1/3.0, 0],
       'Sl': [0.26500000000050006,0,0],
       'S1': [1/3.0, 0, 0],
       'S2': [2/3.0, 0, 0],
       'L':  [0.176666666667, 0.176666666667, 0]}

##################################################################################################################
# 生成计算迁移率所需的K点文件

Kpath_directory = 'D:/PhD_research/Data/Simulation/MoS2/CarrierTransport/4/EffectiveMass/Kpath/'  # 宿舍电脑
# Kpath_directory = 'D:/Projects/PhaseTransistor/Data/Simulation/Conductivity/EffectiveMass/Kpath/'  # 办公室电脑

origin =      [HSP['K'], HSP['K'], HSP['K'], HSP['K'], HSP['G'], HSP['G'], HSP['G'], HSP['L'], HSP['L'], HSP['L']]
destination = [HSP['G'], HSP['M'], HSP['S1'],HSP['S2'],HSP['K'], HSP['M'], HSP['M1'],HSP['K'], HSP['G'], HSP['Sl']]

Kpoints_list = GK.GenerateKpath_segment(origin,destination,30,0.01,reciprocal_lattice)

#GK.GenKPOINTS(Kpath_directory+'K-path_EffectiveMass_22.12.16',Kpoints_list)  # 生成K点路径文件
##################################################################################################################
# 计算电子跟空穴的有效质量
# data_directory = 'D:/Projects/PhaseTransistor/Data/Simulation/CarrierTransport/4/EffectiveMass/Data/'  # 办公室电脑
data_directory = 'D:/PhD_research/Data/Simulation/MoS2/CarrierTransport/4/EffectiveMass/Data/'  # 宿舍电脑

positive_field = '4_PositiveField_EffectiveMass_22.12.16'
negative_field = '4_NegativeField_EffectiveMass_22.12.16'
total = 'EffectiveMass_Total'

Efield_file = ['m0.300','m0.275','m0.250','m0.225','m0.200','m0.175','m0.150','m0.125','m0.100','m0.075','m0.050','m0.025',
               '0.000','0.025','0.050','0.075','0.100','0.125','0.150','0.175','0.200','0.225','0.250','0.275','0.300']

##################################################################################################################
# 检查数据
E_test = 'EffectiveMass_ZeroField_noSOC'
EIGENVAL = data_directory+total+'/'+E_test+'/'+'EIGENVAL'
valence_band, conduction_band = GEB.GetBandEdges(EIGENVAL)
print(conduction_band)

# k = np.array([i for i in range(len(valence_band))])

effective_mass = GA.CalculateEffectiveMass(0.01,conduction_band,10,points_evaluating=12)  # 计算有效质量，取点个数很关键
print(effective_mass)

#plt.plot(k,valence_band)
#plt.plot(k,conduction_band)
#plt.xlim(0,29)
#plt.ylim(1.42,1.44)
##################################################################################################################

# 计算曲率时的取点方向
# 计算曲率时的取点方向
direction = ['K-G','K-M','K-S1','K-S2',
             'G-K','G-M','G-M1',
             'L-K','L-G','L-Sl']

# 批量计算有效质量
EffectiveMass_data = np.zeros((len(Efield_file),len(direction)))  # 用于存放有效质量的矩阵
# Efield_positive = [str(float(Efield_file[i])*5)+' V/nm' for i in range(len(Efield_file))]  # 正电场强度列表
# Efield_negative = [str(-float(Efield_file[i])*5)+' V/nm' for i in range(len(Efield_file))]  # 负电场强度列表
field = np.linspace(-1.5,1.5,25)  # 获取电场强度序列
# print(field)
Efield = [str(float(field[i]))+' V/nm' for i in range(len(Efield_file))]  # 电场强度列表
for i in range(len(Efield_file)):
    # EIGENVAL = data_directory+positive_field+'/'+Efield_file[i]+'/'+'EIGENVAL'  # 正电场数据文件
    # EIGENVAL = data_directory+negative_field+'/m'+Efield_file[i]+'/'+'EIGENVAL'  # 负电场数据文件
    EIGENVAL = data_directory + total + '/' + Efield_file[i] + '/' + 'EIGENVAL'  # 总数据文件
    valence_band, conduction_band = GEB.GetBandEdges(EIGENVAL)  # 提取导带跟价带

    effective_mass = GA.CalculateEffectiveMass(0.01,valence_band,10,points_evaluating=9)  # 计算有效质量

    EffectiveMass_data[i] = effective_mass

#GA.SaveData('C:/Users/13682/OneDrive/桌面/Test/',EffectiveMass_data,
             #file_name='Hole_EffectiveMass_NegativeField',col_index=direction,row_index=Efield_positive)
#GA.SaveData('C:/Users/13682/OneDrive/桌面/Test/',EffectiveMass_data,
            #filename='1',col_index=direction,row_index=Efield)

# saving_directory = 'D:/Projects/PhaseTransistor/Data/Simulation/CarrierTransport/4/EffectiveMass/Data/Summary/9 points_step 0.01/'
# GA.SaveData(saving_directory,EffectiveMass_data,filename='Hole_EffectiveMass_SOC',col_index=direction,row_index=Efield)