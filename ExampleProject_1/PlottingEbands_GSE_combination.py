import re
import codecs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from VaspWheels import GetElectronicBands, GetKpath, Visualization

###############################################################################################################
GBE = GetElectronicBands.vasp()  # 调用GetElectronicbands模块
GK = GetKpath.vasp()            # 调用GetKpath模块
VI = Visualization.plot()       # 调用Visualization模块
###############################################################################################################
# 定义一些可能用到的函数

# 此函数可以用于在列表中定位特定元素的序号引索
def GetIndex(element,target_list):
    return [index for (index,value) in enumerate(target_list) if value == element]

# 此函数可以用于从计算结果中提取费米面
def GetFermiEnergy(Markdown):  # Markdown文件中应记载着准确的费米能级
    pattern = re.compile(r'-?\d+\.?\d+')  # 匹配浮点数的正则表达式
    f = codecs.open(Markdown, 'rb', 'utf-8', 'ignore')
    line = f.readline()
    energy = pattern.findall(line)
    fermi_energy = float(energy[0])
    return fermi_energy

# 此函数可以用于费米面调零
def ShiftFermi(energy, fermi_energy):
    nbands = len(energy)  # number of bands
    nkpoints = len(energy[0])  # number of k points calculated
    shifted_energy = [[] for n in range(nbands)]
    for i in range(nbands):
        for j in range(nkpoints):
            shifted_energy[i].append(energy[i][j]-fermi_energy)
    return shifted_energy

# 此函数可以用于获取可用的能带数据
def Ebands(EIGENVAL,Kpath,fermi_energy,lattice=('HEX', [3.16, 3.16, 12.9, 90, 90, 120], 'primitive')):
    # 提取能带计算结果以及各种参数
    bands = GBE.GetData(EIGENVAL)
    nbands = bands['number']  # 提取能带总数
    nkpoints_total = bands['num kpoints']  # 提取K点总数
    nnodes = len(Kpath)  # 高对称点数
    npoints = int((nkpoints_total - nnodes) / (nnodes - 1))  # 两个高对称点中间的取点数 = (K点总数-高对称点数)/(高对称点数-1)
    energy = bands['energy']  # 能带具体的能量值
    energy = ShiftFermi(energy, fermi_energy)  # 进行费米面调零

    # 确定K点路径（X轴）
    k, knodes = GK.ProjectKpath(Kpath, npoints, LatticeCorrection='True', Lattice=lattice)  # 生成投影到一维的K点路径

    return k,energy,knodes

# 此函数可以用于能带截取
def InterceptingEbands(EIGENVAL,InterceptedKpath,fermi_energy):
    # 提取能带计算结果以及各种参数
    bands = GBE.GetData(EIGENVAL)
    nbands = bands['number']  # 提取能带总数
    kpath = bands['kpath']  # K点路径
    print(kpath)
    energy = bands['energy']  # 能带具体的能量值
    energy = ShiftFermi(energy, fermi_energy)  # 进行费米面调零

    starting_point, end_point = InterceptedKpath  # 从InterceptedKpath中提取起点跟终点
    starting_index = kpath.index(starting_point)
    #end_index = kpath.index(end_point)            # 终点的序号
    end_index = GetIndex(end_point,kpath)[1]

    print(starting_index,end_index)

    # 对能带进行切片
    InterceptedEbands = []
    for i in range(nbands):
        InterceptedEbands.append(energy[i][starting_index:end_index+1])

    K_projected = np.linspace(0,100,end_index-starting_index+1)

    return K_projected, InterceptedEbands, nbands

if __name__=='__main__':
    main_dir = 'D:/Projects/PhaseTransistor/Data/Simulation/GSE/OPTCELL/4_D3BJ_GSE_OPTCELL_1'  # 办公室电脑
    # main_dir = 'D:/PhD_research/办公室电脑/Data/Simulation/GSE/OPTCELL/4_D3BJ_GSE_OPTCELL_1'  # 宿舍电脑
    data_list = ['0.000', '0.100','0.200']
    Efield = ['0 V/nm', '0.5 V/nm', '1.0 V/nm']  # 电场强度

    # 确定晶格参数
    lattice = ('HEX', [3.16, 3.16, 12.9, 90, 90, 120], 'primitive')

    VI.GlobalSetting()  # 全局画图设定

    # fig = plt.figure(figsize=(6.9291,1.6))  # 控制图像大小
    fig = plt.figure(figsize=(7.7952,1.8))  # 控制图像大小
    # grid = plt.GridSpec(2,7,wspace=0)  # 创建柔性网格用于空间分配, wspace调整子图形之间的横向距离

    grid_shape = (3,17)
    num_row, num_col = grid_shape
    bands_fig_1 = plt.subplot2grid(grid_shape,(0,0),rowspan=num_row,colspan=5)
    bands_fig_2 = plt.subplot2grid(grid_shape,(0,5),rowspan=num_row,colspan=5)
    bands_fig_3 = plt.subplot2grid(grid_shape,(0,10),rowspan=num_row,colspan=5)
    DOS_fig = plt.subplot2grid(grid_shape,(0,15),rowspan=num_row,colspan=2)

    for n in [bands_fig_1,bands_fig_2,bands_fig_3,DOS_fig]:
        for m in ['top', 'bottom', 'left', 'right']:
            n.spines[m].set_linewidth(0.5)  # 设置图像边框粗细

    y_major_locator = MultipleLocator(1)
    y_minor_locator = MultipleLocator(0.5)
    bands_fig_1.yaxis.set_major_locator(y_major_locator)
    bands_fig_1.yaxis.set_minor_locator(y_minor_locator)
    bands_fig_2.yaxis.set_major_locator(y_major_locator)
    bands_fig_2.yaxis.set_minor_locator(y_minor_locator)
    bands_fig_3.yaxis.set_major_locator(y_major_locator)
    bands_fig_3.yaxis.set_minor_locator(y_minor_locator)

    plt.subplots_adjust(wspace=0)

    # 画图参数
    fontsize, linewidth = [6,0.5]

    ylim=(-2.1,2.1)
    ymin,ymax = ylim
    # color = [np.array([124,172,247])/255.0,np.array([128,138,248])/255.0,
             #np.array([157,132,249])/255.0,np.array([195,137,250])/255.0]
    # color = [np.array([77,133,189])/255.0,np.array([247,144,61])/255.0,np.array([89,169,90])/255.0]
    color = [VI.CMYK_to_RGB(40,30,8,8),VI.CMYK_to_RGB(60,50,28,8),VI.CMYK_to_RGB(80,70,48,8)]

    grey = np.array([155,165,160])/255.0


    main_Kpath = [[r'$\Gamma$', 'M', 'K', r'$\Gamma$'], [[0, 0, 0], [0.5, 0, 0], [1.0 / 3.0, 1.0 / 3.0, 0], [0, 0, 0]]]


    dos_x, dos_y = [[],[]]  # 用于存放DOS数据的变量
    # 通过循环批量画能带子图并调节子图参数
    bands_fig_list = [bands_fig_1, bands_fig_2, bands_fig_3]
    for n in data_list:
        bands_index = data_list.index(n)  # 获取元素n在data_list中的引索序号，此序号可以用于跟子图编号一一对应

        EIGENVAL = main_dir+'/'+n+'/EIGENVAL'
        Markdown = main_dir+'/'+n+'/Markdown_SCF'

        E_fermi = GetFermiEnergy(Markdown)  # 提取费米能级


        # 利用这个循环提取DOS数据
        DOSCAR = main_dir+'/'+n+'/DOSCAR'
        DOS_data = GBE.GetData(DOSCAR)
        dos_x.append(DOS_data['DOS'])
        dos_energy = DOS_data['energy']
        dos_y.append(GBE.ShiftFermiSurface(dos_energy, E_fermi))
        ###########################################

        # x,y,x_nodes = Ebands(EIGENVAL,main_Kpath[1],E_fermi)
        # 提取能带计算结果以及各种参数
        bands_dict = GBE.GetEbands(EIGENVAL)
        num_bands = bands_dict['num_bands']  # 提取能带总数
        num_kpoints = bands_dict['num_kpoints']  # 提取K点总数
        Kpath = bands_dict['Kpath']  # K点路径
        bands = bands_dict['bands']  # 能带具体的能量值

        # 生成投影到一维的K点路径
        num_segments = 3
        Kpath_projected, Knodes_projected = GK.ProjectKpath(Kpath, num_segments, LatticeCorrection='True',Lattice=lattice)

        # 费米面调零
        Eg, Ev_max, Ec_min, extremum_location = GBE.GetBandgap(EIGENVAL, mode='occupation')

        bands_shifted = GBE.ShiftFermiSurface(bands, Ev_max)

        # 正式开始画图
        bands_plot = bands_fig_list[bands_index]  # 确定此循环中要画的子图

        for i in range(len(bands_shifted)):
            bands_plot.plot(Kpath_projected,bands_shifted[i],color=color[bands_index],linewidth=linewidth)

        bands_plot.set_xlim(min(Kpath_projected),max(Kpath_projected))
        bands_plot.set_ylim(ymin,ymax)

        bands_plot.set_title(r'$\mathcal{E}$ = '+Efield[bands_index],size=fontsize, pad=5)  # 设置标题
        bands_plot.tick_params(labelsize=fontsize)  # 更改刻度大小
        bands_plot.tick_params(axis='x', color='w')  # 更改x轴刻度颜色
        ################################################################################################
        # 设置主次刻度线长度
        plt.tick_params(which='major', length=3, width=0.5)  # 设置主刻度长度
        bands_plot.tick_params(which='minor', length=1.5, width=0.5)  # 设置次刻度长度
        ################################################################################################

        bands_plot.hlines(0,min(Kpath_projected),max(Kpath_projected),linewidth=linewidth, linestyles='dashed', colors=grey)  # 费米能的位置
        for i in range(1,len(Knodes_projected) - 1):  # 第一个高对称点与图的左边界重合，最后一个跟右边界重合，所以不必作分割线
            bands_plot.vlines(Knodes_projected[i], ymin, ymax, linewidth=linewidth, linestyles='dashed', colors=grey)

        # 对个别子图进行调整
        if bands_index >=1:
            bands_plot.set_yticklabels([])
            nodes = [Knodes_projected[i+1] for i in range(3)]
            path = [main_Kpath[0][i+1] for i in range(3)]
        else:
            # bands_plot.set_ylabel('$E-E_{f}$ (eV)',size=20)
            bands_plot.set_ylabel('Energy (eV)', size=fontsize)
            nodes = Knodes_projected
            path = main_Kpath[0]
        bands_plot.set_xticks(nodes)
        bands_plot.set_xticklabels(path,size=fontsize)

    bands_fig_1.text(2.2,0.6,'$\Lambda_\mathrm{min}$',size=fontsize)

    # 画DOS
    DOS_fig.yaxis.set_major_locator(y_major_locator)
    DOS_fig.yaxis.set_minor_locator(y_minor_locator)

    dos_x_min,dos_x_max = (0,20)
    for i in range(3):
        DOS_fig.plot(dos_x[i],dos_y[i],color=color[i],label=r'$\mathcal{E}$ = '+Efield[i],linewidth=linewidth)
    DOS_fig.hlines(0,dos_x_min,dos_x_max,linewidth=linewidth, linestyles='dashed', colors=grey)
    # DOS_fig.legend(loc=(0.1,0.3),fontsize=12,frameon=False)
    DOS_fig.set_yticklabels([])
    DOS_fig.set_xlim(dos_x_min,dos_x_max)
    DOS_fig.set_ylim(ymin,ymax)
    DOS_fig.set_xticks([10])
    #Text properties for the labels. These take effect only if you pass labels. In other cases, please use tick_params.
    DOS_fig.tick_params(axis='x',color='w')  # 对于subplot，要调整刻度样式的话，需要采用tick_params函数
    DOS_fig.set_xticklabels(['DOS (a.u.)'],size=fontsize)

    ################################################################################################
    # 设置主次刻度线长度
    plt.tick_params(which='major', length=3, width=0.5)  # 设置主刻度长度
    plt.tick_params(which='minor', length=1.5, width=0.5)  # 设置次刻度长度
    ################################################################################################
    plt.show()

    # saving_directory = 'D:/Projects/PhaseTransistor/Data/Figures/Giant Stark Effect (GSE)/能带组图/'  # 办公室电脑
    # saving_directory = 'D:/Projects/PhaseTransistor/Data/Figures/CarrierTransportation/Version_22.12.30/'  # 办公室电脑
    saving_directory = 'D:/Projects/PhaseTransistor/Gallery/Figures/All/'  # 办公室电脑汇总
    VI.SavingFigure(saving_directory, filename='GSE_combination', format='pdf')
    VI.SavingFigure(saving_directory, filename='GSE_combination', format='eps')