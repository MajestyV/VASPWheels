import numpy as np
import VaspWheels as vw

if __name__=='__main__':
    control_index = 0  # 控制字

    # JCPGH1
    data_directory = 'D:/Projects/OptoTransition/Data/Homo-structure/Pentalayer/Test/0.0 V-nm/'

    layer_list = ['layer1', 'layer2', 'layer3', 'layer4', 'layer5']

    # 保存目录
    # saving_directory = 'D:/OneDrive/OneDrive - The Chinese University of Hong Kong/Desktop/DataFig_OptoTrans/OrbitalAnalysis'  # MMW502
    # saving_directory = 'D:/Projects/OptoTransition/Data/Figures/临时数据文件夹'  # JCPGH1

    saving_filename = layer_list[control_index]

    # 数据地址
    data_file = data_directory+'/'+layer_list[control_index]+'/PBAND_SUM_SOC.dat'

    # Fermi_factor = 0.20  # 费米面调零参数
    Fermi_factor = 0.185

    num_segments = 3  # 2D

    Kpath, Kpath_nodes = vw.API_vaspkit.GetProjectedKpath(data_file, num_segment=num_segments)  # 获取K空间轨迹的一维投影
    print(Kpath_nodes)

    K_range = (1.8184,3.14956)

    # 获取能带数据
    x_band, y_band, w_band = vw.API_vaspkit.GetProjectedBands(data_file,'tot',Fermi_adjust=Fermi_factor)

    # cmap = vw.colormap.iColarmap['Coolwarm']
    # cmap = vw.colormap.iColarmap['Blue_n_Red']
    # cmap = vw.colormap.iColarmap['Purple_n_Green']
    # cmap = 'seismic'  # 采用matplotlib标准色谱
    cmap = 'viridis'

    # 画图模块
    vw.VisualElectronic_vaspkit.VisualizePartialProjectedBands(x_band,y_band,w_band,K_range=K_range,
                                                         colormap=cmap,
                                                         colormap_norm=(0,0.4),size_band=np.abs(w_band)*5,
                                                         color_background='#4E2271')

    # 5层MoS2，每层对能带的贡献最多为1.0/5=0.2
    #vw.VisualElec_vaspkit.VisualizeProjectedBands(x_band,y_band,w_band,Knodes_projected=Kpath_nodes,
                                                  #colormap=cmap,size_band=np.abs(w_band)*4,
                                                  #y_major_tick=2,colormap_norm=(-1,1),HighSymPath=HighSymPath)

    #vw.SavingFigure(saving_directory=saving_directory, file_name=saving_filename)
    #vw.SavingFigure(saving_directory=saving_directory, file_name=saving_filename, format='eps')

    # 可视化并保存scalebar
    #vw.CustomizingColormap.ShowColorbar(cmap,(-1,1))
    #vw.SavingFigure(saving_directory=saving_directory, file_name=saving_filename+'_scalebar')
    #vw.SavingFigure(saving_directory=saving_directory, file_name=saving_filename+'_scalebar', format='eps')