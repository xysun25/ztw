import numpy as np
import pandas as pd
import scipy 
import re
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

def Init_data()
    element_dict = {'atoms': 0
               ,'bonds': 0
               ,'angles': 0
               ,'torsions': 0
               ,'improper': 0 
               ,'dihedrals': 0}
    

    styl_dict = {'atoms_style': 0
                ,'bonds_style': 0
                ,'angles_style': 0
                ,'torsions_style': 0
                ,'improper_style': 0
                ,'dihedrals_style':0}
                
    return element_dict, styl_dict
    

def Input_Param(input_file):
    input_value = pd.read_table(input_file, header = None)
    input_param = {'input_file_counter': int(input_value[0][0])
                  ,'file_name': input_value[0][1].split()
                  ,'file_rounds': input_value[0][2].split()
                  ,'xyz_file': input_value[0][3]
                  ,'X_limit': input_value[0][4].split()
                  ,'Y_limit': input_value[0][5].split()
                  ,'Z_limit': input_value[0][6].split()
                  ,'Condition': input_value[0][7].split()
                  ,'Bond_Distance_and_Angle': input_value[0][8].split()}
    return input_param
    
def Input_Param_save_molecula_range(input_file):
    input_value = pd.read_table(input_file, header = None)
    input_param = {'input_file_counter': int(input_value[0][0])
                  ,'file_name': input_value[0][1].split()
                  ,'file_rounds': input_value[0][2].split()
                  ,'xyz_file': input_value[0][3]
                  ,'X_limit': input_value[0][4].split()
                  ,'Y_limit': input_value[0][5].split()
                  ,'Z_limit': input_value[0][6].split()}
    return input_param 
    
def Counter(element_dict, input_param, input_data):
    element_dict_one_list = []
    for j in range(input_param['input_file_counter']):
        for i in range(input_data[j].index.stop):
            for k in range(int(input_param['file_rounds'][j])):
                value = input_data[j][0][i]
                if 'of atoms' in value:
                    element_dict['atoms'] += int(input_data[j][0][i+1].split()[0])
                elif 'of bonds' in value:
                    element_dict['bonds'] += int(input_data[j][0][i+1].split()[0])
                elif 'of angles' in value:
                    element_dict['angles'] += int(input_data[j][0][i+1].split()[0])
                elif 'of torsions' in value:
                    element_dict['dihedrals'] += int(input_data[j][0][i+1].split()[0])
                elif 'of dihedrals' in value:
                    element_dict['dihedrals'] += int(input_data[j][0][i+1].split()[0])
                elif '  improper' in value:
                    element_dict['improper'] = element_dict['improper'] + 1
    for j in range(input_param['input_file_counter']):
        element_dict_one = {'atoms': 0
                           ,'bonds': 0
                           ,'angles': 0
                           ,'torsions': 0
                           ,'improper': 0 
                           ,'dihedrals': 0}
        for i in range(input_data[j].index.stop):
            value = input_data[j][0][i]
            if 'of atoms' in value:
                element_dict_one['atoms'] = int(input_data[j][0][i+1].split()[0])
            elif 'of bonds' in value:
                element_dict_one['bonds'] = int(input_data[j][0][i+1].split()[0])
            elif 'of angles' in value:
                element_dict_one['angles'] = int(input_data[j][0][i+1].split()[0])
            elif 'of torsions' in value:
                element_dict_one['dihedrals'] = int(input_data[j][0][i+1].split()[0])
            elif 'of dihedrals' in value:
                element_dict_one['dihedrals'] = int(input_data[j][0][i+1].split()[0])
            elif '  improper' in value:
                element_dict_one['improper'] = element_dict_one['improper'] + 1
        element_dict_one_list.append(element_dict_one)
    return element_dict_one_list
    
    
def Save_Atoms_data(element_dict_one_list, input_param, input_data):
    atoms_value_all = []
    for k in range(input_param['input_file_counter']):  
        atoms_value = []
        for i in range(input_data[k].index.stop):
            value = input_data[k][0][i]
            if element_dict_one_list[k]['atoms'] != 0:
                if 'of atoms' in value:
                    atoms_index = i
                    if 'epsilon' in input_data[k][0][atoms_index + 2]:
                        for j in range(element_dict_one_list[k]['atoms']):
                            atoms_value.append(input_data[k][0][atoms_index + 4 + j])
                    else:
                        for j in range(element_dict_one_list[k]['atoms']):
                            atoms_value.append(input_data[k][0][atoms_index + 2 + j])
        atoms_value_all.append(pd.DataFrame(atoms_value))
    return atoms_value_all 

def Save_Bonds_data(element_dict_one_list, input_param, input_data):
    bonds_value_all = []
    for k in range(input_param['input_file_counter']):  
        bonds_value = []
        for i in range(input_data[k].index.stop):
            value = input_data[k][0][i]
            if element_dict_one_list[k]['bonds'] != 0:
                if 'of bonds' in value:
                    bonds_index = i
                    if 'kr(kcal' in input_data[k][0][bonds_index + 2]:
                        for j in range(element_dict_one_list[k]['bonds']):
                            bonds_value.append(input_data[k][0][bonds_index + 3 + j])
                    else:
                        for j in range(element_dict_one_list[k]['bonds']):
                            bonds_value.append(input_data[k][0][bonds_index + 2 + j])
        bonds_value_all.append(pd.DataFrame(bonds_value))
    return bonds_value_all
    
def Save_Angles_data(element_dict_one_list, input_param, input_data):
    angles_value_all = []
    for k in range(input_param['input_file_counter']):  
        angles_value = []
        for i in range(input_data[k].index.stop):
            value = input_data[k][0][i]
            if element_dict_one_list[k]['angles'] != 0:
                if 'of angles' in value:
                    angles_index = i
                    if 'c1' in input_data[k][0][angles_index + 2]:
                        for j in range(element_dict_one_list[k]['angles']):
                            angles_value.append(input_data[k][0][angles_index + 3 + j])
                    else:
                        for j in range(element_dict_one_list[k]['angles']):
                            angles_value.append(input_data[k][0][angles_index + 2 + j])
        angles_value_all.append(pd.DataFrame(angles_value))
    return angles_value_all
    
    
def Save_Dihedrals_data(element_dict_one_list, input_param, input_data):
    dihedrals_value_all = []
    for k in range(input_param['input_file_counter']):  
        dihedrals_value = []
        for i in range(input_data[k].index.stop):
            value = input_data[k][0][i]
            if element_dict_one_list[k]['dihedrals'] != 0:
                if 'of dihedrals' in value:
                    dihedrals_index = i
                    for j in range(element_dict_one_list[k]['dihedrals']):
                        dihedrals_value.append(input_data[k][0][dihedrals_index + 2 + j])
                elif 'of torsions' in value:
                    dihedrals_index = i
                    for j in range(element_dict_one_list[k]['dihedrals']):
                        dihedrals_value.append(input_data[k][0][dihedrals_index + 2 + j])
        dihedrals_value_all.append(pd.DataFrame(dihedrals_value))
    return dihedrals_value_all


def Read_file(file_name, input_param):
    file = open(file_name, 'r')
    input_param = Input_Param('input.in')
    file_data = pd.Series(file.readlines())
    data_x = []
    data_y = []
    data_z = []
    data_name = []
    for i in file_data:
        if re.findall('\w+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*', i):
            data_name.append(re.findall('\w+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*', i.replace('\n', ''))[0].split()[0])
            data_x.append(re.findall('\w+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*', i.replace('\n', ''))[0].split()[1])
            data_y.append(re.findall('\w+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*', i.replace('\n', ''))[0].split()[2])
            data_z.append(re.findall('\w+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*\-?\d+\.\d+\s*', i.replace('\n', ''))[0].split()[3])
    data = pd.DataFrame({'name':data_name, 'X': data_x, 'Y': data_y, 'Z': data_z})
    return data
    
def Find_all_bond_composition(input_param, bond, postion_limit, data, element_dict_one_list):
    index = 0
    atom_round = 0
    value_bond = []
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            atom_round += 1
            name_index = 0
            for j in range(element_dict_one_list[i]['atoms']):
                if data['name'][index] == bond[0] :
                    name_index += 1
                    if name_index == int(bond[1]):
                        if float(data['X'][index]) >= float(postion_limit['X'][0]) and float(data['X'][index]) <= float(postion_limit['X'][1]) :
                            if float(data['Y'][index]) >= float(postion_limit['Y'][0]) and float(data['Y'][index]) <= float(postion_limit['Y'][1]) :
                                if float(data['Z'][index]) >= float(postion_limit['Z'][0]) and float(data['Z'][index]) <= float(postion_limit['Z'][1]) :
                                    value_bond.append({'Name': data['name'][index], 'molecule_positon': atom_round, 'atoms_positon': index + 1, 'X': data['X'][index], 'Y': data['Y'][index], 'Z': data['Z'][index], 'molecule_style': input_param['file_name'][i]})
                index += 1
    return value_bond
    
def Find_all_bond_composition_no_appoint(input_param, element, postion_limit, data, element_dict_one_list):
    index = 0
    atom_round = 0
    value_bond = []
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            atom_round += 1
            for j in range(element_dict_one_list[i]['atoms']):
                if data['name'][index] == element :
                    if float(data['X'][index]) >= float(postion_limit['X'][0]) and float(data['X'][index]) <= float(postion_limit['X'][1]) :
                        if float(data['Y'][index]) >= float(postion_limit['Y'][0]) and float(data['Y'][index]) <= float(postion_limit['Y'][1]) :
                            if float(data['Z'][index]) >= float(postion_limit['Z'][0]) and float(data['Z'][index]) <= float(postion_limit['Z'][1]) :
                                value_bond.append({'Name': data['name'][index], 'molecule_positon': atom_round, 'atoms_positon': index + 1, 'X': data['X'][index], 'Y': data['Y'][index], 'Z': data['Z'][index], 'molecule_style': input_param['file_name'][i]})
                index += 1
    return value_bond
    
    
def Input_data(input_param):
    input_data = []
    for i in range(input_param['input_file_counter']):
        file = open(input_param['file_name'][i])
        input_data.append(pd.DataFrame(file.readlines()))
    input_data.append(pd.read_table(input_param['xyz_file'], header = None))
    return input_data

def Bond_calculation(bond_composition_1, bond_composition_2, distance_limit):
    bond_suspected = []
    for i in bond_composition_1:
        for j in bond_composition_2:
            distance = ((float(i['X']) - float(j['X']))**2 + (float(i['Y']) - float(j['Y']))**2 + (float(i['Z']) - float(j['Z']))**2)**0.5
            if i['molecule_positon'] != j['molecule_positon']:
                if i['molecule_style'] == 'spc.mmol':
                    if j['molecule_style'] == 'spc.mmol':
                        bond_suspected.append({'bond_positon_molecule':[i['molecule_positon'], j['molecule_positon']]
                                            , 'bond_name' : i['Name'] + '-' + j['Name']
                                            , 'bond_distance' : distance
                                            , 'bond_positon_atoms':[i['atoms_positon'], j['atoms_positon']]})
    bond_true = [x for x in bond_suspected if x['bond_distance'] < distance_limit]
    return bond_suspected, bond_true
    

def output_data_v1(input_param, data, element_dict_one_list, file, distance, number):
    postion_limit = {'X': [float(i) for i in input_param['X_limit']], 'Y': [float(i) for i in input_param['Y_limit']], 'Z': [float(i) for i in input_param['Z_limit']]}
    bond = input_param['Condition']
    if len(input_param['Condition']) == 4:
        output = Bond_calculation(Find_all_bond_composition(input_param, bond[:2], postion_limit, data, element_dict_one_list), Find_all_bond_composition(input_param, bond[2:], postion_limit, data, element_dict_one_list), float(input_param['Bond_Distance_and_Angle'][0]))[1]
    elif len(input_param['Condition']) == 2:
        output = Bond_calculation(Find_all_bond_composition_no_appoint(input_param, bond[0], postion_limit, data, element_dict_one_list), Find_all_bond_composition_no_appoint(input_param, bond[1], postion_limit, data, element_dict_one_list), float(input_param['Bond_Distance_and_Angle'][0]))[1]
    Bond_distance_all = 0
    for i in output:
        Bond_distance_all += i['bond_distance']
        distance += i['bond_distance']
        number += 1
    print('共有{}个满足要求！\t \t 平均键长为{}Å \n'.format(len(output), Bond_distance_all / len(output)))
    file.write('共有{}个满足要求！\t \t 平均键长为{}Å \n'.format(len(output), Bond_distance_all / len(output)))
    return distance, number

    

def Read_file_one_frame_to_dataframe(one_frame):
    data_X = []
    data_Y = []
    data_Z = []
    data_name = []
    for i in one_frame[2:]:
        data_name.append(i.split()[0])
        data_X.append(i.split()[1])
        data_Y.append(i.split()[2])
        data_Z.append(i.split()[3])
    data = pd.DataFrame({'name':data_name, 'X': data_X, 'Y': data_Y, 'Z': data_Z})
    return data


def Read_file_one_frame(one_frame):
    for index, value in enumerate(file_data):
        if 'Atoms. Timestep:' in i:
            flag = index
            data_name = []
            data_X = []
            data_Y = []
            data_Z = []
            for k in range(Atoms_num):
                data_name.append(file_data[index + k].replace('\n', '')[0].split()[0])
                data_X.append(file_data[index + k].replace('\n', '')[0].split()[1])
                data_Y.append(file_data[index + k].replace('\n', '')[0].split()[2])
                data_Z.append(file_data[index + k].replace('\n', '')[0].split()[3])
            data = pd.DataFrame({'name':data_name, 'X': data_X, 'Y': data_Y, 'Z': data_Z})
            data_all.append(data)
    return data_all
    
def Read_file_all(input_param):
    file = open(file_name, 'r')
    file_data = pd.Series(file.readlines())
    Atoms_num = int(file_data[0].replace('\n', ''))
    data_all = []
    for index, value in enumerate(file_data):
        if 'Atoms. Timestep:' in i:
            flag = index
            data_name = []
            data_X = []
            data_Y = []
            data_Z = []
            for k in range(Atoms_num):
                data_name.append(file_data[index + k].replace('\n', '')[0].split()[0])
                data_X.append(file_data[index + k].replace('\n', '')[0].split()[1])
                data_Y.append(file_data[index + k].replace('\n', '')[0].split()[2])
                data_Z.append(file_data[index + k].replace('\n', '')[0].split()[3])
            data = pd.DataFrame({'name':data_name, 'X': data_X, 'Y': data_Y, 'Z': data_Z})
            data_all.append(data)
    return data_all
    
def Read_file_chosed(input_param):
    file = open(input_param['xyz_file'], 'r')
    file_data = pd.Series(file.readlines())
    Atoms_num = int(file_data[0].replace('\n', ''))
    one_Atomes_round = Atoms_num + 2
    frame_chosed_data = []
    frames = 0
    for i in file_data:
        if 'Atoms.' in i:
            frames += 1
    print('the {} has {:4}  frames'.format(input_param['xyz_file'], frames))
    chosed_list = input("please input the frames you want to choose: all, one or range such as all or 1 2 3, or 1-5: ")
    if chosed_list == 'all':
        chosed_list = list(range(1, frames + 1))
    elif '-' in chosed_list:
        chosed_list = list(range(int(chosed_list.split('-')[0]), int(chosed_list.split('-')[1]) + 1))
    else:
        chosed_list = list(chosed_list.split())
    for i in chosed_list:
        frame_chosed = {'index': int(i), 'data': file_data[(int(i) - 1) * one_Atomes_round: int(i) * one_Atomes_round]}
        frame_chosed_data.append(frame_chosed)
    return frame_chosed_data
    
def Read_file_all_frame_to_dataframe(all_frame):
    data = []
    for one_frame in all_frame:
        data.append({'frames': one_frame['index'], 'data': Read_file_one_frame_to_dataframe(one_frame['data'])})
    return data

def Read_file_one_frame_to_dataframe(one_frame):
    data_X = []
    data_Y = []
    data_Z = []
    data_name = []
    for i in one_frame.iloc[2:]:
        data_name.append(i.split()[0])
        data_X.append(i.split()[1])
        data_Y.append(i.split()[2])
        data_Z.append(i.split()[3])
    data = pd.DataFrame({'name':data_name, 'X': data_X, 'Y': data_Y, 'Z': data_Z, 'step': one_frame.iloc[1]})
    return data
    
def output_data_choosed_frames(input_param, data_all_dataframe, element_dict_one_list, out):
    distance = 0
    number = 0
    file = open(out, 'w', encoding='utf-8')
    for one_frame in data_all_dataframe:
        print('No.{}\t frame:'.format(one_frame['frames']))
        file.write('No.{}\t frame:'.format(one_frame['frames']))
        distance, number = output_data_v1(input_param, one_frame['data'], element_dict_one_list, file, distance, number)
    if distance != 0 and number != 0: 
        print('Total frames: {:8} \n Average number per frame:{:2.6f} \t Average distance per frame:{:2.6f}Å'.format(len(data_all_dataframe), number / len(data_all_dataframe), distance / number))
        file.write('Total frames: {:8} \n Average number per frame:{:2.6f} \t Average distance per frame:{:2.6f}Å'.format(len(data_all_dataframe), number / len(data_all_dataframe), distance / number))    
    file.close()
    
def Find_all_bond_composition(input_param, bond, postion_limit, data, element_dict_one_list):
    index = 0
    atom_round = 0
    value_bond = []
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            atom_round += 1
            name_index = 0
            for j in range(element_dict_one_list[i]['atoms']):
                if data['name'][index] == bond[0] :
                    name_index += 1
                    if name_index == int(bond[1]):
                        if float(data['X'][index]) >= float(postion_limit['X'][0]) and float(data['X'][index]) <= float(postion_limit['X'][1]) :
                            if float(data['Y'][index]) >= float(postion_limit['Y'][0]) and float(data['Y'][index]) <= float(postion_limit['Y'][1]) :
                                if float(data['Z'][index]) >= float(postion_limit['Z'][0]) and float(data['Z'][index]) <= float(postion_limit['Z'][1]) :
                                    value_bond.append({'Name': data['name'][index], 'molecule_positon': atom_round, 'atoms_positon': index + 1, 'X': data['X'][index], 'Y': data['Y'][index], 'Z': data['Z'][index], 'molecule_style': input_param['file_name'][i]})
                index += 1
    return value_bond      

def calculate_angle_1_one_frame(input_param, data, element_dict_one_list):
    index = 0
    angle_1_list = []
    postion_limit = {'X': [float(i) for i in input_param['X_limit']], 'Y': [float(i) for i in input_param['Y_limit']], 'Z': [float(i) for i in input_param['Z_limit']]}
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            if float(data['X'][index]) >= float(postion_limit['X'][0]) and float(data['X'][index]) <= float(postion_limit['X'][1]) :
                if float(data['Y'][index]) >= float(postion_limit['Y'][0]) and float(data['Y'][index]) <= float(postion_limit['Y'][1]) :
                    if float(data['Z'][index]) >= float(postion_limit['Z'][0]) and float(data['Z'][index]) <= float(postion_limit['Z'][1]) :
                        if input_param['file_name'][i] == 'spc.mmol':
                            ab = np.array([float(data['X'][index + 1]) - float(data['X'][index]), float(data['Y'][index + 1]) - float(data['Y'][index]), float(data['Z'][index + 1]) - float(data['Z'][index])])
                            ac = np.array([float(data['X'][index + 2]) - float(data['X'][index]), float(data['Y'][index + 2]) - float(data['Y'][index]), float(data['Z'][index + 2]) - float(data['Z'][index])])
                            ad = ab + ac
                            z_vector = np.array([0, 0, 1])
                            cos_ = np.dot(ad, z_vector) / (np.linalg.norm(ad) * np.linalg.norm(z_vector))
                            angle_1_list.append(np.degrees(np.arccos(cos_)))
            index += element_dict_one_list[i]['atoms']
    return angle_1_list
    

def calculate_angle_1_all_frame_test(input_param, data_all, element_dict_one_list):
    angle_1_all_list = []
    angle_1_all = []
    for i in data_all:
        angle_one_frame = calculate_angle_1_one_frame(input_param, i['data'], element_dict_one_list)
        angle_1_all_list.append(angle_one_frame)
        print('No.{}\t frame:'.format(i['frames']))
        print(angle_one_frame)
    for i in angle_1_all_list:
        for j in i:
            angle_1_all.append(j)
    sns.kdeplot(angle_1_all, bw_method = 'silverman', bw_adjust = 0.7)
    plt.xlim(0, 180)
    plt.ylabel('Probability', fontsize=20)
    plt.savefig("angle_1.png")
    return None
    
def calculate_angle_1_all_frame(input_param, data_all, element_dict_one_list, out):
    angle_1_all_list = []
    angle_1_all = []
    file = open(out, 'w', encoding='utf-8')
    for i in data_all:
        angle_one_frame = calculate_angle_1_one_frame(input_param, i['data'], element_dict_one_list)
        file.write('No.{}\t frame:\n'.format(i['frames']))
        for index, j in enumerate(angle_one_frame):
            file.write(str(round(float(j),2)) + ' ')
            if (index + 1) %10 == 0:
                file.write('\n')
        file.write('\n')
        angle_1_all_list.append(angle_one_frame)
        print('No.{}\t frame:'.format(i['frames']))
        print(angle_one_frame)
    for i in angle_1_all_list:
        for j in i:
            angle_1_all.append(j)
    file.write('-----------------------------------------------------------------------------\n')
    file.write('-----------------------------------------------------------------------------\n')
    file.write('-----------------------------------------------------------------------------\n')
    file.write('All_angles for {}  frames:\n'.format(len(data_all)))
    for index, i in enumerate(angle_1_all):
        file.write(str(round(float(i),2)) + ' ')
        if (index + 1)%10 == 0:
            file.write('\n')
    file.write('\n')
    file.close()
    sns.kdeplot(angle_1_all, bw_method = 'silverman', bw_adjust = 0.7)
    plt.xlim(0, 180)
    plt.ylabel('Probability', fontsize=20)
    plt.savefig("angle_1.png")
    return None  
    
def calculate_angle_2_one_frame(input_param, data, element_dict_one_list):
    index = 0
    angle_2_list = []
    postion_limit = {'X': [float(i) for i in input_param['X_limit']], 'Y': [float(i) for i in input_param['Y_limit']], 'Z': [float(i) for i in input_param['Z_limit']]}
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            if float(data['X'][index]) >= float(postion_limit['X'][0]) and float(data['X'][index]) <= float(postion_limit['X'][1]) :
                if float(data['Y'][index]) >= float(postion_limit['Y'][0]) and float(data['Y'][index]) <= float(postion_limit['Y'][1]) :
                    if float(data['Z'][index]) >= float(postion_limit['Z'][0]) and float(data['Z'][index]) <= float(postion_limit['Z'][1]) :
                        if input_param['file_name'][i] == 'spc.mmol':
                            ab = np.array([float(data['X'][index + 1]) - float(data['X'][index]), float(data['Y'][index + 1]) - float(data['Y'][index]), float(data['Z'][index + 1]) - float(data['Z'][index])])
                            ac = np.array([float(data['X'][index + 2]) - float(data['X'][index]), float(data['Y'][index + 2]) - float(data['Y'][index]), float(data['Z'][index + 2]) - float(data['Z'][index])])
                            ad = np.cross(ab, ac)
                            z_vector = np.array([0, 0, 1])
                            cos_ = np.dot(ad, z_vector) / (np.linalg.norm(ad) * np.linalg.norm(z_vector))
                            angle_2_list.append(np.degrees(np.arccos(cos_)))
            index += element_dict_one_list[i]['atoms']
    return angle_2_list
    

def calculate_angle_2_all_frame(input_param, data_all, element_dict_one_list, out):
    angle_2_all_list = []
    angle_2_all = []
    file = open(out, 'w', encoding='utf-8')
    for i in data_all:
        angle_one_frame = calculate_angle_2_one_frame(input_param, i['data'], element_dict_one_list)
        file.write('No.{}\t frame:\n'.format(i['frames']))
        for index, j in enumerate(angle_one_frame):
            file.write(str(round(float(j),2)) + ' ')
            if (index + 1) %10 == 0:
                file.write('\n')
        file.write('\n')
        angle_2_all_list.append(angle_one_frame)
        print('No.{}\t frame:'.format(i['frames']))
        print(angle_one_frame)
    for i in angle_2_all_list:
        for j in i:
            angle_2_all.append(j)
    file.write('-----------------------------------------------------------------------------\n')
    file.write('-----------------------------------------------------------------------------\n')
    file.write('-----------------------------------------------------------------------------\n')
    file.write('All_angles for {}  frames:\n'.format(len(data_all)))
    sns.kdeplot(angle_2_all, bw_method = 'silverman', bw_adjust = 0.7)
    for index, i in enumerate(angle_2_all):
        file.write(str(round(float(i),2)) + ' ')
        if (index + 1)%10 == 0:
            file.write('\n')
    file.write('\n')
    file.close()
    plt.xlim(0, 180)
    plt.ylabel('Probability', fontsize=20)
    plt.savefig("angle_2.png")
    return None


def Atoms_weight(element_dict_one_list, input_data, input_param):
    flag = 0
    molecula_all = []
    for i in range(input_param['input_file_counter']):
        molecula = []
        for index, value in enumerate(input_data[i][0]):
            if '#                     (A)' in value:
                flag = index
        for j in range(element_dict_one_list[i]['atoms']):
            molecula.append({'style': input_data[i][0][flag + 1 + j].split()[0], 'weight': float(input_data[i][0][flag + 1 + j].split()[4])})
        molecula_all.append(molecula)
    return molecula_all
 

def Save_molecules_in_range_one_frame(input_param, data, element_dict_one_list, atoms_weight_list):
    index = 0
    molecules = []
    postion_limit = {'X': [float(i) for i in input_param['X_limit']], 'Y': [float(i) for i in input_param['Y_limit']], 'Z': [float(i) for i in input_param['Z_limit']]}
    for i in range(input_param['input_file_counter']):
        for k in range(int(input_param['file_rounds'][i])):
            x_m = 0.0
            y_m = 0.0
            z_m = 0.0
            mass_sum = 0.0
            atom = []
            for j in range(element_dict_one_list[i]['atoms']):
                x_m = x_m + float(atoms_weight_list[i][j]['weight']) * float(data['X'][index])  
                y_m = y_m + float(atoms_weight_list[i][j]['weight']) * float(data['Y'][index])  
                z_m = z_m + float(atoms_weight_list[i][j]['weight']) * float(data['Z'][index])  
                mass_sum = mass_sum + float(atoms_weight_list[i][j]['weight'])
                atom.append({'name': data['name'][index], 'X': data['X'][index], 'Y': data['Y'][index], 'Z': data['Z'][index]})
                index += 1
            x_com = x_m / mass_sum
            y_com = y_m / mass_sum
            z_com = z_m / mass_sum
            if z_com >= postion_limit['Z'][0] and z_com <= postion_limit['Z'][1]:
                molecules.append({'style': input_param['file_name'][i].split('.')[0], 'number': int(k), 'data': atom})
    return molecules
    
    
def Save_molecules_in_range_all_frame(input_param, data_all, element_dict_one_list, out, atoms_weight_list):
    file = open(out, 'w', encoding='utf-8')
    for one_frame in data_all:
        data = Save_molecules_in_range_one_frame(input_param, one_frame['data'], element_dict_one_list, atoms_weight_list)
        #file.write('No.{:5} frame: There are {} molecules in all\n'.format(int(one_frame['frames']), len(data)))
        num = 0
        for number in data:
            num += len(number['data'])
        file.write('{}\n'.format(num))
        file.write(str(one_frame['data']['step'][0]))
        for molecula in data:
            # file.write('the moulecula is {} the number is {:5}:\n'.format(molecula['style'], int(molecula['number'])))
            for atom in molecula['data']:
                file.write('{} {:2.5f} {:2.5f} {:2.5f}\n'.format(atom['name'], float(atom['X']), float(atom['Y']), float(atom['Z'])))
            print(molecula)
        #file.write('\n\n')
    file.close()
    
    
