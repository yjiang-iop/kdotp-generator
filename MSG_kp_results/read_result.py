# Author: Yi Jiang, <jiangyi14@mails.ucas.ac.cn>, Institute of Physics, Chinese Academy of Sciences
"""
A script used to load the kp results of MSGs. Users can type in desired MSG number, coirrep label, and variable, 
and this script can print corresponding kp results.
"""
import numpy as np 
import os

def load_data(msg_num, variable):
    # variable in list: ['k','E','B','kE','kB', 'EB', 'epsilon']
    filename = 'result_' + variable 
    with open(msg_num + '/' + filename, 'r') as f:
        f = f.readlines()
       #for cnt, line in enumerate(f):
        assert 'Start calculation of msg ' in f[4] and 'Finish calculation of msg ' in f[-1]
        assert f[4].split()[-2] == msg_num
        cnt = 6
        msg_kpdata = []
        for ik in range(8): # at most 8 kpoints
            if '<<<<<<  Finish calculation of msg %s  >>>>>'%msg_num in f[cnt]:
                break
            assert 'Start kp calculation for msg %s  kpoint '%msg_num in f[cnt], (cnt,f[cnt]) 
            klabel = f[cnt].split()[-2]
            
            # read linear coirreps
            cnt += 3
            assert '+++  Irreducible linear co-representations ' in f[cnt]
            linear_coirs = f[cnt]
            cnt += 1
            while '+++++++++++++++++' not in f[cnt]:
                linear_coirs += f[cnt]
                cnt += 1
            cnt += 4

            # read kp results of each coirrep
            kdata_collect = []
            for ith_coir in range(30): # at most 30 coirreps

                assert '***** Start of MSG ' in f[cnt], cnt
                coir_label = f[cnt].split()[-2]
                # read coir data
                coir_rep_data = ''
                while 'Print kp results' not in f[cnt]:
                    coir_rep_data += f[cnt]
                    cnt += 1
                
                cnt += 3
                coir_kp_data = ''
                while '***** End of msg %s  %s  ****'%(msg_num, coir_label) not in f[cnt]:
                    coir_kp_data += f[cnt]
                    cnt += 1

                coir_dict = {'msg_num':msg_num, 'klabel':klabel, 'coir_label':coir_label, 'rep_data':coir_rep_data, 'kp_data':coir_kp_data}
                kdata_collect.append(coir_dict)
                cnt += 4

                if '===== Finish kp calculation for msg %s  kpoint %s ===='%(msg_num, klabel) in f[cnt-2]: 
                    break

            msg_kpdata.append(kdata_collect)

    return msg_kpdata


def read_coir_kp_data(data):
    # read the data for each coir
    data = data.split('\n')
    order_collect = []
    for cnt, line in enumerate(data):
        if '====  Result of msg ' in line:
            order_collect.append([line])
        else:
            order_collect[-1].append(line)

    dict_data = []
    for ith, ith_data in enumerate(order_collect):
        if 'No symmetry-allowed kp models.' in ith_data[1]:
            dict_data.append({'kpmodel':[], 'basis_vec':[]})
            continue
            
        cnt = 1
        while 'Warning' in ith_data[cnt]:
            cnt += 1
        assert 'Number of independent kp models: ' in ith_data[cnt], ith_data
        num_kp = int(ith_data[cnt].strip().split()[-1])

        ith_dict = {'kpmodel':[], 'basis_vec':[]}
        for cnt, line in enumerate(ith_data):
            if '-th kp model:' in line:
                ith_dict['kpmodel'].append(ith_data[cnt+1])

                basis_vec = ith_data[cnt+2][14:]
                ith_dict['basis_vec'].append(basis_vec)

        assert num_kp == len(ith_dict['kpmodel']), ('num not match!', num_kp, ith_dict['kpmodel'])
        dict_data.append(ith_dict)
        
    return dict_data






if __name__ == '__main__':
    msg_num = '222.103'
    coir = 'GM6d'
    variable = 'k'  # load the result for this variable, can be anyone of ['k', 'E', 'B', 'epsilon', 'kE', 'kB', 'EB']

    msg_data = load_data(msg_num, variable)

    coir_data = [ d for kdict in msg_data for d in kdict if d['coir_label'] == coir ][0]
    print('\nmsg_num: %s  coir_label: %s\n'%(coir_data['msg_num'], coir_data['coir_label']))
    print(coir_data['rep_data'])
    print(coir_data['kp_data'])

    #dict_data = read_coir_kp_data(coir_data['kp_data'])
