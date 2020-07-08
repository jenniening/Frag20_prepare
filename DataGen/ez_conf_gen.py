import shutil

import DataGen
import pandas as pd
import logging
import os
import os.path as osp
import torch
import glob

from DataGen.genconfs import runGenerator
from DataGen.genGaus import file_transfer, file_seperate, gaussian_gen
from DataGen.genjob_hpc import gen_job, check_gaussian
from DataGen.analysis import check
from DataGen.prepare_data import prepare_torch, prepare_PhysNet_input, prepare_target_csv
from rdkit import Chem


def gen_header_text(n_cpu, mem, method='B3LYP/6-31G*', path='header.txt'):
    header = "%NProcShared={}\n%Mem={}GB\n\n#p {} opt freq\n\ntest\n".format(n_cpu, mem, method)
    with open(path, 'w') as f:
        f.write(header)
    return


def logger_setup(directory, logger_name):
    logging.basicConfig(filename=osp.join(directory, logger_name), format='%(asctime)s %(message)s', filemode='w')
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    return logger


def gen_conf(smiles_csv_file, source_data_name, out_dir, record_dir, num_conf=300, debug_mode=False, **kwargs):
    """
    This function generate conformations based on SMILES with RDKit
    :param debug_mode: if debug mode, we will only chose 10 molecules for generation
    :param record_dir: a directory for log files
    :param smiles_csv_file: csv file which contains two columns: 'index'(molecule index) and 'smiles'(chemical SMILES)
    :param source_data_name: ccdc | zinc | pubchem
    :param out_dir: output directory
    :param num_conf: number of conformations to be generated
    :param kwargs: other key word arguments used by runGenerator
    :return:
    """
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    if not osp.exists(record_dir):
        os.mkdir(record_dir)

    logger = logger_setup(record_dir, 'genConf.log')

    smiles_file = pd.read_csv(smiles_csv_file)
    index_list = smiles_file['index'].values.tolist()
    smiles_list = smiles_file['smiles'].values.tolist()
    if debug_mode:
        n_mol_chosen = 10
        logger.info('debug enabled, only {} molecules chosen')
        index_list = index_list[:n_mol_chosen]
        smiles_list = smiles_list[:n_mol_chosen]
    logger.info('index and smiles file loaded, generating conformations...')
    logger.info('total number of molecules queried: {}'.format(len(index_list)))
    logger.info('target number of conformations per molecule: {}'.format(num_conf))

    failed_list = runGenerator(index_list, smiles_list, source_data_name, datadir=out_dir, numConfs=num_conf, **kwargs)
    logger.info('{} of {} molecule(s) failed.'.format(len(failed_list), len(index_list)))
    torch.save(failed_list, osp.join(record_dir, 'genConf_failed_ids.pt'))


def gen_gaussian_input(input_dir, out_dir, record_dir, initial_index=0, n_cpu=1, mem=2, method='B3LYP/6-31G*'):
    """
    This function generate Gaussian input
    :param input_dir: directory containing generated conformation
    :param out_dir: directory for generated Gaussian input
    :param record_dir: directory for log files
    :param initial_index: generated input will be named $(n).opt.com starting from $(initial_index+1).opt.log
    :param n_cpu:
    :param mem:
    :param method:
    :return:
    """
    if not osp.exists(out_dir):
        os.mkdir(out_dir)
    logger = logger_setup(record_dir, 'genGaussian.log')

    logger.info('separating conformations...')
    infile_list = glob.glob(osp.join(input_dir, '*confors.sdf'))
    for i, path in enumerate(infile_list):
        infile_list[i] = path.split('/')[-1]
    infile_list.sort(key=lambda string: int(string.split('_')[0]))

    logger.info('{} molecules to be separated'.format(len(infile_list)))
    final_index = file_seperate(infile_list, input_dir, out_dir, initial_index=initial_index)
    logger.info('separation success, {} conformation files generated.'.format(final_index-initial_index-1))
    conf_index = [i for i in range(initial_index + 1, final_index + 1)]

    logger.info('converting *.sdf files into *.com files')
    file_transfer(conf_index, 'conformations', out_dir)
    logger.info('converting complete.')

    logger.info('generating Gaussian input files...')
    header_file = osp.join(record_dir, 'header.txt')
    gen_header_text(n_cpu, mem, method, header_file)
    gaussian_gen(conf_index, header_file, out_dir)
    logger.info('Gaussian files generation complete.')

    logger.info('removing temp files...')
    tmp_files = glob.glob(osp.join(out_dir, 'tmp*'))
    for tmp_file in tmp_files:
        os.remove(tmp_file)
    logger.info('success.')


def gen_slurm_jobs(job_dir, target_dir, num_jobs, time, mem, n_cpu):
    """
    Generate *.pbs jobs which can be directly submitted through slurm
    :param job_dir: directory storing all the jobs
    :param target_dir: directory for Gaussian input
    :param num_jobs: number of Gaussian calculations per *.pbs file
    :param time: number of hours per *.pbs job
    :param mem: memory required (GB)
    :param n_cpu: number of cpus per job
    :return:
    """
    if not osp.exists(job_dir):
        os.mkdir(job_dir)

    index_list = glob.glob(osp.join(target_dir, '*.opt.com'))
    for i, file in enumerate(index_list):
        index_list[i] = int(file.split('/')[-1].split('.')[0])
    index_list.sort()

    gen_job(index_list, 0, job_dir, target_dir, num_jobs, time, mem, n_cpu)


def prepare_dataset(data_dir, output_dir, record_dir, reference, file_name='data'):
    """
    remove Gaussian errors, performing post analysis and prepare dataset(PhysNet, target and mol files)
    :param data_dir: directory containing all Gaussian input and output
    :param output_dir:
    :param record_dir:
    :param reference: atomic reference energies for QM method, which is used to get the atomization energies
    :param file_name: name which will be part of prepared files
    :return:
    """
    if not osp.exists(output_dir):
        os.mkdir(output_dir)

    conf_indexes = glob.glob(osp.join(data_dir, '*.opt.log'))
    for i, file in enumerate(conf_indexes):
        conf_indexes[i] = int(file.split('/')[-1].split('.')[0])

    check_gaussian(conf_indexes, data_dir)
    shutil.move(osp.join(data_dir, 'Gaussian_error.csv'),
                osp.join(record_dir, '{}_Gaussian_error.csv'.format(file_name)))
    gaussian_error_list = []
    with open(osp.join(record_dir, '{}_Gaussian_error.csv'.format(file_name))) as f:
        for line in f.readlines():
            gaussian_error_list.append(int(line.strip()))

    '''
    Removing error files
    '''
    removed_dir = '../../removed'
    if not osp.exists(removed_dir):
        os.mkdir(removed_dir)
    for error_id in gaussian_error_list:
        shutil.move(osp.join(data_dir, '{}.sdf'.format(error_id)), osp.join(removed_dir, '{}.sdf'.format(error_id)))
        shutil.move(osp.join(data_dir, '{}.com'.format(error_id)), osp.join(removed_dir, '{}.com'.format(error_id)))
        shutil.move(osp.join(data_dir, '{}.opt.com'.format(error_id)),
                    osp.join(removed_dir, '{}.opt.com'.format(error_id)))
        shutil.move(osp.join(data_dir, '{}.opt.log'.format(error_id)),
                    osp.join(removed_dir, '{}.opt.log'.format(error_id)))

    '''
    Post analysis
    '''
    conf_indexes = glob.glob(osp.join(data_dir, '*.opt.log'))
    for i, file in enumerate(conf_indexes):
        conf_indexes[i] = int(file.split('/')[-1].split('.')[0])
    if not os.path.exists(osp.join(output_dir, 'data_consistent_strict_rmrpc.csv')):
        ana_check = check(data_dir, conf_indexes, output_dir, confs=True)
        ana_check.built_initialdata()
        ana_check.check_consistency('strict')
        ana_check.check_others('data_consistent_strict.csv')

    '''
    Calculate max number of atoms for PhysNet dataset preparation
    '''
    checked_data = pd.read_csv(osp.join(output_dir, 'data_consistent_strict_rmrpc.csv'))
    index_list = checked_data['index'].values.tolist()
    smiles_list = checked_data['SMILES']
    largest_num_atoms = 0
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        num_atom = mol.GetNumAtoms()
        if num_atom > largest_num_atoms:
            largest_num_atoms = num_atom

    '''
    Preparing QM and MMFF geometry dataset
    '''
    index_list = [str(i) for i in index_list]
    for method in ['QM', 'MMFF']:
        prepare_torch(index_list, osp.join(output_dir, '{}_{}.pt'.format(file_name, method)), data_dir, reference,
                      ftype='confs', method=method)
        prepare_PhysNet_input(index_list, osp.join(output_dir, '{}_{}_PhysNet'.format(file_name, method))
                              , data_dir, reference, method=method, ftype='conf', largest_num_atoms=largest_num_atoms)
    prepare_target_csv(index_list, osp.join(output_dir, 'target.csv'), data_dir, reference, 'conf')

    '''
    change names
    '''
    shutil.move(osp.join(output_dir, 'data_consistent_strict.csv'),
                osp.join(output_dir, '{}_consistent_strict.csv'.format(file_name)))
    shutil.move(osp.join(output_dir, 'data_consistent_strict_rmrpc.csv'),
                osp.join(output_dir, '{}_consistent_strict_rmrpc.csv'.format(file_name)))
    shutil.move(osp.join(output_dir, 'initial_dataset.csv'),
                osp.join(output_dir, '{}_initial_dataset.csv'.format(file_name)))
    shutil.move(osp.join(output_dir, 'target.csv'),
                osp.join(output_dir, '{}_target.csv'.format(file_name)))
    return


if __name__ == '__main__':
    # gen_conf('input/Index_ccdc_20_overlapwithmol20_removeempty.csv', 'ccdc', 'conformations', 'record',
    #          num_conf=10, debug_mode=True)
    # gen_gaussian_input('conformations', 'ginput', 'record')
    # gen_slurm_jobs('gjobs', 'ginput', 2, 5, 10, 1)
    prepare_dataset('ginput', 'output', 'record', 'input/atomref.B3LYP_631Gd.10As.npz', file_name='conf20_batch1')
