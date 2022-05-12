# For using this module, you have to create a virtual environment.
# Once created, you run the following command:
# pip install git+https://github.com/zhanglab/psamm-import.git
# 
# Or you can follow this tutorial (useful for further information):
# https://psamm.readthedocs.io/en/latest/index.html
#
# If attribute error due to fractions has no attribute "gcd", then modify the file
# fractions.py importing gcd as follows:
# from math import gcd
#
# Author: Sergio Pachón Dotor 
# Github: https://github.com/SergioPachonDotor

import os

model_name = 'iAF1260_mutant_delta_f6p_c'
model_path_sbml = f'./genome/Mutants/{model_name}.xml'
model_path_yaml= f'./genome/{model_name}_yaml/'
results_path= f'../../results/{model_name}/'

def import_existing_model(ex_model_path=model_path_sbml, destination_path=model_path_yaml):
    os.system(f'psamm-import sbml --source {ex_model_path} --dest {destination_path}')


def gapfind_and_gapcheck(solver='gurobi', model_path=model_path_yaml, deep_check=False, results_path=results_path, file_name='gapcheck', file_extension='txt'):
    
    os.chdir(model_path)

    if deep_check == True:
        os.system(f'psamm-model gapcheck --method prodcheck --solver name={solver} > {results_path}{model_name}_{file_name}_prodcheck.{file_extension}', objective= 'thr__L_c')
        os.system(f'psamm-model gapcheck --method sinkcheck --unrestricted-exchange --solver name={solver} > {results_path}{model_name}_{file_name}_sinkcheck_unrestricted-exchange.{file_extension}')
        os.system(f'psamm-model gapcheck --method gapfind --unrestricted-exchange --solver name={solver} > {results_path}{model_name}_{file_name}_gapfind_unrestricted-exchange.{file_extension}')

    if deep_check == False:
        os.system(f'psamm-model gapcheck --method gapfind --unrestricted-exchange --solver name={solver} > {results_path}{model_name}_{file_name}_gapfind_unrestricted-exchange.{file_extension}')


def gapfill_psamm(model_path=model_path_yaml, solver='gurobi', results_path=results_path, file_name='gapcheck', file_extension='txt'):
    os.chdir(model_path)    
    os.system(f'psamm-model gapfill --solver name={solver} --epsilon 0.0001 --compound thr__L_c > {results_path}{model_name}_{file_name}_gapfill.{file_extension}')

def export_model():
    pass


if __name__ == '__main__':
    # import_existing_model()
    # gapfind_and_gapcheck()
    # gapfill_psamm()
    pass
