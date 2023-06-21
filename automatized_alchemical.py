"""
This script is designed to run alchemical through the command-line.
"""

import os
import argparse as ap
import numpy as np

import peleffy
from peleffy.topology import Molecule

from peleffy.forcefield import OpenForceField
from peleffy.topology import Topology

from peleffy.topology import Alchemizer
from peleffy.template import Impact

DEFAULT_OFF_FORCEFIELD = 'openff_unconstrained-2.0.0.offxml'

def parse_args(args):
    """
    It parses the command-line arguments.
    Parameters
    ----------
    args : list[str]
        List of command-line arguments to parse
    Returns
    -------
    parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """

    parser = ap.ArgumentParser(prog="PELE Force Field Yielder for alchemistry (alchemistry peleffy)")
    parser.add_argument("pdb_file1", metavar="PDB FILE", type=str,
                        help="Path PDB file to parameterize",
                        default=None)
    parser.add_argument("pdb_file2", metavar="PDB FILE", type=str,
                        help="Path PDB file to parameterize",
                        default=None)
    parser.add_argument("lambdas", metavar="NATURAL NUMBER",type=int,
                        help="Number of lambdas (number of intermediate" +
                        "steps)", default=9)
    parser.add_argument('-b', '--balance',dest="balance",
                        action='store_true', help="Activate balance method. Defaut is linear FEP")
    parser.add_argument('-l', '--leo',dest="leo",
                        action='store_true', help="Activate leo method. Defaut is linear FEP")
    parser.add_argument("-f", "--forcefield", metavar="NAME",
                        type=str, help="OpenForceField's forcefield name. " +
                        "Default is " + str(DEFAULT_OFF_FORCEFIELD),
                        default=DEFAULT_OFF_FORCEFIELD),
    

    parser.set_defaults(balance=False)
    parser.set_defaults(leo=False)

    parsed_args = parser.parse_args(args)

    # Check force field
    from peleffy.forcefield import ForceFieldSelector

    selector = ForceFieldSelector()
    available_ffs = selector.get_list()

    if parsed_args.forcefield.lower() not in available_ffs:
        raise ValueError('Force field ' +
                         '\'{}\' '.format(parsed_args.forcefield) +
                         'is unknown')

    return parsed_args



def run_alchemical(pdb_file1, 
            pdb_file2,
            lambdas,
            balance,
            leo,
            forcefield_name = DEFAULT_OFF_FORCEFIELD,
            ):

    """
    It runs alchemistry.
    
    Parameters
    ----------
    pdb_file1 : str
        The path to the pdb_file to the reference ligand
    pdb_file2 : str
        The path to the pdb_file to final ligand
    lambdas : int
        The value to define the number of different lambdas. 
        This lambda affects all the parameters. Default is 9.
    balance : boolean
        Define if we want to use balance transformation.
    leo : boolean
        Define if we want to use leo transformation.
    forcefield_name : str
        The name of an OpenForceField's forcefield
    
    """

    molecule1 = Molecule(pdb_file1)
    molecule2 = Molecule(pdb_file2)
    openff = OpenForceField(forcefield_name)

    molecule1_params = openff.parameterize(molecule1)
    molecule2_params = openff.parameterize(molecule2)

    molecule1_top = Topology(molecule1, molecule1_params)
    molecule2_top = Topology(molecule2, molecule2_params)

    alchemistry = Alchemizer(molecule1_top, molecule2_top)

    alchemistry.molecule1_to_pdb('state1.pdb')
    alchemistry.molecule2_to_pdb('state2.pdb')
    alchemistry.hybrid_to_pdb('hybrid.pdb')


    # If no method is selected, the default is linear FEP
    if balance == False and leo == False:
        fep_lambdas = np.linspace(0.0, 1.0, lambdas+2)

        print("fep_lambdas = " + str(fep_lambdas))

        for fep in range(len(fep_lambdas)):
            alchemical_top = alchemistry.get_alchemical_topology(fep_lambda=round(fep_lambdas[fep],2))
    
            template = Impact(alchemical_top)
            template.to_file(f'hybz_{fep}')

            alchemistry.rotamer_library_to_file(
                        f'HYB.rot.assign{fep}',
                        fep_lambda=round(fep_lambdas[fep],2))

            alchemistry.obc_parameters_to_file(
                        f'ligandParams_{fep}.txt',
                        fep_lambda=round(fep_lambdas[fep],2))

       
    # If balance method is selected. Returns a balanced transformation beginning with coul1, then fep and finally coul2
    elif balance and leo == False:
        coul1_lambdas = [0.0] * (lambdas+2)
        fep_lambdas = [0.0] * (lambdas+2)
        coul2_lambdas = [0.0] * (lambdas+2)

        if (lambdas+2) % 3 == 1:
            segment = ((lambdas+2)//3) + 1 
            val = np.linspace(0, 1, segment)

            for i in range(0,segment):
                coul1_lambdas[i] = round(val[i],2)
                fep_lambdas[i+segment-1] = round(val[i],2)
                coul2_lambdas[i+segment+segment-2] = round(val[i],2)
            
            for i in range(segment,len(coul1_lambdas)):
                coul1_lambdas[i] = 1.0
            
            for i in range(2*segment-1, len(fep_lambdas)):
                fep_lambdas[i] = 1.0
                

        else:
            segment = ((lambdas+2)//3) + 1
            segment2 = (lambdas+2)-(2*segment-2)
            val1 = np.linspace(0, 1, segment)
            val2 = np.linspace(0, 1, segment2)

            for i in range(0,segment):
                coul1_lambdas[i] = round(val1[i],2)
                coul2_lambdas[i+segment+segment2-2] = round(val1[i],2)
            
            for i in range(0,segment2):
                fep_lambdas[i+segment-1] = round(val2[i],2)
            
            for i in range(segment,len(coul1_lambdas)):
                coul1_lambdas[i] = 1.0
            
            for i in range(segment+segment2-1, len(fep_lambdas)):
                fep_lambdas[i] = 1.0

        print("coul1_lambdas = " + str(coul1_lambdas))
        print("fep_lambdas = " + str(fep_lambdas))
        print("coul2_lambdas = " + str(coul2_lambdas))


        for i, (fep_lambda, coul1_lambda, coul2_lambda) \
                in enumerate(zip(fep_lambdas, coul1_lambdas, coul2_lambdas)):
            alchemical_top = alchemistry.get_alchemical_topology(fep_lambda=fep_lambda,
                                                                    coul1_lambda= coul1_lambda,
                                                                    coul2_lambda=coul2_lambda)
            
            template = Impact(alchemical_top)
            template.to_file(f'hybz_{i}')

            alchemistry.rotamer_library_to_file(
                                            f'HYB.rot.assign{i}',
                                            fep_lambda=fep_lambda,
                                            coul1_lambda=coul1_lambda,
                                            coul2_lambda=coul2_lambda)

            alchemistry.obc_parameters_to_file(
                                                f'ligandParams_{i}.txt',
                                                fep_lambda=fep_lambda,
                                                coul1_lambda=coul1_lambda,
                                                coul2_lambda=coul2_lambda)
            
    # If leo method is selected. Returns a transformation with the half of the lambdas will be on leonnard jones parameters.
    elif leo and balance == False:
        coul1_lambdas = [0.0] * (lambdas+2)
        fep_lambdas = [0.0] * (lambdas+2)
        coul2_lambdas = [0.0] * (lambdas+2)

        coul_segment = round((lambdas+2)*0.25) 
        leo_segment = lambdas+2 - 2*coul_segment + 2

        val1 = np.linspace(0, 1, coul_segment)
        val2 = np.linspace(0, 1, leo_segment)

        for i in range(0,coul_segment):
            coul1_lambdas[i] = round(val1[i],2)
            coul2_lambdas[i+coul_segment+leo_segment-2] = round(val1[i],2)
                    
        for i in range(0,leo_segment):
            fep_lambdas[i+coul_segment-1] = round(val2[i],2)
                    
        for i in range(coul_segment,len(coul1_lambdas)):
            coul1_lambdas[i] = 1.0
                    
        for i in range(coul_segment+leo_segment-1, len(fep_lambdas)):
            fep_lambdas[i] = 1.0

        print("coul1_lambdas = " + str(coul1_lambdas))
        print("fep_lambdas = " + str(fep_lambdas))
        print("coul2_lambdas = " + str(coul2_lambdas))

        for i, (fep_lambda, coul1_lambda, coul2_lambda) \
                in enumerate(zip(fep_lambdas, coul1_lambdas, coul2_lambdas)):
            alchemical_top = alchemistry.get_alchemical_topology(fep_lambda=fep_lambda,
                                                                    coul1_lambda= coul1_lambda,
                                                                    coul2_lambda=coul2_lambda)
            
            template = Impact(alchemical_top)
            template.to_file(f'hybz_{i}')

            alchemistry.rotamer_library_to_file(
                                            f'HYB.rot.assign{i}',
                                            fep_lambda=fep_lambda,
                                            coul1_lambda=coul1_lambda,
                                            coul2_lambda=coul2_lambda)

            alchemistry.obc_parameters_to_file(
                                                f'ligandParams_{i}.txt',
                                                fep_lambda=fep_lambda,
                                                coul1_lambda=coul1_lambda,
                                                coul2_lambda=coul2_lambda)
            
    else:
        print("Select just one option: balance or leo")


def main(args):
    """
    It reads the command-line arguments and runs peleffy.

    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    
    Examples
    --------
    
    From the command-line:
    
    >>> python script.py molecule1.pdb molecule2.pdb 14 
    -f openff_unconstrained-1.2.0.offxml -c
    
    """
    run_alchemical(pdb_file1=args.pdb_file1,
                pdb_file2=args.pdb_file2,
                lambdas=args.lambdas,
                balance=args.balance,
                leo=args.leo,
                forcefield_name=args.forcefield,
                )


if __name__ == '__main__':
    import sys
    args = parse_args(sys.argv[1:])
    main(args)

