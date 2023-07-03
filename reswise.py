#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from itertools import product

"""Usage Ex.

    ./reswise.py -x equilibration.gro -y analysis.xtc -cut 0.6 -r1 POP2
"""

def parse_gro(gro):
    ''' Parse .gro file '''

    with open(gro, 'r') as fp:
        lines = fp.readlines()
        lines = [i.split(' ') for i in lines]        
    
    # Loop through and get unique residue names
    
    protein = ('MET',
                'ALA',
                'LEU',
                'ASP',
                'GLY',
                'ILE',
                'ARG',
                'PRO',
                'CYS',
                'TYR',
                'THR',
                'TRP',
                'GLU',
                'SER',
                'VAL',
                'HIS',
                'ASN',
                'LYS',
                'GLN',
                'PHE'
                )
    
    residues = []
    res_ids = []
    for n, i in enumerate(lines):
        if n > 1 and n < len(lines)-1:
            i = [j for j in i if j != '']
            res_str = ''.join([j for j in i[0] if j.isdigit() == False])
            res_num = ''.join([j for j in i[0] if j.isdigit() == True])
            if res_str in protein:
                res_ids.append(int(res_num))
                residues.append(res_str)
    
    # Check if protein residues are sequential first

    broken = False
    for n, i in enumerate(residues):
        if n > 1:
            if residues[n] != residues[n-1]:
                if res_ids[n] - res_ids[n-1] != 1:
                    print('Warning: protein residues are not sequential.')
                    broken = True
                    break
    
    if broken == False:
        print(f'Identified residues {min(res_ids)} to {max(res_ids)} as Protein')
        return res_ids
    else:
        return res_ids

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', help='.gro file')
    parser.add_argument('-y', nargs='+', help='.xtc file(s)')
    parser.add_argument('-r1', help='name of non-protein residue to count contacts for (e.g. "POP2")')
    parser.add_argument('-cut', type=float, default=0.6, help='threshold distance for counting contacts between atom pairs')
    
    args = parser.parse_args()
    res_ids = parse_gro(args.x)
    
    # Initialize dictionary of contact count 0
    contacts = dict()
    for res in res_ids:
        contacts[res] = 0
    
    # Loop through all protein residue ids and compute number of contacts using cutoff specified   
    for res in res_ids:
        for k in range(len(args.y)):
            iterator = md.iterload(args.y[k], top=args.x, chunk=1000)
            for n, traj in enumerate(iterator):
                sele = traj.topology.select('resid {}'.format(res))
                other = traj.topology.select('resname {}'.format(args.r1))
                pairs = np.array([i for i in product(sele, other)])
                d = md.compute_distances(traj, pairs)
                num_contacts = len(np.where(d < args.cut)[0])
                print(f"Residue ID {res} has {num_contacts} contacts at trajectory {k} chunk {n}")
                contacts[res] += num_contacts

    # Write data to file
    with open('contacts.dat', 'w') as fp:
        for key, value in contacts.items():
            fp.write('{} {}\n'.format(key, value))
