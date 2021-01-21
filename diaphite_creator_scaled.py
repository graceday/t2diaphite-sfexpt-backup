# -*- coding: utf-8 -*-
"""
Diaphite Creator

A program to create layers of diamond-graphite mixes, known as diaphite.
Adapted from a Pascal program of unknown provenance.

@author: Matt Bailey
@date: 2020-01-06
"""

import argparse
from typing import Iterable, Tuple

import numpy as np

from layer_utils import LayerTypes, read_sequence_from_file, sequence_from_string
from writer_utils import (
    write_cif,
    write_gro,
    write_lammpsdata,
    write_lammpstrj,
    write_xyz,
)

parser = argparse.ArgumentParser(
    description="Generate input files for diaphite simulations."
)
group = parser.add_mutually_exclusive_group(required=False)
group.add_argument(
    "--seq_file",
    type=str,
    default="seq.txt",
    help="Name of the file containing the sequence",
)
group.add_argument("--seq", type=str, help="Sequence of integers 1-4 to use")

parser.add_argument(
    "--out_file",
    type=str,
    default="diaphite.data",
    help="Name of the file to output to. Accepts *.xyz, *.cif, *.lammpstrj, *.gro or *.data",
)

parser.add_argument(
    "--nx", type=int, default=1, help="Number of unit cell repeats in x direction"
)

parser.add_argument(
    "--ny", type=int, default=1, help="Number of unit cell repeats in y direction"
)

parser.add_argument(
    "--nz", type=int, default=1, help="Number of unit cell repeats in z direction"
)

# Add scale factor command line option
parser.add_argument(
    "--sfa", type=float, default=1, help="Scale factor a"
)

parser.add_argument(
    "--sfb", type=float, default=1, help="Scale factor b"
)

parser.add_argument(
    "--sfc", type=float, default=1, help="Scale factor c"
)


args = parser.parse_args()

# Default cell constants.
CELL_A = 6.06885 * args.sfa
CELL_B = 2.50020 * args.sfb
CELL_C = 4.4279 * args.sfc

def generate_diamond_layer(
    cell_a: float = CELL_A,
    cell_b: float = CELL_B,
    cell_c: float = CELL_C,
    zstack: float = 0.0,
) -> np.array:
    """
    Generate atomic positions for one layer of type 1.

    A layer of type 1 is a diamond-type layer, and can be followed
    by a layer of type 1 (another diamond) or type 2 (a diamond / graphite
                                                      interface)
    Parameters
    ----------
    cell_a
        The length of the a unit cell vector, in Angstroms
    cell_b
        The length of the b unit cell vector, in Angstroms
    cell_c
        The length of the c unit cell vector, in Angstroms
    zstack
        The height we have stacked up to so far, in Angstroms
        i.e. the bottom of this layer
    Returns
    -------
    positions
        The positions of atoms within this layer, a 12x3 numpy array.
    """
    positions = np.array(
        [
            [0.470145, -0.499993, 0.082607],
            [0.219914, -0.499993, 0.085752],
            [-0.107809, 0.000007, 0.244577],
            [0.142256, 0.000007, 0.251337],
            [-0.184199, -0.499992, 0.408416],
            [-0.435601, -0.499992, 0.411123],
            [0.485617, 0.000008, 0.575337],
            [0.237729, 0.000008, 0.584477],
            [-0.089248, -0.499992, 0.736977],
            [0.163055, -0.499992, 0.753949],
            [-0.171446, 0.000008, 0.890393],
            [-0.419689, 0.000008, 0.899108],
        ]
    )
    positions[:, 0] *= cell_a
    positions[:, 1] *= cell_b
    positions[:, 2] *= cell_c
    positions[:, 2] += zstack

    return positions


def generate_diamond_graphene_layer(
    cell_a: float = CELL_A,
    cell_b: float = CELL_B,
    cell_c: float = CELL_C,
    zstack: float = 0.0,
) -> np.array:
    """
    Generate atomic positions for one layer of type 2.

    A layer of type 1 is a diamond-graphene interface layer, and can be followed
    only by a layer of type 3 (graphene)
    Parameters
    ----------
    cell_a
        The length of the a unit cell vector, in Angstroms
    cell_b
        The length of the b unit cell vector, in Angstroms
    cell_c
        The length of the c unit cell vector, in Angstroms
    zstack
        The height we have stacked up to so far, in Angstroms
        i.e. the bottom of this layer
    Returns
    -------
    positions
        The positions of atoms within this layer, a 5x3 numpy array.
    """

    positions = np.array(
        [
            [-0.477407, -0.499992, 0.080036],
            [0.25363, -0.499992, 0.102710],
            [-0.110947, 0.000007, 0.223845],
            [0.159511, 0.000008, 0.261057],
            [-0.235384, -0.499993, 0.335202],
        ]
    )

    positions[:, 0] *= cell_a
    positions[:, 1] *= cell_b
    positions[:, 2] *= cell_c
    positions[:, 2] += zstack

    return positions


def generate_graphene_layer(
    cell_a: float = CELL_A,
    cell_b: float = CELL_B,
    cell_c: float = CELL_C,
    zstack: float = 0.0,
) -> np.array:
    """
    Generate atomic positions for one layer of type 3.

    A layer of type 3 is a graphene layer, and can only be followed
    by another graphene layer or a graphene-diamond interface.
    Parameters
    ----------
    cell_a
        The length of the a unit cell vector, in Angstroms
    cell_b
        The length of the b unit cell vector, in Angstroms
    cell_c
        The length of the c unit cell vector, in Angstroms
    zstack
        The height we have stacked up to so far, in Angstroms
        i.e. the bottom of this layer
    Returns
    -------
    positions
        The positions of atoms within this layer, an 8x3 numpy array.
    """

    positions = np.array(
        [
            [0.220055, 0.000008, 0.168452],
            [-0.262165, -0.499992, 0.240856],
            [0.234004, -0.499992, 0.323156],
            [-0.265033, 0.000008, 0.396307],
            [0.253444, -0.499991, 0.649590],
            [-0.259395, 0.000008, 0.722740],
            [0.259402, 0.000009, 0.807303],
            [-0.253437, -0.499992, 0.880449],
        ]
    )

    positions[:, 0] *= cell_a
    positions[:, 1] *= cell_b
    positions[:, 2] *= cell_c
    positions[:, 2] += zstack
    return positions


def generate_graphene_diamond_layer(
    cell_a: float = CELL_A,
    cell_b: float = CELL_B,
    cell_c: float = CELL_C,
    zstack: float = 0.0,
) -> np.array:
    """
    Generate atomic positions for one layer of type 4.

    A layer of type 3 is a graphene  diamond interface layer,
    and can only be followed by a diamond layer
    Parameters
    ----------
    cell_a
        The length of the a unit cell vector, in Angstroms
    cell_b
        The length of the b unit cell vector, in Angstroms
    cell_c
        The length of the c unit cell vector, in Angstroms
    zstack
        The height we have stacked up to so far, in Angstroms
        i.e. the bottom of this layer
    Returns
    -------
    positions
        The positions of atoms within this layer, a 9x3 numpy array.
    """

    positions = np.array(
        [
            [0.26504, 0.000009, 0.163204],
            [-0.233995, -0.499992, 0.236350],
            [0.262172, -0.499991, 0.318650],
            [-0.220045, 0.000008, 0.391059],
            [0.23539, -0.499991, 0.650848],
            [-0.159505, 0.000008, 0.724989],
            [0.110953, 0.000008, 0.762205],
            [-0.253625, -0.499993, 0.883341],
            [0.477411, -0.499993, 0.906009],
        ]
    )

    positions[:, 0] *= cell_a
    positions[:, 1] *= cell_b
    positions[:, 2] *= cell_c
    positions[:, 2] += zstack
    return positions


def generate_diaphite_from_seq(
    layer_sequence: Iterable[LayerTypes],
) -> Tuple[np.array, np.array]:
    """
    Generate atomic positions for diaphite from a list of layer types.

    Parameters
    ----------
    layer_sequence
        An iterable of integers, each of which represents a layer type.

    Returns
    -------
    positions
        A Nx3 numpy array representing atomic positions in angstroms.
    simulation cell
        A 3x2 numpy array in the form [[xlo, xhi], [ylo, yhi], [zlo, zhi]].
    """
    layers = []
    zstack = 0.0
    for item in layer_sequence:
        if item == LayerTypes.Diamond:
            layers.append(generate_diamond_layer(zstack=zstack))
            zstack += CELL_C
        elif item == LayerTypes.DiamondGraphene:
            layers.append(generate_diamond_graphene_layer(zstack=zstack))
            zstack += CELL_C * 0.42655
        elif item == LayerTypes.Graphene:
            layers.append(generate_graphene_layer(zstack=zstack))
            zstack += CELL_C * 0.97053
        elif item == LayerTypes.GrapheneDiamond:
            layers.append(generate_graphene_diamond_layer(zstack=zstack))
            zstack += CELL_C * 0.99663
    positions = np.vstack(layers)
    simulation_cell = np.array([[0.0, CELL_A], [0.0, CELL_B], [0.0, zstack]])
    return positions, simulation_cell


def repeat_unit_cell(
    positions: np.array, simulation_cell: np.array, nx: int, ny: int, nz: int
):
    """
    Repeat the positions nx, ny and nz times in x, y, z directions.

    Leaves xlo, ylo and zlo unchanged and extends in the direction of xhi, yhi, zhi.
    Parameters
    ----------
    positions
        An Nx3 array of atomic positions.
    simulation cell
        A 3x2 numpy array in the form [[xlo, xhi], [ylo, yhi], [zlo, zhi]].
    nx
        Number of repeats in the x direction
    ny
        Number of repeats in the y direction
    nz
        Number of repeats in the z direction
    Returns
    -------
    positions
        An Nx3 array of repeated atomic positions
    simulation cell
        A 3x2 numpy array in the form [[xlo, xhi], [ylo, yhi], [zlo, zhi]].
    """
    assert nx >= 1, "nx must be a positive integer"
    assert ny >= 1, "ny must be a positive integer"
    assert nz >= 1, "nz must be a positive integer"
    cell_a = simulation_cell[0, 1] - simulation_cell[0, 0]
    cell_b = simulation_cell[1, 1] - simulation_cell[1, 0]
    cell_c = simulation_cell[2, 1] - simulation_cell[2, 0]

    positions = np.vstack(
        [positions + np.array([i * cell_a, 0.0, 0.0]) for i in range(nx)]
    )
    positions = np.vstack(
        [positions + np.array([0.0, j * cell_b, 0.0]) for j in range(ny)]
    )
    positions = np.vstack(
        [positions + np.array([0.0, 0.0, k * cell_c]) for k in range(nz)]
    )

    # Re-calculate the size of the simulation box
    simulation_cell[0, 1] = simulation_cell[0, 0] + (nx * cell_a)
    simulation_cell[1, 1] = simulation_cell[1, 0] + (ny * cell_b)
    simulation_cell[2, 1] = simulation_cell[2, 0] + (nz * cell_c)
    return positions, simulation_cell


def main():
    """
    Generate a diaphite sequence and write to a file.
    """

    if args.seq:
        sequence = sequence_from_string(args.seq)
    elif args.seq_file:
        sequence = read_sequence_from_file(args.seq_file)
    positions, simulation_cell = generate_diaphite_from_seq(sequence)

    positions, simulation_cell = repeat_unit_cell(
        positions, simulation_cell, args.nx, args.ny, args.nz
    )

    out_file = args.out_file
    if out_file.endswith(".xyz"):
        write_xyz(out_file, positions, simulation_cell, sequence)
    elif out_file.endswith(".cif"):
        write_cif(out_file, positions, simulation_cell)
    elif out_file.endswith(".data"):
        write_lammpsdata(out_file, positions, simulation_cell)
    elif out_file.endswith(".lammpstrj"):
        write_lammpstrj(out_file, positions, simulation_cell)
    elif out_file.endswith(".gro"):
        write_gro(out_file, positions, simulation_cell)
    else:
        raise RuntimeError(
            "Did not add a suffix of the form .xyz, .cif or .data to the output file."
        )
    print(f"Successfully wrote the diaphite positions to {out_file}")


if __name__ == "__main__":
    main()
