#!/bin/bash

#SBATCH --partition=q128
#SBATCH --mem=750GB
#SBATCH --time=2-00:00

##SBATCH --qos=qosqfat
##SBATCH --partition=qfat
##SBATCH --mem=0GB
##SBATCH --time=14-00:00

#SBATCH --job-name=Matrix-CNOT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaco@cs.au.dk
#SBATCH --exclusive

# Defaults:

QUBITS=6    # Number of Qubits
NAUTY=1     # using Nauty
EXTRA=1     # size of hash-tables for levels + EXTRA bits
DIST=-1     # maximum distance to try (-1 is unlimited)
POLY=0      # Compute polynomial coefficients
SWAP=0      # Swaps-for-free
BEAT=0      # Heart-beat every BEAT seconds
MAX=34      # max table size 2^MAX

while getopts B:D:E:M:N:P:Q:S:T:h flag
do
    case "${flag}" in
        B) BEAT=${OPTARG};;
        D) DIST=${OPTARG};;
        E) EXTRA=${OPTARG};;
        M) MAX=${OPTARG};;
        N) NAUTY=${OPTARG};;
        P) POLY=${OPTARG};;
        Q) QUBITS=${OPTARG};;
        S) SWAP=${OPTARG};;
        T) export OMP_NUM_THREADS=${OPTARG};;
        h) echo "Usage: matrix_cnot.sh [options] [goal]: optimal CNOT synthesis" 
           echo
           echo "Compile-time Options:"
           echo "  -B beat    : heart-beat every BEAT seconds (0=no beat) (default $BEAT)"
           echo "  -E extra   : size of hash-tables for levels + EXTRA bits (default $EXTRA)"
           echo "  -M max     : max table size 2^MAX (default $MAX)"
           echo "  -N nauty   : using Nauty (0 no, 1 yes) (default $NAUTY)"
           echo "  -P poly    : compute polynomial coefficients (0 no, 1 yes) (default $POLY)"
           echo "  -Q qubits  : number of Qubits (default $QUBITS)"
           echo "  -S swap    : swaps-are-for-free, requires nauty (0 no, 1 yes) (default $SWAP)"
           echo "  -T threads : number of OpenMP threads to use (\"\" is all cores) (default \"$OMP_NUM_THREADS\")"
           echo "  -h         : this help"
           echo
           echo "Run-time Options:"
           echo "  -D dist    : maximum distance to try (-1 is unlimited) (default $DIST)"
           echo "  goal       : filename for Goal matrix (see 'Inputs/' for examples)"
           echo
           echo "Note: all options are compile-time options, except for dist and goal."
           echo "Compiles a binary \"./matrix_cnot<QUBITS>\" and runs it."
           echo "Use ./matrix_cnot<Q> -<dist> goal" to rerun with different goal/distance
           exit 0;;
    esac
done

exec=matrix_cnot${QUBITS}
opts="-DN=$QUBITS -DE=$EXTRA -DMAX=$MAX -DPOLY=$POLY -DNAUTY=$NAUTY -DSWAP=$SWAP -DBEAT=$BEAT"
args="-fopenmp -O3 -DNDEBUG -march=native"
if [ $NAUTY -eq 1 ]; then
    nauty_args="-Inauty nauty/nautyW1.a -DWORDSIZE=32 -DMAXN=WORDSIZE -Wno-attributes"
fi
shift $((OPTIND - 1))
goal=$1

\rm -f $exec
set -x
g++ -o $exec matrix_cnot.cpp $opts $args $nauty_args
./$exec -$DIST $goal | tee matrix_cnot$QUBITS.txt