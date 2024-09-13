#!/bin/bash

###SBATCH -J
#SBATCH -N 1
##SBATCH -p intel_q
#SBATCH -p intel_preemptable_q
##SBATCH --cpus-per-task=2
##SBATCH --mem=220G
#SBATCH --tasks-per-node=8
#SBATCH -t 90:00:00
#SBATCH --account=personal
#SBATCH --account=nmayhall_group
##SBATCH --account=nmayhall_group-paid
#SBATCH --exclusive


#export OMP_NUM_THREADS=4

sleep 10
#source /projects/nmayhall_lab/arnab/dmrg/bimetallics/myenv/bin/activate
hostname

module reset
module load Python/3.10.8-GCCcore-12.2.0
#module load Python/3.8.6-GCCcore-10.2.0
#module load OpenMPI/4.0.5-GCC-10.2.0
module load OpenMPI/4.1.4-GCC-12.2.0

#module load intel/2022b
#module load Julia/1.7.2-linux-x86_64
#module purge
#module load gcc/9.2.0
export MKL_DEBUG_CPU_TYPE=5
export MKL_ENABLE_INSTRUCTIONS=AVX2
export PYSCF_TMPDIR=/projects/nmayhall_lab/arnab/dmrg/block2-example-data/05-Fe2OCl6/runs/stdmrg-5
export OMP_NUM_THREADS=4
which python3
python3 --version
python3 -c "import pyscf; print(pyscf.__version__)"
python3 -c "import pyscf; print(pyscf.__file__)"
python3 -c "import pyscf; print(pyscf.lib.param.TMPDIR)"
python3 -c "import block2; print(block2.__file__)"
python3 -c "import pyblock2; print(pyblock2.__file__)"

echo SLURM_TASKS_PER_NODE=$SLURM_TASKS_PER_NODE
echo OMP_NUM_THREADS=$OMP_NUM_THREADS
echo SLURM_JOBID=$SLURM_JOBID
echo SLURM_JOB_NAME=$SLURM_JOB_NAME
echo HOST_NAME = $(hostname)
echo PWD = $(pwd)
echo SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

export XRUN=orterun
export PYSCF_MPIPREFIX="$XRUN --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$SLURM_CPUS_PER_TASK"
export CPUTYPE=$(lscpu | grep 'Vendor ID' | awk '{print $3}')
if [ "$CPUTYPE" = "AuthenticAMD" ]; then
    echo "int mkl_serv_intel_cpu_true() { return 1; }" > fixcpu.c
    $CC -shared -fPIC -o libfixcpu.so fixcpu.c
    export LD_PRELOAD=$PWD/libfixcpu.so
fi

if [ "${SLURM_CPUS_PER_TASK}" != "" ]; then
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
fi

if [ "0" = "1" ]; then
    SCPT=dmrg
    if [ "0" = "1" ]; then
        SCPT=dmrg-rev
    fi
else
    SCPT=29_cr2_morokuma_highspin2
fi

TJ=$(echo ${SCPT}.out.* | tr ' ' '\n' | grep '\*$' -v | wc -l)
export TJ=$(expr ${TJ} + 1)
echo ${SCPT}.out.${TJ} >> OUTFILE
echo $SLURM_JOBID >> JOBIDS

which $XRUN

if [ "$?" = "1" ] || [ "${SLURM_TASKS_PER_NODE}" = "" ] || [ "0" = "1" ]; then
    if [ "0" = "1" ]; then
        [ -f ./FCIDUMP ] && rm ./FCIDUMP
        ln -s /projects/nmayhall_lab/arnab/dmrg/bimetallics/cr2_morokuma/FCIDUMP ./FCIDUMP
        cp /projects/nmayhall_lab/arnab/dmrg/bimetallics/cr2_morokuma/stdmrg-5/${SCPT}.conf ${SCPT}.conf.${TJ}
        python3 -u $(which block2main) ${SCPT}.conf.${TJ} > ${SCPT}.out.${TJ}
    else
        python3 -u ${SCPT}.py 0 > ${SCPT}.out.${TJ}
    fi
else
    if [ "0" = "1" ]; then
        [ -f ./FCIDUMP ] && rm ./FCIDUMP
        ln -s /projects/nmayhall_lab/arnab/dmrg/bimetallics/cr2_morokuma/FCIDUMP ./FCIDUMP
        cp /projects/nmayhall_lab/arnab/dmrg/bimetallics/cr2_morokuma/${SCPT}.conf ${SCPT}.conf.${TJ}
        if [ "$XRUN" = "srun" ]; then
            srun python3 -u $(which block2main) ${SCPT}.conf.${TJ} > ${SCPT}.out.${TJ}
        else
            $XRUN --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS \
                python3 -u $(which block2main) ${SCPT}.conf.${TJ} > ${SCPT}.out.${TJ}
        fi
    else
        if [ "$XRUN" = "srun" ]; then
            srun python3 -u ${SCPT}.py 0 > ${SCPT}.out.${TJ}
        else
            $XRUN --map-by ppr:$SLURM_TASKS_PER_NODE:node:pe=$OMP_NUM_THREADS \
                python3 -u ${SCPT}.py 0 > ${SCPT}.out.${TJ}
        fi
    fi
fi

if [ "$?" = "0" ]; then
    echo "SUCCESSFUL TERMINATION"
else
    echo "ERROR TERMINATION"
fi

if [ "0" = "1" ]; then
    cp /projects/nmayhall_lab/arnab/dmrg/bimetallics/cr2_morokuma/node0/1pdm.npy ${SCPT}.1pdm.${TJ}.npy
fi
