#PBS -N 00Prep
#PBS -l select=1:ncpus=4:mem=4gb:interconnect=hdr,walltime=72:00:00

cd ${PBS_O_WORKDIR}
module add anaconda3/5.1.0-gcc/8.3.1

source activate deep-learning
BASEDIR=$(pwd)
python ${BASEDIR}/00-DataPrep/Prepare_inputs.py -i ${BASEDIR}/Test/Test_exp.csv -p ${BASEDIR}/Test/Test