
## Run Moments on Gawain
#TODO - Somehow threads go crazy when projection is more than 50 & 50
#TODO - fix R entering or just call a script - prefer entering R
#TODO - Fix conda path as variable or remove
#TODO - Create second script for fitting best model(s)
#		1. Take in SFS & pop names as well as model name and top parameters
#		2. Fit GOF for model
#		3. Refit model with BS of SFS to generate CIs
#		4. Make plot from moments package of specified model diagram

#Inputs
condENV=${1} #name of conda environment to use
sfsPath=${2} #path to easySFS.py
sfsPlot=${3} #path to 2dAFS_fold.py
momentsPipelinePath=${4} #path containing all the *.py optimization model scripts

outDir=$(pwd)/${5}

vcf=$(pwd)/${6} #vcf file to use
popSet=$(pwd)/${7} #file of population to individual translation

#Population names must match - extract from popset?
pop1=${8} #name of population 1
pop2=${9} #name of population 2
project_pop1=${10} #Number of indivuals in pop1 to project to
project_pop2=${11} #Number of indivuals in pop1 to project to

#Would like to pass through number of rounds/replicates to include (?)

THREADS=${12}

#Set up
current_dir=$(pwd)
mkdir -p ${outDir}/p${project_pop1}.p${project_pop2}
#conda activate ${condENV}

#Create SFS and project to specified number of individuals
${sfsPath} -i ${vcf} \
  -p ${popSet} \
  --proj ${project_pop1},${project_pop2} \
  -o ${outDir}/p${project_pop1}.p${project_pop2}/sfs
the_sfs=${outDir}/p${project_pop1}.p${project_pop2}/sfs/dadi/${pop1}-${pop2}.sfs

#Plot SFS
cd ${outDir}/p${project_pop1}.p${project_pop2}
${sfsPlot} ${the_sfs} ${pop1} ${pop2} ${project_pop1} ${project_pop2}
cd ${current_dir}

#Run moments
mkdir -p ${outDir}/p${project_pop1}.p${project_pop2}/momentsOptimize
cd ${outDir}/p${project_pop1}.p${project_pop2}/momentsOptimize
ls ${momentsPipelinePath}/moments_Run_2D_??_*py | parallel --no-notice -j $THREADS "python {} $the_sfs $pop1 $pop2" > optimizations.log 2> optimizations.err
cd ${current_dir}

#Process moments models and output summary graphs
Rscript scripts/momentsSummarise_ML.R ${outDir}/p${project_pop1}.p${project_pop2}/momentsOptimize ${outDir}/p${project_pop1}.p${project_pop2}
