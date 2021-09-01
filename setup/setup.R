# to run R 4.0.3:
# conda activate arabidopsis_test
# sudo /home/paulati/anaconda3/envs/arabidopsis_test/bin/R

R.version
#conda install -c conda-forge r-base=4.0.3

install.packages("BiocManager", version = '3.12')

BiocManager::install("RSQLite")
BiocManager::install("GO.db", version = '3.12')
BiocManager::install("GOSemSim")
BiocManager::install("topGO")
BiocManager::install("biomaRt")
BiocManager::install("ViSEAGO")

