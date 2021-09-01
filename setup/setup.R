R.version
#conda install -c conda-forge r-base=4.0.3

install.packages("BiocManager", version = '3.12')

BiocManager::install("RSQLite")
BiocManager::install("GO.db", version = '3.12')
BiocManager::install("GOSemSim")
BiocManager::install("topGO")

# to run R 4.0.3:
# conda activate arabidopsis_test
# sudo /home/paulati/anaconda3/envs/arabidopsis_test/bin/R


install.packages("httr")
install.packages("curl")


install.packages("XML")

BiocManager::install("BiocFileCache")



BiocManager::install("biomaRt")
BiocManager::install("ViSEAGO")



In .inet_warning(msg) :
  installation of package ‘curl’ had non-zero exit status
2: In .inet_warning(msg) :
  installation of package ‘XML’ had non-zero exit status
3: In .inet_warning(msg) :
  installation of package ‘httr’ had non-zero exit status
4: In .inet_warning(msg) :
  installation of package ‘BiocFileCache’ had non-zero exit status
5: In .inet_warning(msg) :
  installation of package ‘biomaRt’ had non-zero exit status

