
git clone https://github.com/paulati/arabidopsis_phospho

cd ./arabidopsis_phospho/setup

conda remove --yes --name arabidopsis_test --all

conda create --yes --name arabidopsis_test python=3.7

conda install --yes --name arabidopsis_test --file requirements.txt

conda install -c conda-forge --yes --name arabidopsis_test --file requirements_conda-forge.txt

conda activate arabidopsis_test


