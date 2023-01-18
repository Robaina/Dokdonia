# Download diffexpr source repo
echo "Downloading diffexpr repo"
git clone https://github.com/wckdouglas/diffexpr.git
cd diffexpr

# Create conda environment containing all dependencies
echo "Creating conda environment"
conda env create -q -n dokdonia -f environment.yml
source activate dokdonia

# Install diffexpr
echo "Installing diffexpr"
Rscript setup.R
python setup.py install

# Install R packages
echo "Installing R packages"
cd ..
Rscript setup.R

# Delete diffexr repo
echo "Deleting diffepr repo"
rm -fr diffexpr

echo "All done!"
