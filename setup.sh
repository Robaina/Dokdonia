# Create conda environment containing all dependencies
echo "Creating conda environment"
conda env create -q -n dokdonia -f environment.yml
source activate dokdonia

# Install R packages
echo "Installing R packages"
Rscript setup.R

# Download diffexpr source repo
echo "Downloading diffexpr repo"
git clone https://github.com/wckdouglas/diffexpr.git

# Install diffexpr
echo "Installing diffexpr"
cd diffexpr
Rscript setup.R
python setup.py -q install
cd ..

# Delete diffexr repo
echo "Deleting diffepr repo"
rm -fr diffexpr

# Install Dokdonia's package
echo "Installing Dokdonia"
poetry build
pip install dist/dokdonia*.whl
rm -fr dist

echo "All done!"
