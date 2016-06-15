##install python
sudo apt-get install build-essential checkinstall
sudo apt-get install libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev
wget https://www.python.org/ftp/python/2.7.11/Python-2.7.11.tgz
tar -xvf Python-2.7.11.tgz
cd Python-2.7.11
./configure
make
sudo checkinstall
######################################################################
##install biopytho
sudo apt-get install python-dev libxml2-dev libxslt1-dev zlib1g-dev
sudo apt-get install python-pip
pip install --user biopython

