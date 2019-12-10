### Setting up CUDA Drivers

To run `signatureanalyzer` using available GPUs with Pytorch, one needs to install and set up CUDA. We provide the following set up code for setting up Cuda 10.1 using `ubuntu1804`. Once this is set up and a GPU is attached to the instance, an easy way to check is running `watch nvidia-smi`. For more questions, please see the NVIDIA website.

```
# -----------------------------
# Install CUDA for Tensorflow/Pytorch
# -----------------------------

# Dependencies
sudo apt-get install -y build-essential \
    dkms \
    freeglut3 \
    freeglut3-dev \
    libxi-dev \
    libxmu-dev

# Install CUDA
sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/cuda-ubuntu1804.pin
sudo mv cuda-ubuntu1804.pin /etc/apt/preferences.d/cuda-repository-pin-600
sudo apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
sudo add-apt-repository "deb http://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/ /"
sudo apt-get update
sudo apt-get -y install cuda

# Add to path
export PATH=$PATH:/usr/local/cuda-10.1/bin
export CUDADIR=/usr/local/cuda-10.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64
```
