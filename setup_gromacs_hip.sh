#!/usr/bin/env bash
#######################################################################
# GROMACS HIP Environment Setup Script
#
# This script sets up a complete conda environment for GROMACS
# with AMD GPU (HIP) support for protein-protein interaction analysis.
#
# What it does:
# 1. Creates a conda environment with all dependencies
# 2. Downloads and builds GROMACS 2025.4 with HIP support
# 3. Installs Python analysis packages
# 4. Configures environment for AMD MI210 GPU
#
# Requirements:
# - ROCm 6.x or 7.x installed
# - Conda/Mamba installed
# - ~10GB disk space
#
# Usage:
#   ./setup_gromacs_hip.sh              # Full setup
#   ./setup_gromacs_hip.sh --deps-only  # Only install dependencies
#   ./setup_gromacs_hip.sh --build-only # Only build GROMACS
#
# After setup:
#   conda activate gromacs_HIP
#   gmx_mpi --version
#
#######################################################################

set -e

#------------------------------------------------------------------------------
# CONFIGURATION
#------------------------------------------------------------------------------

# Environment name
ENV_NAME="gromacs_HIP"

# GROMACS version
GROMACS_VERSION="2025.4"
GROMACS_URL="https://ftp.gromacs.org/gromacs/gromacs-${GROMACS_VERSION}.tar.gz"

# Build directory
BUILD_DIR="/tmp/gromacs-build-hip"

# Installation prefix (inside conda env)
# Will be set after conda env is created

# Number of parallel jobs for compilation
NJOBS=$(nproc)

# ROCm path (adjust if different)
ROCM_PATH="${ROCM_PATH:-/opt/rocm}"

# Target GPU architecture (MI210 = gfx90a)
GPU_TARGETS="gfx90a"

#------------------------------------------------------------------------------
# HELPER FUNCTIONS
#------------------------------------------------------------------------------

log() { echo -e "\n\033[1;32m[$(date '+%H:%M:%S')]\033[0m $1"; }
log_warn() { echo -e "\n\033[1;33m[WARNING]\033[0m $1"; }
log_error() { echo -e "\n\033[1;31m[ERROR]\033[0m $1"; exit 1; }

check_rocm() {
    log "Checking ROCm installation..."
    
    if [ ! -d "$ROCM_PATH" ]; then
        log_error "ROCm not found at $ROCM_PATH. Please install ROCm first."
    fi
    
    if command -v rocm-smi &> /dev/null; then
        log "ROCm version: $(cat $ROCM_PATH/.info/version 2>/dev/null || echo 'unknown')"
        rocm-smi --showproductname 2>/dev/null | head -3 || true
    else
        log_warn "rocm-smi not found, but ROCm directory exists"
    fi
    
    # Check for HIP compiler
    if [ ! -f "$ROCM_PATH/bin/hipcc" ]; then
        log_error "hipcc not found. Please install ROCm HIP."
    fi
    
    log "HIP compiler: $($ROCM_PATH/bin/hipcc --version 2>&1 | head -1)"
}

check_conda() {
    log "Checking conda..."
    
    if ! command -v conda &> /dev/null; then
        log_error "Conda not found. Please install Miniconda or Anaconda."
    fi
    
    log "Conda: $(conda --version)"
    
    # Check for mamba, install if not present
    if ! command -v mamba &> /dev/null; then
        log "Installing mamba for faster package resolution..."
        conda install -y -c conda-forge mamba || log_warn "Could not install mamba, using conda"
    fi
    
    # Set package manager (prefer mamba)
    if command -v mamba &> /dev/null; then
        PKG_MGR="mamba"
    else
        PKG_MGR="conda"
    fi
    log "Package manager: $PKG_MGR"
}

#------------------------------------------------------------------------------
# INSTALL DEPENDENCIES
#------------------------------------------------------------------------------

install_dependencies() {
    log "Creating/updating conda environment: $ENV_NAME"
    
    # Check if environment exists
    if conda env list | grep -q "^${ENV_NAME} "; then
        log "Environment exists, updating..."
        CONDA_CMD="install"
    else
        log "Creating new environment..."
        conda create -n "$ENV_NAME" python=3.11 -y
        CONDA_CMD="install"
    fi
    
    # Activate environment
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    
    # Set package manager (prefer mamba if available)
    if command -v mamba &> /dev/null; then
        PKG_MGR="mamba"
    else
        PKG_MGR="conda"
    fi
    
    # Install build dependencies
    log "Installing build dependencies with $PKG_MGR..."
    $PKG_MGR install -y -c conda-forge \
        cmake \
        ninja \
        make \
        pkg-config \
        fftw \
        openmpi \
        openmpi-mpicc \
        libiconv
    
    # hwloc is optional - try to install, skip if not available
    log "Attempting to install hwloc (optional)..."
    $PKG_MGR install -y -c conda-forge hwloc 2>/dev/null || \
        pip install hwloc 2>/dev/null || \
        log_warn "hwloc not available, GROMACS will use system hwloc or disable"
    
    # Install Python dependencies for analysis and visualization
    log "Installing Python packages for analysis and visualization..."
    $PKG_MGR install -y -c conda-forge \
        numpy \
        scipy \
        pandas \
        matplotlib \
        seaborn \
        biopython \
        networkx \
        pillow \
        imageio \
        imageio-ffmpeg \
        scikit-learn \
        scikit-image \
        tqdm \
        pyyaml \
        h5py
    
    # Install visualization and analysis tools
    log "Installing visualization and analysis packages..."
    pip install --quiet \
        MDAnalysis \
        prolif \
        py3Dmol \
        nglview \
        moviepy \
        plotly \
        bokeh \
        pdb-tools \
        biotite
    
    # Install gnuplot for quick plotting
    log "Installing gnuplot and ffmpeg..."
    $PKG_MGR install -y -c conda-forge gnuplot ffmpeg 2>/dev/null || \
        log_warn "gnuplot/ffmpeg not available from conda, install manually if needed"
    
    # Get conda env path
    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME
    
    log "Dependencies installed to: $CONDA_PREFIX"
}

#------------------------------------------------------------------------------
# BUILD GROMACS
#------------------------------------------------------------------------------

build_gromacs() {
    log "Building GROMACS ${GROMACS_VERSION} with HIP support..."
    
    # Activate environment
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    
    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME
    
    # Create build directory
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"
    
    # Download GROMACS if not present
    if [ ! -f "gromacs-${GROMACS_VERSION}.tar.gz" ]; then
        log "Downloading GROMACS ${GROMACS_VERSION}..."
        wget -q "$GROMACS_URL" -O "gromacs-${GROMACS_VERSION}.tar.gz"
    fi
    
    # Extract
    if [ ! -d "gromacs-${GROMACS_VERSION}" ]; then
        log "Extracting..."
        tar -xzf "gromacs-${GROMACS_VERSION}.tar.gz"
    fi
    
    # Create build directory
    mkdir -p "gromacs-${GROMACS_VERSION}/build_hip"
    cd "gromacs-${GROMACS_VERSION}/build_hip"
    
    # Set up compiler paths
    export CC=$ROCM_PATH/bin/amdclang
    export CXX=$ROCM_PATH/bin/amdclang++
    export HIP_PATH=$ROCM_PATH
    export ROCM_PATH=$ROCM_PATH
    
    # Configure with CMake
    log "Configuring CMake..."
    cmake .. \
        -GNinja \
        -DCMAKE_INSTALL_PREFIX="$CONDA_PREFIX" \
        -DCMAKE_BUILD_TYPE=Release \
        -DGMX_BUILD_OWN_FFTW=OFF \
        -DGMX_FFT_LIBRARY=fftw3 \
        -DFFTWF_INCLUDE_DIR="$CONDA_PREFIX/include" \
        -DFFTWF_LIBRARY="$CONDA_PREFIX/lib/libfftw3f.so" \
        -DGMX_GPU=HIP \
        -DGMX_HIP_TARGET_ARCH="$GPU_TARGETS" \
        -DCMAKE_HIP_COMPILER="$ROCM_PATH/bin/amdclang++" \
        -DCMAKE_HIP_ARCHITECTURES="$GPU_TARGETS" \
        -DHIP_PATH="$ROCM_PATH" \
        -DGMX_MPI=ON \
        -DGMX_OPENMP=ON \
        -DGMX_SIMD=AVX2_256 \
        -DREGRESSIONTEST_DOWNLOAD=OFF \
        -DGMX_BUILD_HELP=OFF
    
    # Build
    log "Building (this may take 15-30 minutes)..."
    ninja -j${NJOBS}
    
    # Install
    log "Installing..."
    ninja install
    
    log "GROMACS installed to: $CONDA_PREFIX"
}

#------------------------------------------------------------------------------
# POST-INSTALL CONFIGURATION
#------------------------------------------------------------------------------

configure_environment() {
    log "Configuring environment..."
    
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    
    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME
    
    # Create activation script
    ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
    DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"
    
    mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"
    
    # Activation script
    cat > "$ACTIVATE_DIR/gromacs_hip.sh" << 'ACTIVATE_EOF'
# GROMACS HIP environment activation

# ROCm settings
export ROCM_PATH="${ROCM_PATH:-/opt/rocm}"
export PATH="$ROCM_PATH/bin:$PATH"
export LD_LIBRARY_PATH="$ROCM_PATH/lib:$LD_LIBRARY_PATH"

# AMD GPU settings (optimized for MI210)
export GMX_ENABLE_DIRECT_GPU_COMM=1
export GPU_MAX_HW_QUEUES=8
export HIP_VISIBLE_DEVICES=0
export ROCR_VISIBLE_DEVICES=0

# DO NOT set HSA_OVERRIDE_GFX_VERSION - causes kernel mismatch!

# GROMACS settings
export GMX_MAXBACKUP=-1  # Disable backup files

# OpenMP settings
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Source GROMACS environment if available
if [ -f "$CONDA_PREFIX/bin/GMXRC.bash" ]; then
    source "$CONDA_PREFIX/bin/GMXRC.bash"
fi

echo "GROMACS HIP environment activated"
echo "  GPU: AMD HIP (target: gfx90a)"
echo "  Use: gmx_mpi or mpirun -np 1 gmx_mpi mdrun ..."
ACTIVATE_EOF

    # Deactivation script
    cat > "$DEACTIVATE_DIR/gromacs_hip.sh" << 'DEACTIVATE_EOF'
# Cleanup GROMACS HIP environment
unset GMX_ENABLE_DIRECT_GPU_COMM
unset GPU_MAX_HW_QUEUES
DEACTIVATE_EOF

    chmod +x "$ACTIVATE_DIR/gromacs_hip.sh"
    chmod +x "$DEACTIVATE_DIR/gromacs_hip.sh"
    
    log "Activation scripts created"
}

verify_installation() {
    log "Verifying installation..."
    
    eval "$(conda shell.bash hook)"
    conda activate "$ENV_NAME"
    
    CONDA_PREFIX=$(conda info --base)/envs/$ENV_NAME
    
    # Check GROMACS binary
    if [ -f "$CONDA_PREFIX/bin/gmx_mpi" ]; then
        log "GROMACS binary found: $CONDA_PREFIX/bin/gmx_mpi"
        "$CONDA_PREFIX/bin/gmx_mpi" --version 2>&1 | head -5
    else
        log_warn "gmx_mpi not found!"
        return 1
    fi
    
    # Check GPU detection
    log "Checking GPU detection..."
    "$CONDA_PREFIX/bin/gmx_mpi" mdrun -version 2>&1 | grep -i "gpu\|hip\|rocm" || true
    
    # Check Python imports
    log "Checking Python packages..."
    python3 -c "
import numpy
import scipy
import pandas
import matplotlib
import MDAnalysis
print('All Python packages OK')
"
    
    log "✓ Installation verified successfully!"
    
    echo ""
    echo "============================================================"
    echo " Setup Complete!"
    echo "============================================================"
    echo ""
    echo "To use GROMACS:"
    echo "  conda activate $ENV_NAME"
    echo "  gmx_mpi --version"
    echo ""
    echo "To run simulations with GPU:"
    echo "  mpirun -np 1 gmx_mpi mdrun -deffnm md -nb gpu -pme gpu -update gpu"
    echo ""
    echo "HIP GPU limitations:"
    echo "  - Use '-bonded cpu' (GPU bonded not implemented)"
    echo "  - Do NOT set HSA_OVERRIDE_GFX_VERSION"
    echo ""
}

#------------------------------------------------------------------------------
# MAIN
#------------------------------------------------------------------------------

main() {
    echo ""
    echo "============================================================"
    echo " GROMACS HIP Environment Setup"
    echo "============================================================"
    echo ""
    echo "Environment: $ENV_NAME"
    echo "GROMACS: $GROMACS_VERSION"
    echo "ROCm: $ROCM_PATH"
    echo "GPU Target: $GPU_TARGETS"
    echo ""
    
    case "${1:-}" in
        --deps-only)
            check_conda
            install_dependencies
            ;;
        --build-only)
            check_rocm
            build_gromacs
            configure_environment
            verify_installation
            ;;
        --verify)
            verify_installation
            ;;
        *)
            check_conda
            check_rocm
            install_dependencies
            build_gromacs
            configure_environment
            verify_installation
            ;;
    esac
}

main "$@"
