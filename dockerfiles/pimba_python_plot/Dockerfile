# Use Miniconda with an older Python version
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Copy the Python script into the container
COPY pimba_plot_graphs.py /app/pimba_plot_graphs.py

# Create a Conda environment with Python 3.9
RUN conda create -n plot_env python=3.9 -y && \
    conda install -n plot_env -y conda-forge::matplotlib \
    anaconda::scipy anaconda::scikit-bio anaconda::seaborn && \
    conda clean --all -y

# Activate the environment
SHELL ["conda", "run", "-n", "pimba_env", "/bin/bash", "-c"]

# Set the default command to open a shell
CMD ["/bin/bash"]
