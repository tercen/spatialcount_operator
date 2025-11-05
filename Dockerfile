FROM tercen/runtime-python39:0.1.0

COPY . /operator
WORKDIR /operator

ENV PYTHONPATH "${PYTHONPATH}:~/.pyenv/versions/3.9.0/bin/python3"

# Install scimap first (this will install compatible numpy, pandas, scipy, etc.)
RUN python3 -m pip install --no-cache-dir scimap>=2.0.0

# Install tercen client's direct dependencies that aren't already installed
RUN python3 -m pip install --no-cache-dir git+https://github.com/tercen/pytson@1.8.10

# Now install tercen client (it will use the already-installed numpy/pandas/scipy from scimap)
RUN python3 -m pip install --no-cache-dir --no-deps git+https://github.com/tercen/tercen_python_client@0.12.12

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT ["python3", "main.py"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
