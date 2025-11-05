FROM tercen/runtime-python39:0.1.0

COPY . /operator
WORKDIR /operator

ENV PYTHONPATH "${PYTHONPATH}:~/.pyenv/versions/3.9.0/bin/python3"

# Install scimap first (this will install compatible numpy, pandas, scipy, etc.)
RUN python3 -m pip install --no-cache-dir scimap>=2.0.0

# Install tercen client's direct dependencies (without their sub-dependencies)
RUN python3 -m pip install --no-cache-dir \
    git+https://github.com/tercen/pytson@1.8.10 \
    polars \
    pytz \
    python-dateutil

# Now install tercen client with --no-deps (uses already-installed packages)
RUN python3 -m pip install --no-cache-dir --no-deps git+https://github.com/tercen/tercen_python_client@0.12.12

ENV TERCEN_SERVICE_URI https://tercen.com

ENTRYPOINT ["python3", "main.py"]
CMD ["--taskId", "someid", "--serviceUri", "https://tercen.com", "--token", "sometoken"]
