# Python runtime

The executed environment was an isolated venv created with Python 3.12.13. `CONDA_PYTHON_EXPLICIT.txt` fixes the 26 installed Conda packages in the interpreter's dependency closure; unrelated tools in the original base environment are excluded. `CONDA_PYTHON_PACKAGE_RECORDS.json` and `CONDA_PYTHON_SHA256.tsv` record original source URLs, builds and package hashes. `PYTHON_REQUIREMENTS_FREEZE_ALL.txt` fixes every installed venv package version, including packaging tools.

One portable restoration route, from this directory and using new directories, is:

```sh
export MAMBA_ROOT_PREFIX="$PWD/mamba_root"
micromamba create -y -p ./python_base --file CONDA_PYTHON_EXPLICIT.txt
./python_base/bin/python -m venv ./pydeseq2_env
./pydeseq2_env/bin/python -m pip install -r PYTHON_REQUIREMENTS_FREEZE_ALL.txt
./pydeseq2_env/bin/python -c 'import sys,pydeseq2; assert sys.version_info[:3] == (3,12,13); assert pydeseq2.__version__ == "0.5.4"'
```

This clean restoration was not executed during the revision. The existing environment was used and its versions checked before fitting. Conda hashes are installed-package metadata; original pip wheels were not preserved, so the pip freeze locks versions rather than exact wheel bytes. No dependency installation was performed for this sensitivity analysis.
