#!/bin/sh
# run-all.sh -- R, then Python, then R, then R.
#
# Usage:  RSA_PYTHON=/path/to/rsaenv/bin/python ./run-all.sh
# The default assumes a venv named `rsaenv` beside this script.
set -e
here=$(cd "$(dirname "$0")" && pwd)
python_bin=${RSA_PYTHON:-"$here/rsaenv/bin/python"}

if [ ! -x "$python_bin" ]; then
  echo "No Python at $python_bin." >&2
  echo "Create it with:" >&2
  echo "  uv venv -p 3.12 \"$here/rsaenv\"" >&2
  echo "  uv pip install --python \"$here/rsaenv/bin/python\" -r \"$here/requirements.txt\"" >&2
  echo "or point RSA_PYTHON at an interpreter with rsatoolbox==0.3.2." >&2
  exit 1
fi

Rscript "$here/01-fixture.R"
"$python_bin" "$here/02-rsatoolbox.py"
Rscript "$here/03-compare.R"
Rscript "$here/04-extension.R"
Rscript "$here/05-manifest.R"
