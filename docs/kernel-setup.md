# Running the notebook: kernel setup

## The symptom

Opening `notebooks/scope_map_walkthrough.ipynb` in VS Code or JupyterLab and running it,
the kernel dies on the first code cell. The front end reports only that the kernel crashed
and is restarting. There is no traceback and no Python error.

## The cause

NumPy **delay-loads** its BLAS library, `libblas.dll`, which lives in the environment's
`Library\bin`. Windows resolves a delay-loaded DLL through PATH. Activating the
environment is what puts `Library\bin` on PATH; launching `<env>\python.exe` directly does
not, so the load fails as a structured exception rather than a Python error, and the
process is terminated before it can raise.

Run with `-X faulthandler`, the abort is precise:

```
Windows fatal exception: code 0xc06d007f
  File "...\numpy\__init__.py", line 881 in blas_fpe_check
  File "...\numpy\__init__.py", line 890 in <module>
```

`0xc06d007f` is a delay-load dependency failure. Measured on this machine, from a shell
with no conda directory on PATH:

| Command | Result |
|---|---|
| `envs\spatiotemporal\python.exe -c "import numpy"` | exits `-1066598273` (0xC0000409), no output |
| same for `scipy`, `matplotlib` | identical abort |
| `import zmq`, `import ipykernel` | succeed |
| the same imports after `activate.bat` | succeed |
| after prepending `Library\bin` to `os.environ["PATH"]` | succeed |

That `zmq` and `ipykernel` load but NumPy does not is why the kernel appears to start
normally and then dies: it survives until the first cell imports NumPy.

This is the failure the workspace `CLAUDE.md` §9 warns about — "calling `env\python.exe`
directly leaves `Library\bin` off PATH, and NumPy, SciPy and OpenMM then abort the process
with no traceback".

**A fix that does not work.** Adding the directory with `os.add_dll_directory` leaves the
abort unchanged. That call affects `LoadLibraryEx` searches, whereas delay-loading goes
through PATH. Setting `os.environ["PATH"]` does work.

## The fix

```
conda activate spatiotemporal
python tools/install_kernel.py
```

This applies two independent fixes, because a front end can start a kernel in two
different ways.

**1. A kernelspec that activates the environment.** Written per user to
`%APPDATA%\jupyter\kernels\spatiotemporal\`, displayed as **Python (spatiotemporal)**. Its
`argv` points at a generated `launch_kernel.bat` that calls `activate.bat` and then starts
`ipykernel_launcher`.

The activation lives in a batch file rather than directly in `argv` because a kernelspec's
`argv` is passed to CreateProcess, which re-quotes each element; quotes around a path then
reach `cmd.exe` escaped as `\"` and it answers "is not recognized as an internal or
external command". Inside a batch file the quoting is ordinary, so paths with spaces work.

**2. A `.pth` inside the environment** that prepends `Library\bin` to PATH at interpreter
startup, before any import. This covers front ends that ignore the kernelspec and launch
the interpreter themselves — which is what VS Code does when the kernel is picked under
*Python Environments*, shown there as `spatiotemporal (Python 3.11.17)`.

The kernelspec lives outside the environment and survives rebuilding its packages. The
`.pth` lives inside it and does not, so re-run the script after recreating the
environment. `python tools/install_kernel.py --remove` undoes both.

### In VS Code

With both fixes in place either choice works. **Select Kernel → Python (spatiotemporal)**
under *Jupyter Kernel* is the one to prefer, because it does not depend on the `.pth`
surviving an environment rebuild.

`tools/make_notebook.py` writes `spatiotemporal` into the notebook's kernelspec metadata,
so a freshly generated notebook already asks for that kernel.

### Checking it worked

```
conda activate spatiotemporal
python tools/make_notebook.py --execute
```

This runs the whole notebook headless through the same kernelspec. It prints the cell that
failed if one does, and the notebook's last cell raises if any status check failed.

## Known limitation

Because the kernel runs under `cmd.exe`, "Interrupt Kernel" may not reach the Python
process on Windows; use "Restart Kernel" if a long cell has to be stopped. The stock
kernelspec has the same limitation whenever it is wrapped for activation, and it is the
usual trade for having the environment set up correctly.
