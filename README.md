# sin_wave_wing_tail

## Project Overview

This project generates and analyzes sinusoidal wing/tail curves using Python.  
All code relies on common scientific computing libraries.

## How to Clone the Repository

To get a copy of this project on your local machine, run:

```bash
git clone https://github.com/Fgkl10/sin_wave_wing_tail.git
cd sin_wave_wing_tail
```

## Setting Up the Environment and Installing Dependencies

It is recommended to use a virtual environment to isolate and manage project dependencies. You can use either [venv](https://docs.python.org/3/library/venv.html) or [conda](https://docs.conda.io/).

### 1. Using venv (Recommended)

```bash
# 1. Create a virtual environment
python -m venv venv

# 2. Activate the virtual environment
# On Windows
venv\Scripts\activate
# On macOS/Linux
source venv/bin/activate

# 3. Install dependencies
pip install -r requirements.txt
```

### 2. Using conda

```bash
# 1. Create and activate a virtual environment
conda create -n myenv python=3.11
conda activate myenv

# 2. Install dependencies
pip install -r requirements.txt
```

### 3. Deactivating the Virtual Environment

- venv:  
  ```bash
  deactivate
  ```
- conda:  
  ```bash
  conda deactivate
  ```

---

For more help, please refer to the official documentation of each library or open an Issue in this repository.
