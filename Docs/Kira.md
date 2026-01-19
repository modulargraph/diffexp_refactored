# Kira - A Feynman Integral Reduction Program

Kira reduces Feynman integrals to a basis of master integrals using Integration-By-Parts (IBP) identities and Laporta's algorithm.

**Repository**: https://gitlab.com/kira-pyred/kira
**License**: GNU GPLv3
**Local installation**: `./kira/`

## Installation

### Pre-built Binaries (Easiest)

Download statically linked executables for Linux x86_64 from https://kira.hepforge.org

### Building from Source (macOS/Linux)

The Kira source code is already cloned in `./kira/`. To build:

#### Prerequisites

- **Platform**: Linux x86_64 or macOS
- **Compilers**: C++17 (C++) and C11 (C) compliant compilers
- **Build system**: Meson (0.46+) and Ninja (recommended)
- **Dependencies**:
  - GiNaC and CLN
  - zlib
  - yaml-cpp (auto-downloaded if missing)
  - **Fermat** (required at runtime) - set `FERMATPATH` environment variable

#### Build with Meson

```bash
cd kira
pip3 install --user meson  # if meson not installed
meson setup --prefix=$HOME/.local builddir
cd builddir
ninja
ninja install
```

#### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `-Dfirefly=false` | `true` | Disable FireFly finite field reconstruction |
| `-Dflint=false` | `true` | Disable FLINT (requires FireFly as subproject) |
| `-Dmpi=true` | `false` | Enable MPI for cluster parallelization |
| `-Djemalloc=true` | `false` | Enable jemalloc (often 20%+ speedup with FireFly) |
| `-Dweight_width=128` | `64` | Use 128-bit integers (creates `kira128` executable) |

## Usage

Kira uses YAML configuration files. A typical project structure:

```
my_reduction/
├── config/
│   ├── kinematics.yaml        # Kinematic variables and relations
│   └── integralfamilies.yaml  # Integral family definitions
├── jobs.yaml                  # Reduction job specification
└── myintegrals                # (optional) List of integrals to reduce
```

### Running Kira

```bash
kira jobs.yaml
```

### Configuration Files

#### kinematics.yaml

Defines momenta, invariants, and scalar product rules:

```yaml
kinematics:
  incoming_momenta: [p1, p2, p3, p4]
  outgoing_momenta: []
  momentum_conservation: [p4, -p1-p2-p3]
  kinematic_invariants:
    - [s, 2]   # [name, mass dimension]
    - [t, 2]
    - [m2, 2]
  scalarproduct_rules:
    - [[p1, p1], m2]
    - [[p2, p2], m2]
    - [[p1+p2, p1+p2], s]
  symbol_to_replace_by_one: m2  # optional: set one variable to 1
```

#### integralfamilies.yaml

Defines integral families (topologies):

```yaml
integralfamilies:
  - name: "box"
    loop_momenta: [k1]
    top_level_sectors: [15]  # binary: 1111 = all 4 propagators
    propagators:
      - ["k1", m2]           # propagator 1: 1/(k1^2 - m2)
      - ["k1+p1", 0]         # propagator 2: 1/(k1+p1)^2
      - ["k1+p1+p2", m2]     # propagator 3
      - ["k1+p1+p2+p3", 0]   # propagator 4
    # cut_propagators: [3, 4]  # optional: apply cuts
```

#### jobs.yaml

Specifies the reduction job:

```yaml
jobs:
  - reduce_sectors:
      reduce:
        - {sectors: [15], r: 6, s: 2}  # r = sum of positive indices, s = sum of negative
      select_integrals:
        select_mandatory_recursively:
          - {sectors: [15], r: 6, s: 2, d: 2}  # d = max dots
      run_initiate: true
      run_triangular: true
      run_back_substitution: true
```

## Key Job Options

### Reduction Control

| Option | Description |
|--------|-------------|
| `reduce_sectors` | Initiate IBP reduction |
| `run_symmetries` | Generate symmetry relations |
| `run_initiate` | Generate system of equations |
| `run_triangular` | Forward elimination |
| `run_back_substitution` | Back substitution |
| `run_firefly: true` | Use FireFly for finite field reconstruction |

### Integral Selection

| Option | Description |
|--------|-------------|
| `select_mandatory_recursively` | Select integrals by sector/r/s/d bounds |
| `select_mandatory_list` | Read integrals from file |
| `preferred_masters` | Specify preferred master integrals |

### Performance Options

| Option | Default | Description |
|--------|---------|-------------|
| `integral_ordering` | 1 | Choose from 8 orderings (1-8) |
| `iterative_reduction` | false | Reduce masterwise/sectorwise for memory efficiency |
| `factor_scan` | auto | Factor univariate factors (auto-enabled for 3+ variables) |
| `conditional` | true | Resume interrupted reductions |

### Output Options

| Option | Description |
|--------|-------------|
| `alt_dir` | Alternative output directory |
| `data_file: true` | Write human-readable results |

## Examples

The `kira/examples/` directory contains several examples:

- `1-loop-box/` - Simple 1-loop box integral
- `sunrise/` - Sunrise/banana diagram
- `double_box/` - 2-loop double box
- `double_pentagon/` - Double pentagon
- `tennis_court/` - Tennis court diagram

To run an example:

```bash
cd kira/examples/1-loop-box
kira jobs.yaml
```

## Integration with DiffExp

Kira outputs master integral reductions that can be used with DiffExp for differential equation solving:

1. Use Kira to reduce your Feynman integrals to master integrals
2. Derive differential equations for the master integrals
3. Use DiffExp to solve these differential equations via series expansion

## References

- Kira 1.0: [arXiv:1705.05610](https://arxiv.org/abs/1705.05610)
- Kira 2.0: [arXiv:2106.13989](https://arxiv.org/abs/2106.13989)
- Laporta's algorithm: [arXiv:hep-ph/0102033](https://arxiv.org/abs/hep-ph/0102033)
