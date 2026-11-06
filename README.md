# cmbx-python-parser

Standalone Python library for reading sequence metadata and signal data from
Chromeleon `.cmbx` archives. The parser uses only the Python standard library.

## Requirements

Python 3.10 or later.

## Install

Clone this repository and install it into the active environment:

```powershell
python -m pip install .
```

Alternatively, copy `cmbx_parser.py` into a project. No additional runtime
packages are required.

## Quick start

```python
from cmbx_parser import CmbxParseError, read_cmbx, read_cmbx_channel

try:
    sequence = read_cmbx("sequence.cmbx")
    print(sequence.name)

    for injection in sequence.injections:
        print(injection.name, injection.signal_channels)

        channel = next(
            (
                item
                for item in injection.channels
                if item.item_type.endswith(".Signal") and item.raw_filename
            ),
            None,
        )
        if channel is not None:
            chromatogram = read_cmbx_channel(
                sequence, injection, channel.name
            )
            print(chromatogram.times_min[:3])
            print(chromatogram.values[:3])
            break
except CmbxParseError as error:
    print(f"Could not read CMBX data: {error}")
```

## API

- `read_cmbx(path)` reads sequence metadata and returns a `CmbxSequence`.
  `path` may be a string or `pathlib.Path`. Signal samples are decoded on
  demand.
- `read_cmbx_channel(sequence, injection, channel_name, *, signal_unit=None)`
  reads a channel's RAW data and returns a `Chromatogram`. Set `signal_unit`
  when the physical unit is known but not one of the built-in units.
- `read_cmbx_uv_vis_1`, `read_cmbx_emission_1`,
  `read_cmbx_pump_pressure`, and `read_cmbx_cc_temp` read those named channels
  with their known units.
- `CmbxSequence` exposes `name`, `source`, `injections`, and
  `signal_channels`. Each `CmbxInjection` exposes `name`, `url`, `channels`,
  and `signal_channels`. Each `CmbxChannel` describes the channel's name,
  item type, RAW filename, declared size, and URL.
- `Chromatogram` exposes `times_min`, `values`, `steps_s`, `source`, and
  `metadata`, along with `injection_name`, `channel`, and `signal_unit`.
- `CmbxParseError` is raised for invalid archives, missing data, unsupported
  encodings, and malformed or oversized signals.

`injection.signal_channels` lists signal names in header order, including
channels that may not have RAW data. Check `CmbxChannel.raw_filename` before
reading a channel. `sequence.signal_channels` lists sequence-wide
RAW-backed signal names.

Supported RAW point encodings are `PtsLL2Df`, `PtsLDiff`, and `PtsDDCmp`.
Built-in signal units are `mAU` for `UV_VIS_1`, `counts` for `Emission_1`,
`bar` for `Pump_Pressure`, and `°C` for `CC_Temp`. Units for other channels
are left unspecified unless supplied by the caller.

## Tests

Run the test suite with:

```powershell
python -m unittest discover -s tests -v
```
