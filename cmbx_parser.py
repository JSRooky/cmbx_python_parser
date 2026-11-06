"""Standalone parser for Chromeleon CMBX archives and RAW signal channels.

A CMBX file is a ZIP archive. Its ``header.xml`` maps injections and channels
to RAW members; each RAW member contains a signal header and encoded point
sections. This module reads that metadata and decodes the supported signal
formats using only the Python standard library.
"""

from __future__ import annotations

from dataclasses import dataclass
import lzma
import math
from pathlib import Path
import struct
from typing import Protocol
from zipfile import BadZipFile, ZipFile
from xml.etree import ElementTree

_SIGNAL_HEADER_TAG = b"SignHdr\0"
_PTS_LL2_TAG = b"PtsLL2Df"
_PTS_LDIFF_TAG = b"PtsLDiff"
_PTS_DDCMP_TAG = b"PtsDDCmp"
_PTS_LDIFF_HEADER_SIZE = 28
_PTS_LDIFF_VALUE_SCALE = 0.01
_MAX_RAW_BYTES = 256 << 20
_MAX_LZMA2_DICTIONARY = 64 << 20
_MAX_SIGNAL_POINTS = 1 << 25
_MAX_SIGNED_INTEGER = (1 << 63) - 1
_KNOWN_SIGNAL_UNITS = {
    "UV_VIS_1": "mAU",
    "Emission_1": "counts",
    "Pump_Pressure": "bar",
    "CC_Temp": "°C",
}


class CmbxParseError(ValueError):
    """Raised when a CMBX archive or its RAW signal is invalid."""


@dataclass(frozen=True)
class CmbxChannel:
    """Channel metadata from ``header.xml`` and its optional RAW member."""

    name: str
    item_type: str
    raw_filename: str | None
    size: int | None
    url: str


@dataclass(frozen=True)
class CmbxInjection:
    """An injection and its channels as listed in the CMBX header."""

    name: str
    url: str
    channels: tuple[CmbxChannel, ...]

    @property
    def signal_channels(self) -> tuple[str, ...]:
        """Return this injection's unique signal channel names in header order."""
        return tuple(
            dict.fromkeys(
                channel.name
                for channel in self.channels
                if channel.item_type.endswith(".Signal")
            )
        )


@dataclass(frozen=True)
class CmbxSequence:
    """Parsed sequence metadata; signal values are read on demand."""

    source: Path
    name: str
    injections: tuple[CmbxInjection, ...]
    archive_entries: frozenset[str]

    @property
    def signal_channels(self) -> tuple[str, ...]:
        """Return unique RAW-backed signal names in first-seen order."""
        return tuple(
            dict.fromkeys(
                channel.name
                for injection in self.injections
                for channel in injection.channels
                if channel.item_type.endswith(".Signal") and channel.raw_filename
            )
        )


@dataclass(frozen=True)
class Chromatogram:
    """A decoded signal with time in minutes and values in channel units."""

    source: Path
    metadata: dict[str, str]
    times_min: tuple[float, ...]
    values: tuple[float, ...]
    steps_s: tuple[float | None, ...]

    @property
    def injection_name(self) -> str:
        return self.metadata.get("Name", self.source.stem)

    @property
    def channel(self) -> str:
        return self.metadata.get("Channel", "Unknown channel")

    @property
    def signal_unit(self) -> str:
        return self.metadata.get("Signal Unit", "")


class _ChromatogramLike(Protocol):
    metadata: dict[str, str]

    @property
    def channel(self) -> str: ...


def _sequence_sample_key(url: str) -> tuple[str, str] | None:
    """Extract stable sequence/sample IDs from a Chromeleon object URL."""
    parts = url.split("/")
    for index, part in enumerate(parts[:-1]):
        if part.endswith(".seq") and parts[index + 1].endswith(".smp"):
            return part, parts[index + 1]
    return None


def export_matches_injection(
    chromatogram: _ChromatogramLike, injection: CmbxInjection
) -> bool:
    """Match exports across Chromeleon hosts by sequence and sample IDs."""
    export_url = chromatogram.metadata.get("URL", "")
    export_key = _sequence_sample_key(export_url)
    if export_key is None or export_key != _sequence_sample_key(injection.url):
        return False
    return any(
        channel.name == chromatogram.channel
        and _sequence_sample_key(channel.url) == export_key
        for channel in injection.channels
    )


def _read_signed_varint(data: bytes, offset: int) -> tuple[int, int]:
    """Read Chromeleon's sign-and-magnitude variable-length integer."""
    if offset >= len(data):
        raise CmbxParseError("RAW block ends inside an integer")

    first = data[offset]
    offset += 1
    negative = bool(first & 0x40)
    magnitude = first & 0x3F
    shift = 6
    while first & 0x80:
        if offset >= len(data):
            raise CmbxParseError("RAW block ends inside an integer")
        first = data[offset]
        offset += 1
        chunk = first & 0x7F
        if shift >= 63 or chunk > (_MAX_SIGNED_INTEGER - magnitude) >> shift:
            raise CmbxParseError("RAW integer is too large")
        magnitude |= chunk << shift
        shift += 7

    return (-magnitude if negative else magnitude), offset


def _decode_ddcmp_block(
    body: bytes, *, existing_point_count: int
) -> list[tuple[float, float]]:
    """Decode one length-prefixed raw LZMA2 block of little-endian float pairs."""
    offset = 0
    header: list[int] = []
    while len(header) < 4:
        value, offset = _read_signed_varint(body, offset)
        header.append(value)
    version, reserved, decompressed_size, compressed_size = header
    if (version, reserved) != (1, 0):
        raise CmbxParseError(
            f"Unsupported PtsDDCmp block header: {version}, {reserved}"
        )
    if (
        decompressed_size <= 0
        or decompressed_size % 16
        or decompressed_size > _MAX_SIGNAL_POINTS * 16
    ):
        raise CmbxParseError(
            "Invalid decompressed data size in PtsDDCmp block"
        )
    point_count = decompressed_size // 16
    if existing_point_count + point_count > _MAX_SIGNAL_POINTS:
        raise CmbxParseError("RAW file contains too many points")
    if compressed_size != len(body) - offset or compressed_size < 2:
        raise CmbxParseError(
            "Compressed data size in PtsDDCmp block does not match the header"
        )

    dictionary_property = body[offset]
    if dictionary_property > 40:
        raise CmbxParseError("Invalid LZMA2 dictionary parameter in PtsDDCmp")
    # LZMA2 stores its dictionary size in one property byte before the stream.
    dictionary_size = (
        0xFFFFFFFF
        if dictionary_property == 40
        else (2 | (dictionary_property & 1)) << (dictionary_property // 2 + 11)
    )
    if dictionary_size > _MAX_LZMA2_DICTIONARY:
        raise CmbxParseError("LZMA2 dictionary in PtsDDCmp exceeds the allowed size")

    try:
        decoder = lzma.LZMADecompressor(
            format=lzma.FORMAT_RAW,
            filters=[{"id": lzma.FILTER_LZMA2, "dict_size": dictionary_size}],
        )
        decompressed = decoder.decompress(
            body[offset + 1 :], max_length=decompressed_size + 1
        )
    except (lzma.LZMAError, ValueError) as error:
        raise CmbxParseError(f"Failed to decompress PtsDDCmp block: {error}") from error

    if (
        len(decompressed) != decompressed_size
        or not decoder.eof
        or decoder.unused_data
    ):
        raise CmbxParseError(
            "Decompressed PtsDDCmp block size does not match the declared size"
        )

    points = list(struct.iter_unpack("<dd", decompressed))
    if not all(
        math.isfinite(time_min) and math.isfinite(value)
        for time_min, value in points
    ):
        raise CmbxParseError("Non-numeric values found in PtsDDCmp block")
    return points


def _decode_signal(
    raw_data: bytes,
    injection: CmbxInjection,
    source: Path,
    channel_name: str,
    signal_unit: str,
) -> Chromatogram:
    """Parse a RAW member, decode its single point encoding, and normalize units."""
    if len(raw_data) > _MAX_RAW_BYTES:
        raise CmbxParseError("RAW file is too large to parse")

    offset = 0
    signal_name: str | None = None
    time_units: float | None = None
    value_units: float | None = None
    signal_encoding: bytes | None = None
    ticks: list[int] = []
    values_raw: list[int] = []
    diff_header: bytes | None = None
    diff_blocks: list[tuple[float, list[int]]] = []
    diff_point_count = 0
    ddcmp_points: list[tuple[float, float]] = []

    # RAW members are concatenated sections with an 8-byte tag and 8-byte size.
    while offset < len(raw_data):
        if offset + 16 > len(raw_data):
            raise CmbxParseError(f"RAW file ends in the middle of a section at offset {offset}")

        tag = raw_data[offset : offset + 8]
        section_size = struct.unpack_from("<Q", raw_data, offset + 8)[0]
        if section_size < 16 or offset + section_size > len(raw_data):
            raise CmbxParseError(
                f"Invalid RAW section size at offset {offset}: {section_size}"
            )
        body = raw_data[offset + 16 : offset + section_size]

        if tag == _SIGNAL_HEADER_TAG:
            if signal_name is not None:
                raise CmbxParseError("RAW file contains multiple SignHdr headers")
            if len(body) < 15 or body[:3] != b"\x01\x02\x07":
                raise CmbxParseError("Invalid SignHdr header in RAW file")
            name_length = body[14]
            try:
                signal_name = body[15 : 15 + name_length].decode("utf-8")
            except UnicodeDecodeError as error:
                raise CmbxParseError(
                    "Channel name in SignHdr has an invalid encoding"
                ) from error
            if len(signal_name.encode("utf-8")) != name_length:
                raise CmbxParseError("Channel name in SignHdr is truncated")
            if signal_name != channel_name:
                raise CmbxParseError(
                    f"RAW file contains channel {signal_name!r}, "
                    f"expected {channel_name!r}"
                )
        elif tag == _PTS_LL2_TAG:
            if signal_name is None:
                raise CmbxParseError("RAW file has point section data before SignHdr")
            if signal_encoding is not None and signal_encoding != tag:
                raise CmbxParseError("RAW file mixes point section formats")
            signal_encoding = tag
            if len(body) < 16:
                raise CmbxParseError(
                    "PtsLL2Df section does not contain data parameters"
                )

            block_time_units, block_value_units = struct.unpack_from("<dd", body)
            if (
                not math.isfinite(block_time_units)
                or not math.isfinite(block_value_units)
                or block_time_units <= 0
                or block_value_units <= 0
            ):
                raise CmbxParseError(
                    "Invalid time or signal units in PtsLL2Df"
                )
            if time_units is not None and (
                block_time_units != time_units or block_value_units != value_units
            ):
                raise CmbxParseError(
                    "RAW file changes units between PtsLL2Df blocks"
                )
            time_units = block_time_units
            value_units = block_value_units

            encoded: list[int] = []
            block_offset = 16
            while block_offset < len(body):
                integer, block_offset = _read_signed_varint(body, block_offset)
                encoded.append(integer)
            if not encoded or encoded.pop() != 0 or len(encoded) % 2:
                raise CmbxParseError(
                    "PtsLL2Df section must contain (time, value) pairs and a zero terminator"
                )
            point_count = len(encoded) // 2
            if len(ticks) + point_count > _MAX_SIGNAL_POINTS:
                raise CmbxParseError("RAW file contains too many points")

            # PtsLL2Df stores pairs of second differences after the first pair.
            tick = value = tick_delta = value_delta = 0
            for pair_index in range(point_count):
                time_part, value_part = encoded[pair_index * 2 : pair_index * 2 + 2]
                if pair_index == 0:
                    tick, value = time_part, value_part
                elif pair_index == 1:
                    tick_delta, value_delta = time_part, value_part
                    tick += tick_delta
                    value += value_delta
                else:
                    tick_delta += time_part
                    value_delta += value_part
                    tick += tick_delta
                    value += value_delta
                if not -(1 << 63) <= tick_delta <= _MAX_SIGNED_INTEGER:
                    raise CmbxParseError(
                        "RAW time step exceeds the allowed range"
                    )
                if not -(1 << 63) <= value_delta <= _MAX_SIGNED_INTEGER:
                    raise CmbxParseError(
                        "RAW signal step exceeds the allowed range"
                    )
                if not -(1 << 63) <= tick <= _MAX_SIGNED_INTEGER:
                    raise CmbxParseError(
                        "RAW time exceeds the allowed range"
                    )
                if not -(1 << 63) <= value <= _MAX_SIGNED_INTEGER:
                    raise CmbxParseError(
                        "RAW signal value exceeds the allowed range"
                    )
                ticks.append(tick)
                values_raw.append(value)
        elif tag == _PTS_LDIFF_TAG:
            if signal_name is None:
                raise CmbxParseError("RAW file has point section data before SignHdr")
            if signal_encoding is not None and signal_encoding != tag:
                raise CmbxParseError("RAW file mixes point section formats")
            signal_encoding = tag
            if len(body) <= _PTS_LDIFF_HEADER_SIZE:
                raise CmbxParseError("PtsLDiff section does not contain a signal value")

            # PtsLDiff stores a block start time, shared parameters, and
            # first differences of the signal values.
            start_time_min = struct.unpack_from("<d", body)[0]
            if not math.isfinite(start_time_min):
                raise CmbxParseError("Invalid start time in PtsLDiff block")
            block_header = body[8:_PTS_LDIFF_HEADER_SIZE]
            if diff_header is not None and block_header != diff_header:
                raise CmbxParseError("RAW file changes PtsLDiff block parameters")
            diff_header = block_header

            block_offset = _PTS_LDIFF_HEADER_SIZE
            first_value, block_offset = _read_signed_varint(body, block_offset)
            block_values = [first_value]
            if diff_point_count + len(block_values) > _MAX_SIGNAL_POINTS:
                raise CmbxParseError("RAW file contains too many points")
            terminated = False
            while block_offset < len(body):
                if body[block_offset] == 0x40:
                    block_offset += 1
                    terminated = True
                    break
                delta, block_offset = _read_signed_varint(body, block_offset)
                value = block_values[-1] + delta
                if not -(1 << 63) <= value <= _MAX_SIGNED_INTEGER:
                    raise CmbxParseError(
                        "RAW signal value exceeds the allowed range"
                    )
                block_values.append(value)
                if diff_point_count + len(block_values) > _MAX_SIGNAL_POINTS:
                    raise CmbxParseError(
                        "RAW file contains too many points"
                    )
            if not terminated or block_offset != len(body):
                raise CmbxParseError(
                    "PtsLDiff section is truncated or contains data after the terminator"
                )
            diff_blocks.append((start_time_min, block_values))
            diff_point_count += len(block_values)
        elif tag == _PTS_DDCMP_TAG:
            if signal_name is None:
                raise CmbxParseError("RAW file has point section data before SignHdr")
            if signal_encoding is not None and signal_encoding != tag:
                raise CmbxParseError("RAW file mixes point section formats")
            signal_encoding = tag
            ddcmp_points.extend(
                _decode_ddcmp_block(body, existing_point_count=len(ddcmp_points))
            )
        else:
            section_name = tag.decode("ascii", errors="replace").rstrip("\0")
            raise CmbxParseError(
                f"Unsupported RAW section {section_name!r}; "
                f"for {channel_name} the supported formats are "
                "PtsLL2Df, PtsLDiff, and PtsDDCmp"
            )

        offset += section_size

    if signal_name is None or signal_encoding is None:
        raise CmbxParseError(
            f"RAW file is missing the header or points for {channel_name}"
        )
    if signal_encoding == _PTS_LDIFF_TAG:
        if not diff_blocks:
            raise CmbxParseError(f"RAW file is missing points for {channel_name}")
        if len(diff_blocks) == 1 and len(diff_blocks[0][1]) > 1:
            raise CmbxParseError(
                "Unable to determine the time interval for a single PtsLDiff block"
            )

        # Each PtsLDiff block starts after the samples in the preceding block;
        # use adjacent starts to recover the common sample interval.
        interval_candidates = [
            (next_start - start) * 60 / len(block_values)
            for (start, block_values), (next_start, _) in zip(
                diff_blocks, diff_blocks[1:]
            )
        ]
        if interval_candidates:
            step_s = interval_candidates[0]
            if not math.isfinite(step_s) or step_s <= 0 or any(
                not math.isclose(step_s, candidate, rel_tol=1e-9, abs_tol=1e-12)
                for candidate in interval_candidates[1:]
            ):
                raise CmbxParseError(
                    "PtsLDiff block time intervals are incompatible"
                )
        else:
            step_s = 0.0

        times_min = tuple(
            start_time_min + point_index * step_s / 60
            for start_time_min, block_values in diff_blocks
            for point_index in range(len(block_values))
        )
        values = tuple(
            value * _PTS_LDIFF_VALUE_SCALE
            for _, block_values in diff_blocks
            for value in block_values
        )
        steps_s = (None,) + (step_s,) * (len(values) - 1)
    elif signal_encoding == _PTS_DDCMP_TAG:
        # PtsDDCmp carries each time/value pair explicitly, so its sampling
        # interval may vary between adjacent points.
        times_min = tuple(time_min for time_min, _ in ddcmp_points)
        values = tuple(value for _, value in ddcmp_points)
        steps_s = (None,) + tuple(
            (current - previous) * 60
            for previous, current in zip(times_min, times_min[1:])
        )
    else:
        if time_units is None or value_units is None or not ticks:
            raise CmbxParseError(
                f"RAW file is missing the header or points for {channel_name}"
            )
        if any(current <= previous for previous, current in zip(ticks, ticks[1:])):
            raise CmbxParseError(
                f"Timestamps for {channel_name} in the RAW file are not increasing"
            )
        times_min = tuple(tick / time_units for tick in ticks)
        values = tuple(value / value_units for value in values_raw)
        steps_s = (None,) + tuple(
            (current - previous) / time_units * 60
            for previous, current in zip(ticks, ticks[1:])
        )

    if not all(math.isfinite(value) for value in (*times_min, *values)):
        raise CmbxParseError(
            "RAW file contains non-numeric values outside the valid range"
        )
    if any(current <= previous for previous, current in zip(times_min, times_min[1:])):
        raise CmbxParseError(f"Timestamps for {channel_name} in the RAW file are not increasing")
    return Chromatogram(
        source=source,
        metadata={
            "Name": injection.name,
            "Channel": channel_name,
            "Signal Unit": signal_unit,
            "Data Points": str(len(values)),
            "URL": injection.url,
        },
        times_min=times_min,
        values=values,
        steps_s=steps_s,
    )


def read_cmbx_channel(
    sequence: CmbxSequence,
    injection: CmbxInjection,
    channel_name: str,
    *,
    signal_unit: str | None = None,
) -> Chromatogram:
    """Read a channel stored using a supported RAW point encoding.

    Units are known for UV_VIS_1, Emission_1, Pump_Pressure, and CC_Temp.
    Supply ``signal_unit`` for other channels when their physical unit is known
    from external metadata. PtsLL2Df, PtsLDiff, and PtsDDCmp point encodings
    are supported.

    Raises:
        CmbxParseError: If the channel/member is missing or its RAW data is
            malformed, too large, or uses an unsupported encoding.
    """
    channel = next(
        (item for item in injection.channels if item.name == channel_name),
        None,
    )
    if channel is None:
        raise CmbxParseError(
            f"Injection {injection.name!r} does not contain channel {channel_name}"
        )
    if not channel.raw_filename:
        raise CmbxParseError(
            f"Channel {channel_name} in injection {injection.name!r} has no RAW file"
        )

    try:
        with ZipFile(sequence.source) as archive:
            try:
                info = archive.getinfo(channel.raw_filename)
            except KeyError as error:
                raise CmbxParseError(
                    f"CMBX archive is missing RAW file {channel.raw_filename}"
                ) from error
            if info.file_size > _MAX_RAW_BYTES:
                raise CmbxParseError(
                    f"RAW file {channel.raw_filename} is too large"
                )
            if channel.size is not None and info.file_size != channel.size:
                raise CmbxParseError(
                    f"RAW file size for {channel.raw_filename} "
                    "does not match header.xml"
                )
            raw_data = archive.read(info)
    except CmbxParseError:
        raise
    except (OSError, BadZipFile, RuntimeError) as error:
        raise CmbxParseError(
            f"Failed to read RAW file {channel.raw_filename}: {error}"
        ) from error

    resolved_unit = (
        _KNOWN_SIGNAL_UNITS.get(channel_name, "")
        if signal_unit is None
        else signal_unit
    )
    return _decode_signal(
        raw_data, injection, sequence.source, channel_name, resolved_unit
    )


def read_cmbx_uv_vis_1(
    sequence: CmbxSequence, injection: CmbxInjection
) -> Chromatogram:
    """Read the UV_VIS_1 signal in mAU from an injection's RAW member."""
    return read_cmbx_channel(sequence, injection, "UV_VIS_1")


def read_cmbx_emission_1(
    sequence: CmbxSequence, injection: CmbxInjection
) -> Chromatogram:
    """Read the Emission_1 signal in counts from an injection's RAW member."""
    return read_cmbx_channel(sequence, injection, "Emission_1")


def read_cmbx_pump_pressure(
    sequence: CmbxSequence, injection: CmbxInjection
) -> Chromatogram:
    """Read the Pump_Pressure signal in bar from an injection's RAW member."""
    return read_cmbx_channel(sequence, injection, "Pump_Pressure")


def read_cmbx_cc_temp(
    sequence: CmbxSequence, injection: CmbxInjection
) -> Chromatogram:
    """Read the CC_Temp signal in degrees Celsius from an injection's RAW member."""
    return read_cmbx_channel(sequence, injection, "CC_Temp")


def read_cmbx(path: Path | str) -> CmbxSequence:
    """Read sequence, injection, and channel metadata from a CMBX ZIP file.

    Signal samples are not decoded here; use :func:`read_cmbx_channel` to load
    a channel when needed.

    Raises:
        CmbxParseError: If the archive or its ``header.xml`` is invalid, or the
            sequence contains no injections.
    """
    source = Path(path)
    try:
        with ZipFile(source) as archive:
            entries = frozenset(archive.namelist())
            try:
                header = ElementTree.fromstring(archive.read("header.xml"))
            except KeyError as error:
                raise CmbxParseError("Archive is missing header.xml") from error
            except ElementTree.ParseError as error:
                raise CmbxParseError("Failed to parse header.xml") from error
    except (OSError, BadZipFile) as error:
        raise CmbxParseError(
            f"Failed to open CMBX archive {source}: {error}"
        ) from error

    sequence_element = header.find("./ChromeleonElement")
    if sequence_element is None:
        raise CmbxParseError("header.xml is missing the sequence object")

    # Injections are direct children of the sequence; their children describe
    # channels and point to optional RAW members.
    injections: list[CmbxInjection] = []
    for element in sequence_element.findall("./ChromeleonElement"):
        if not element.attrib.get("ItemType", "").endswith(".Injection"):
            continue
        channels: list[CmbxChannel] = []
        for child in element.findall("./ChromeleonElement"):
            raw_size = child.attrib.get("Size")
            try:
                size = int(raw_size) if raw_size else None
            except ValueError as error:
                raise CmbxParseError(
                    "Invalid RAW file size for channel "
                    f"{child.attrib.get('Name', '')}"
                ) from error
            channels.append(
                CmbxChannel(
                    name=child.attrib.get("Name", "Unknown channel"),
                    item_type=child.attrib.get("ItemType", ""),
                    raw_filename=child.attrib.get("RawDataFilename"),
                    size=size,
                    url=child.attrib.get("Url", ""),
                )
            )
        injections.append(
            CmbxInjection(
                name=element.attrib.get("Name", "Unnamed injection"),
                url=element.attrib.get("Url", ""),
                channels=tuple(channels),
            )
        )

    if not injections:
        raise CmbxParseError("No injections were found in the sequence")

    return CmbxSequence(
        source=source,
        name=sequence_element.attrib.get("Name", source.stem),
        injections=tuple(injections),
        archive_entries=entries,
    )


__all__ = [
    "Chromatogram",
    "CmbxChannel",
    "CmbxInjection",
    "CmbxParseError",
    "CmbxSequence",
    "export_matches_injection",
    "read_cmbx",
    "read_cmbx_channel",
    "read_cmbx_cc_temp",
    "read_cmbx_emission_1",
    "read_cmbx_pump_pressure",
    "read_cmbx_uv_vis_1",
]
