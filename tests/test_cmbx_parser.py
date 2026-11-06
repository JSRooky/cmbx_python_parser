import struct
import tempfile
import unittest
from pathlib import Path
from zipfile import ZipFile

from cmbx_parser import CmbxParseError, read_cmbx, read_cmbx_channel


def _encode_signed_varint(value: int) -> bytes:
    negative = value < 0
    magnitude = abs(value)
    first = magnitude & 0x3F
    magnitude >>= 6
    if negative:
        first |= 0x40
    if magnitude:
        first |= 0x80
    result = bytearray((first,))
    while magnitude:
        byte = magnitude & 0x7F
        magnitude >>= 7
        if magnitude:
            byte |= 0x80
        result.append(byte)
    return bytes(result)


def _raw_section(tag: bytes, body: bytes) -> bytes:
    return tag + struct.pack("<Q", 16 + len(body)) + body


def _signal_archive(path: Path) -> None:
    channel_name = b"Conductivity"
    signal_header = (
        b"\x01\x02\x07"
        + b"\x02\x0a"
        + struct.pack("<Q", 134_000_000_000_000_000)
        + b"\x03"
        + bytes((len(channel_name),))
        + channel_name
        + bytes(16)
    )
    units = struct.pack("<dd", 60_000_000_000.0, 10.0)
    encoded_points = b"".join(
        _encode_signed_varint(value)
        for value in (0, 20, 60_000_000_000, 10, 0)
    )
    raw_data = (
        _raw_section(b"SignHdr\0", signal_header)
        + _raw_section(b"PtsLL2Df", units + encoded_points)
    )
    header = (
        '<Chromeleon><ChromeleonElement Name="test-sequence">'
        '<ChromeleonElement Name="sample-1" Url="test-sequence/1.smp" '
        'ItemType="Dionex.Chromeleon.Data.Injection">'
        '<ChromeleonElement Name="Conductivity" '
        'Url="test-sequence/1.smp/Conductivity.channel" '
        'ItemType="Dionex.Chromeleon.Data.Signal" '
        'RawDataFilename="1_signal.raw" '
        f'Size="{len(raw_data)}"/>'
        "</ChromeleonElement></ChromeleonElement></Chromeleon>"
    )
    with ZipFile(path, "w") as archive:
        archive.writestr("header.xml", header)
        archive.writestr("1_signal.raw", raw_data)


class CmbxParserTests(unittest.TestCase):
    def test_reads_sequence_metadata_and_signal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "sequence.cmbx"
            _signal_archive(path)

            sequence = read_cmbx(path)
            injection = sequence.injections[0]
            chromatogram = read_cmbx_channel(
                sequence, injection, "Conductivity", signal_unit="uS/cm"
            )

        self.assertEqual(sequence.name, "test-sequence")
        self.assertEqual(injection.name, "sample-1")
        self.assertEqual(injection.signal_channels, ("Conductivity",))
        self.assertEqual(sequence.signal_channels, ("Conductivity",))
        self.assertEqual(chromatogram.times_min, (0.0, 1.0))
        self.assertEqual(chromatogram.values, (2.0, 3.0))
        self.assertEqual(chromatogram.steps_s, (None, 60.0))
        self.assertEqual(chromatogram.signal_unit, "uS/cm")

    def test_rejects_non_zip_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "invalid.cmbx"
            path.write_text("not a zip archive")

            with self.assertRaises(CmbxParseError):
                read_cmbx(path)

    def test_rejects_missing_channel(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "sequence.cmbx"
            _signal_archive(path)
            sequence = read_cmbx(path)

            with self.assertRaises(CmbxParseError):
                read_cmbx_channel(
                    sequence, sequence.injections[0], "Missing"
                )


if __name__ == "__main__":
    unittest.main()
