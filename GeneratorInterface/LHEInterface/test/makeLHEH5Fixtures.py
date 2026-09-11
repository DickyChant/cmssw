"""Small independent LHEH5 2.0.0 fixtures; h5py is a test-only dependency.

Column conventions cross-checked against pylhe f7889ce52b1d93dd0262656a624b53b20885d5d3.
Generate in a caller-owned temporary directory; never modify checked-in data.
"""
import argparse
from pathlib import Path

import h5py
import numpy as np

INIT = "beamA beamB energyA energyB PDFgroupA PDFgroupB PDFsetA PDFsetB weightingStrategy numProcesses".split()
PROC = "procId npLO npNLO xSection error unitWeight".split()
EVENT = "pid nparticles start trials scale fscale rscale aqed aqcd NOMINAL event_num".split()
PARTICLE = "id status mother1 mother2 color1 color2 px py pz e m lifetime spin".split()
PARTICLES = np.array([
    [11, -1, 0, 0, 0, 0, 0, 0, 45.5, 45.5, 0, 0, 9],
    [-11, -1, 0, 0, 0, 0, 0, 0, -45.5, 45.5, 0, 0, 9],
    [23, 2, 1, 2, 0, 0, 0, 0, 0, 91, 91, 0, 9],
    [13, 1, 3, 3, 0, 0, 45.5, 0, 0, 45.5, 0, 0, 9],
    [-13, 1, 3, 3, 0, 0, -45.5, 0, 0, 45.5, 0, 0, 9],
], dtype="f8")


def fixture(directory, name, count, *, permute=False, strategy=-4, run_energy=45.5, number="valid", xml_weights=False):
    with h5py.File(directory / f"{name}.h5", "w") as file, (directory / f"{name}.lhe").open("w") as xml:
        file.create_dataset("version", data=[2, 0, 0], dtype="i8")
        init = [11, -11, run_energy, run_energy, 0, 0, 0, 0, strategy, 1]

        def table(label, data, names, **kwargs):
            if permute:
                order = list(reversed(range(len(names))))
                data = np.asarray(data)[..., order]
                names = [names[i] for i in order]
            dataset = file.create_dataset(label, data=data, dtype="f8", **kwargs)
            attr = label if permute else "properties"
            dataset.attrs[attr] = np.asarray(names, dtype=h5py.string_dtype())
            return dataset

        table("init", init, INIT)
        table("procInfo", [[1, np.nan, np.nan, 1, 0.01, 2]], PROC)
        # Keep generated fixture memory bounded too; append small batches below.
        ev = table("events", np.empty((0, len(EVENT))), EVENT, maxshape=(None, len(EVENT)), chunks=(64, len(EVENT)))
        pt = table("particles", np.empty((0, len(PARTICLE))), PARTICLE, maxshape=(None, len(PARTICLE)), chunks=(320, len(PARTICLE)))
        xml.write('<LesHouchesEvents version="3.0">\n<init>\n')
        xml.write(" ".join(map(str, init)) + "\n1 0.01 2 1\n</init>\n")
        for first in range(0, count, 64):
            stop = min(first + 64, count)
            rows = []
            for i in range(first, stop):
                weight = -2.0 if i % 2 else 1.5
                if strategy > 0:
                    weight = abs(weight)
                rows.append([1, 5, i * 5, np.nan, 91, np.nan, np.nan, 0.007297, 0.118, weight, i + 100])
                xml.write(f"<event>\n5 1 {weight} 91 0.007297 0.118\n")
                for particle in PARTICLES:
                    xml.write(" ".join(format(x, ".17g") for x in particle) + "\n")
                if xml_weights:
                    xml.write('<rwgt><wgt id="scale_up">2</wgt><wgt id="scale_down">1</wgt></rwgt>\n')
                if number == "missing":
                    xml.write('<event_num/>\n')
                elif number == "invalid":
                    xml.write('<event_num num="9999999999999999999999999"/>\n')
                elif number != "none":
                    xml.write(f'<event_num num="{i + 100}"/>\n')
                xml.write('</event>\n')
            particles = np.tile(PARTICLES, (stop - first, 1))
            rows = np.asarray(rows)
            if permute:
                rows = rows[:, ::-1]
                particles = particles[:, ::-1]
            ev.resize((stop, len(EVENT)))
            ev[first:stop] = rows
            pt.resize((stop * 5, len(PARTICLE)))
            pt[first * 5:stop * 5] = particles
        xml.write("</LesHouchesEvents>\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--large", type=int, default=0)
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    if args.large:
        fixture(args.directory, "large", args.large)
        return
    for count in (0, 1, 3, 1599, 1600, 1601):
        fixture(args.directory, f"events{count}", count)
    fixture(args.directory, "permuted", 3, permute=True)
    fixture(args.directory, "signed_unit", 3, strategy=-3)
    fixture(args.directory, "different_run", 3, run_energy=50)
    fixture(args.directory, "xml_weights", 3, xml_weights=True)
    fixture(args.directory, "missing_num", 1, number="missing")
    fixture(args.directory, "invalid_num", 1, number="invalid")
    for name in ("bad_version", "bad_offset", "bad_count", "bad_mother", "duplicate_label",
                 "missing_label", "nonfinite", "unknown_process", "metadata", "counterterms", "negative_positive_strategy", "large_chunk"):
        fixture(args.directory, name, 1)
        with h5py.File(args.directory / f"{name}.h5", "r+") as file:
            if name == "bad_version":
                file["version"][0] = 99
            elif name == "bad_offset":
                file["events"][0, EVENT.index("start")] = 0.5
            elif name == "bad_count":
                file["events"][0, EVENT.index("nparticles")] = 100001
            elif name == "bad_mother":
                file["particles"][3, PARTICLE.index("mother1")] = 99
            elif name == "duplicate_label":
                labels = EVENT.copy()
                labels[-1] = "pid"
                file["events"].attrs["properties"] = np.asarray(labels, dtype=h5py.string_dtype())
            elif name == "missing_label":
                del file["events"].attrs["properties"]
            elif name == "nonfinite":
                file["events"][0, EVENT.index("NOMINAL")] = np.inf
            elif name == "unknown_process":
                file["events"][0, EVENT.index("pid")] = 999
            elif name == "metadata":
                file["events"][0, EVENT.index("rscale")] = 45.5
            elif name == "counterterms":
                file.create_dataset("ctevents", data=[[1.]])
            elif name == "negative_positive_strategy":
                file["init"][INIT.index("weightingStrategy")] = 4
                file["events"][0, EVENT.index("NOMINAL")] = -1
            elif name == "large_chunk":
                del file["particles"]
                dataset = file.create_dataset("particles", shape=(0, 13), maxshape=(None, 13), chunks=(1000000, 13), dtype="f8")
                dataset.attrs["properties"] = np.asarray(PARTICLE, dtype=h5py.string_dtype())


if __name__ == "__main__":
    main()
