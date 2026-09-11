"""Compare actual persisted LHE products, including the event-number regression."""
import argparse

from DataFormats.FWLite import Events, Handle, Runs


def event_records(path, shower, expect_loss=False):
    events = Events(path)
    handle = Handle("io_v1::LHEEventProduct")
    generator = Handle("io_v1::GenEventInfoProduct")
    for event in events:
        event.getByLabel("source", handle)
        assert handle.isValid(), "Missing persisted LHEEventProduct"
        product = handle.product()
        record = product.hepeup()
        if expect_loss:
            assert any("LHEH5 core-only input; omitted:" in str(product.getComment(i))
                       for i in range(product.comments_size())), "Missing persisted loss declaration"
        if shower:
            event.getByLabel("generator", generator)
            assert generator.isValid(), "Missing shower output"
        assert product.evtnum() >= 100, "Original event number was not persisted"
        yield (
            (event.eventAuxiliary().id().run(), event.eventAuxiliary().id().event()),
            record.NUP, record.IDPRUP, record.XWGTUP, product.originalXWGTUP(),
            record.SCALUP, record.AQEDUP, record.AQCDUP,
            tuple(record.IDUP), tuple(record.ISTUP),
            tuple((x.first, x.second) for x in record.MOTHUP),
            tuple((x.first, x.second) for x in record.ICOLUP),
            tuple(tuple(record.PUP[i][j] for j in range(5)) for i in range(record.NUP)),
            tuple(record.VTIMUP), tuple(record.SPINUP),
            product.npLO(), product.npNLO(), product.evtnum(), tuple(product.scales()),
            tuple((w.id, w.wgt) for w in product.weights()),
        )


def run_records(path):
    handle = Handle("io_v1::LHERunInfoProduct")
    for run in Runs(path):
        run.getByLabel("source", handle)
        assert handle.isValid(), "Missing persisted LHERunInfoProduct"
        record = handle.product().heprup()
        yield (
            (record.IDBMUP.first, record.IDBMUP.second), (record.EBMUP.first, record.EBMUP.second),
            (record.PDFGUP.first, record.PDFGUP.second), (record.PDFSUP.first, record.PDFSUP.second),
            record.IDWTUP, record.NPRUP, tuple(record.XSECUP), tuple(record.XERRUP),
            tuple(record.XMAXUP), tuple(record.LPRUP),
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xml")
    parser.add_argument("hdf5")
    parser.add_argument("count", type=int)
    parser.add_argument("--shower", action="store_true")
    parser.add_argument("--loss", action="store_true")
    args = parser.parse_args()
    from itertools import zip_longest
    count = 0
    for a, b in zip_longest(event_records(args.xml, args.shower), event_records(args.hdf5, args.shower, args.loss)):
        assert a is not None and a == b, f"Persisted event mismatch at {count}"
        count += 1
    assert count == args.count, (count, args.count)
    assert list(run_records(args.xml)) == list(run_records(args.hdf5)), "Run metadata mismatch"
    print(f"Persisted XML/HDF5 parity PASS: {count} events; shower={args.shower}")


if __name__ == "__main__":
    main()
