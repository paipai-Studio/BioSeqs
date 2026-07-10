#!/usr/bin/env python3
import json
from Bio.Seq import Seq
from Bio.SeqFeature import (
    SeqFeature, FeatureLocation, CompoundLocation,
    ExactPosition, WithinPosition, BeforePosition, AfterPosition,
    UncertainPosition, UnknownPosition
)

def json_serial(obj):
    if isinstance(obj, Seq):
        return str(obj)
    if isinstance(obj, SeqFeature):
        return {
            "type": obj.type,
            "id": obj.id,
            "location": json_serial(obj.location) if obj.location else None,
            "strand": obj.strand,
            "qualifiers": {k: v for k, v in obj.qualifiers.items()}
        }
    if isinstance(obj, FeatureLocation):
        return {
            "start": json_serial(obj.start),
            "end": json_serial(obj.end),
            "strand": obj.strand,
            "ref": obj.ref,
            "ref_db": obj.ref_db,
            "start_int": int(obj.start) if obj.start else None,
            "end_int": int(obj.end) if obj.end else None
        }
    if isinstance(obj, CompoundLocation):
        return {
            "operator": obj.operator,
            "parts": [json_serial(p) for p in obj.parts]
        }
    if isinstance(obj, (ExactPosition, WithinPosition, BeforePosition,
                        AfterPosition, UncertainPosition, UnknownPosition)):
        return {
            "type": type(obj).__name__,
            "position": int(obj) if obj.position is not None else None,
            "fuzzy_start": obj.fuzzy_start if hasattr(obj, 'fuzzy_start') else None,
            "fuzzy_end": obj.fuzzy_end if hasattr(obj, 'fuzzy_end') else None
        }
    if isinstance(obj, type(None)):
        return None
    return str(obj)

def test_feature_location_basic():
    loc = FeatureLocation(ExactPosition(5), ExactPosition(10))
    return {
        "start": int(loc.start),
        "end": int(loc.end),
        "strand": loc.strand,
        "str": str(loc),
        "repr": repr(loc)
    }

def test_feature_location_strand():
    loc_plus = FeatureLocation(7, 110, strand=1)
    loc_minus = FeatureLocation(9, 108, strand=-1)
    return {
        "plus_strand": {
            "start": int(loc_plus.start),
            "end": int(loc_plus.end),
            "strand": loc_plus.strand,
            "str": str(loc_plus)
        },
        "minus_strand": {
            "start": int(loc_minus.start),
            "end": int(loc_minus.end),
            "strand": loc_minus.strand,
            "str": str(loc_minus)
        }
    }

def test_seqfeature_basic():
    f1 = SeqFeature(FeatureLocation(5, 10), type="domain")
    f2 = SeqFeature(FeatureLocation(7, 110, strand=1), type="CDS")
    f3 = SeqFeature(FeatureLocation(9, 108, strand=-1), type="CDS")
    return {
        "f1_type": f1.type,
        "f1_str": str(f1),
        "f2_type": f2.type,
        "f2_str": str(f2),
        "f3_type": f3.type,
        "f3_str": str(f3),
        "f1_location_strand": f1.location.strand if f1.location else None,
        "f2_location_strand": f2.location.strand if f2.location else None,
        "f3_location_strand": f3.location.strand if f3.location else None
    }

def test_seqfeature_extract():
    seq = Seq("MKQHKAMIVALIVICITAVVAAL")
    f = SeqFeature(FeatureLocation(8, 15), type="domain")
    extracted = f.extract(seq)
    return {
        "parent_seq": str(seq),
        "extracted_seq": str(extracted),
        "extracted_length": len(extracted)
    }

def test_seqfeature_extract_minus_strand():
    seq = Seq("ATGCCGATCGATCGATCGAT")
    f = SeqFeature(FeatureLocation(0, 6, strand=-1), type="CDS")
    extracted = f.extract(seq)
    return {
        "parent_seq": str(seq),
        "extracted_seq": str(extracted),
        "strand": f.location.strand if f.location else None
    }

def test_seqfeature_qualifiers():
    f = SeqFeature(FeatureLocation(0, 10), type="gene")
    f.qualifiers["gene"] = ["abc"]
    f.qualifiers["note"] = ["interesting gene"]
    return {
        "qualifiers": {k: v for k, v in f.qualifiers.items()},
        "has_gene": "gene" in f.qualifiers,
        "gene_value": f.qualifiers.get("gene", [None])[0]
    }

def test_compound_location():
    loc1 = FeatureLocation(0, 10, strand=1)
    loc2 = FeatureLocation(20, 30, strand=1)
    compound = CompoundLocation([loc1, loc2], operator="join")
    return {
        "operator": compound.operator,
        "num_parts": len(compound.parts),
        "str": str(compound),
        "repr": repr(compound)
    }

def test_position_types():
    exact = ExactPosition(10)
    within = WithinPosition(5, 5, 15)
    before = BeforePosition(10)
    after = AfterPosition(10)
    uncertain = UncertainPosition(10)
    unknown = UnknownPosition()
    return {
        "exact": {"type": type(exact).__name__, "position": int(exact)},
        "within": {"type": type(within).__name__, "position": int(within)},
        "before": {"type": type(before).__name__, "position": int(before)},
        "after": {"type": type(after).__name__, "position": int(after)},
        "uncertain": {"type": type(uncertain).__name__, "position": int(uncertain)},
        "unknown": {"type": type(unknown).__name__}
    }

def test_seqfeature_none_location():
    f = SeqFeature(None, type="domain")
    seq = Seq("MKQHKAMIVALIVICITAVVAAL")
    try:
        extracted = f.extract(seq)
        return {"error": None, "extracted": str(extracted)}
    except ValueError as e:
        return {"error": "ValueError", "message": str(e)}

def run_test(name, func):
    try:
        result = func()
        print(f"=== CASE: {name} ===")
        print(json.dumps(result, indent=2, default=json_serial))
        print("=== END ===")
    except Exception as e:
        print(f"=== CASE: {name} ===")
        print(json.dumps({"error": type(e).__name__, "message": str(e)}))
        print("=== END ===")

if __name__ == "__main__":
    run_test("feature_location_basic", test_feature_location_basic)
    run_test("feature_location_strand", test_feature_location_strand)
    run_test("seqfeature_basic", test_seqfeature_basic)
    run_test("seqfeature_extract", test_seqfeature_extract)
    run_test("seqfeature_extract_minus_strand", test_seqfeature_extract_minus_strand)
    run_test("seqfeature_qualifiers", test_seqfeature_qualifiers)
    run_test("compound_location", test_compound_location)
    run_test("position_types", test_position_types)
    run_test("seqfeature_none_location", test_seqfeature_none_location)