use super::fixture::*;
use crate::jpl::bodies::{EARTH_MOON_BARYCENTER as EMB, MOON};
use celestial_core::constants::J2000_JD;

#[test]
fn summary_chain_that_loops_is_rejected() {
    let bytes = corrupted(|b| {
        let at = summary_record(b);
        put_f64(b, at, 3.0)
    });
    assert_eq!(
        rejected(&bytes),
        "summary record 3 is reached twice; the record chain loops"
    );
    let bytes = corrupted(|b| {
        let at = summary_record(b);
        put_f64(b, at, 3.0);
        put_f64(b, at + 16, 0.0);
    });
    rejected(&bytes);
}

#[test]
fn summary_count_must_be_a_whole_number_that_fits() {
    for nsum in [13.7, 26.0, -1.0, f64::NAN] {
        let bytes = corrupted(|b| {
            let at = summary_record(b) + 16;
            put_f64(b, at, nsum)
        });
        assert_eq!(
            rejected(&bytes),
            format!(
                "summary record 3 holds {} summaries; expected a whole number from 0 to 25",
                nsum
            )
        );
    }
}

#[test]
fn record_pointers_must_name_a_summary_record() {
    for fward in [-1, 0, 1, 99999] {
        let bytes = corrupted(|b| put_i32(b, 76, fward));
        assert_eq!(
            rejected(&bytes),
            format!(
                "first summary record {} is not a record of this file (2 to 10640)",
                fward
            )
        );
    }
    let bytes = corrupted(|b| {
        let at = summary_record(b);
        put_f64(b, at, 2.5)
    });
    assert_eq!(
        rejected(&bytes),
        "next summary record 2.5 is not a record of this file (2 to 10640)"
    );
}

#[test]
fn summary_format_must_be_an_spk() {
    for (nd, ni) in [(2, -1), (3, 6), (2, 5)] {
        let bytes = corrupted(|b| {
            put_i32(b, 8, nd);
            put_i32(b, 12, ni);
        });
        assert_eq!(
            rejected(&bytes),
            format!(
                "summary format ND = {}, NI = {}; an SPK has ND = 2, NI = 6",
                nd, ni
            )
        );
    }
}

#[test]
fn binary_format_must_be_ieee() {
    let bytes = corrupted(|b| b[88..96].copy_from_slice(b"VAX-GFLT"));
    assert_eq!(
        rejected(&bytes),
        "binary format \"VAX-GFLT\" is not LTL-IEEE or BIG-IEEE"
    );
    let bytes = corrupted(|b| b[88..96].copy_from_slice(b"        "));
    let pristine = load(&de432s()).unwrap();
    let blank = load(&bytes).unwrap();
    assert_eq!(
        state(&blank, MOON, EMB, J2000_JD, 0.37).unwrap(),
        state(&pristine, MOON, EMB, J2000_JD, 0.37).unwrap()
    );
}

#[test]
fn id_word_must_name_an_spk() {
    let bytes = corrupted(|b| b[0..8].copy_from_slice(b"DAF/CK  "));
    assert_eq!(
        rejected(&bytes),
        "ID word \"DAF/CK  \" is not DAF/SPK or NAIF/DAF"
    );
    let bytes = corrupted(|b| b[0..8].copy_from_slice(b"NAIF/DAF"));
    assert!(load(&bytes).is_ok());
}

#[test]
fn only_a_damaged_ftp_string_is_rejected() {
    let bytes = corrupted(|b| b[699..727].fill(0));
    assert!(load(&bytes).is_ok());
    let bytes = corrupted(|b| b[706] = b'\n');
    assert_eq!(
        rejected(&bytes),
        "FTP validation string is damaged; the file was probably transferred in text mode"
    );
}

#[test]
fn truncated_file_is_rejected_at_open() {
    let bytes = de432s();
    assert_eq!(
        rejected(&bytes[..8192]),
        "segment 0 (body 1, center 0): data addresses 513 to 201508 are not \
         an ordered range inside 1 to 1024"
    );
    assert_eq!(
        rejected(&bytes[..1000]),
        "the file is shorter than one 1024-byte DAF record"
    );
}
