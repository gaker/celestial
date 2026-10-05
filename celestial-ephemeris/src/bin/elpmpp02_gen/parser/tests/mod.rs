use super::records::*;
use super::*;
use celestial_core::constants::{HALF_PI, PI, TWOPI};
use std::io::ErrorKind;
use Coordinate::{Distance, Latitude, Longitude};

mod files;
mod records;

const MAIN_LINE: &str = "  0  2  0  0     -411.60287      168.48   -18433.81     -121.62        0.40       -0.18        0.00";
const PERT_LINE: &str =
    "    1-0.1274921554086D+02 0.6368794709728D+01  0  0  1  0  0-18 16  0  0  0  0  0  0  0  0  0";

fn with_field(line: &str, at: usize, text: &str) -> String {
    format!("{}{}{}", &line[..at], text, &line[at + text.len()..])
}

fn pert(sin_coeff: f64, cos_coeff: f64) -> PertTerm {
    PertTerm {
        sin_coeff,
        cos_coeff,
        multipliers: [0; 16],
    }
}

#[test]
fn coordinate_names() {
    assert_eq!(Longitude.to_string(), "Longitude");
    assert_eq!(Latitude.to_string(), "Latitude");
    assert_eq!(Distance.to_string(), "Distance");
}

#[test]
fn main_amplitude_is_the_size_of_the_first_coefficient() {
    for a0 in [22639.55, -12345.67] {
        let term = MainTerm {
            delaunay: [0, 0, 1, 0],
            coeffs: [a0, 1.0, 2.0, 412529.61, 3.0, 4.0, 5.0],
        };
        assert_eq!(term.amplitude(), a0.abs());
    }
}

#[test]
fn pert_amplitude_and_phase() {
    assert_eq!(pert(3.0, -4.0).amplitude(), 5.0);
    // Negative phases are stored in the first turn, as the authors' code does.
    let phases = [
        (pert(1.0, 0.0), 0.0),
        (pert(0.0, 1.0), HALF_PI),
        (pert(-1.0, 0.0), PI),
        (pert(0.0, -1.0), TWOPI - HALF_PI),
    ];
    for (term, phase) in phases {
        assert_eq!(term.phase(), phase);
    }
}

#[test]
fn totals() {
    let main = |coordinate, n| MainSeries {
        coordinate,
        terms: vec![parse_main_term(MAIN_LINE).unwrap(); n],
    };
    let pert_series = |coordinate, sizes: &[usize]| PertSeries {
        coordinate,
        blocks: (sizes.iter().enumerate())
            .map(|(power, &n)| PertBlock {
                time_power: power as u8,
                terms: vec![pert(1.0, 0.0); n],
            })
            .collect(),
    };
    let data = ElpData {
        main: [main(Longitude, 2), main(Latitude, 1), main(Distance, 3)],
        pert: [
            pert_series(Longitude, &[2]),
            pert_series(Latitude, &[1, 0]),
            pert_series(Distance, &[]),
        ],
    };
    assert_eq!(data.total_main_terms(), 6);
    assert_eq!(data.total_pert_terms(), 3);
    assert_eq!(data.total_terms(), 9);
}

#[test]
fn parse_error_messages() {
    let io = ParseError::from(std::io::Error::new(ErrorKind::PermissionDenied, "denied"));
    assert!(matches!(io, ParseError::IoError(ref e) if e.kind() == ErrorKind::PermissionDenied));
    let cases = [
        (io, "IO error: denied"),
        (ParseError::InvalidHeader("a".into()), "Invalid header: a"),
        (
            ParseError::InvalidMainTerm("b".into()),
            "Invalid main term: b",
        ),
        (
            ParseError::InvalidPertTerm("c".into()),
            "Invalid pert term: c",
        ),
        (ParseError::InvalidFormat("d".into()), "Invalid format: d"),
    ];
    for (error, message) in cases {
        assert_eq!(error.to_string(), message);
    }
    fn assert_std_error<T: std::error::Error>() {}
    assert_std_error::<ParseError>();
}
