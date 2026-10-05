mod api;
mod authors;

use celestial_time::julian::JulianDate;
use celestial_time::scales::tdb::TDB;

fn tdb(jd: f64) -> TDB {
    TDB::from_julian_date(JulianDate::new(jd, 0.0))
}
