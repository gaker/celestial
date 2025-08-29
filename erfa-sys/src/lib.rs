use std::ffi::c_int;

extern "C" {
    pub fn eraDat(iy: c_int, im: c_int, id: c_int, fd: f64, deltat: *mut f64) -> c_int;

    pub fn eraUtctai(utc1: f64, utc2: f64, tai1: *mut f64, tai2: *mut f64) -> c_int;

    pub fn eraTaiutc(tai1: f64, tai2: f64, utc1: *mut f64, utc2: *mut f64) -> c_int;

    pub fn eraCal2jd(iy: c_int, im: c_int, id: c_int, djm0: *mut f64, djm: *mut f64) -> c_int;

    pub fn eraJd2cal(
        dj1: f64,
        dj2: f64,
        iy: *mut c_int,
        im: *mut c_int,
        id: *mut c_int,
        fd: *mut f64,
    ) -> c_int;

    pub fn eraTttdb(tt1: f64, tt2: f64, ut: f64, tdb1: *mut f64, tdb2: *mut f64) -> c_int;

    pub fn eraTdbtt(tdb1: f64, tdb2: f64, ut: f64, tt1: *mut f64, tt2: *mut f64) -> c_int;

    pub fn eraDtdb(date1: f64, date2: f64, ut: f64, elong: f64, u: f64, v: f64) -> f64;

    pub fn eraTaitt(tai1: f64, tai2: f64, tt1: *mut f64, tt2: *mut f64) -> c_int;

    pub fn eraTttai(tt1: f64, tt2: f64, tai1: *mut f64, tai2: *mut f64) -> c_int;

    // Geodetic/Geocentric coordinate conversions
    pub fn eraGd2gc(n: c_int, elong: f64, phi: f64, height: f64, xyz: *mut f64) -> c_int;
    pub fn eraGc2gd(n: c_int, xyz: *const f64, elong: *mut f64, phi: *mut f64, height: *mut f64) -> c_int;
    pub fn eraGd2gce(a: f64, f: f64, elong: f64, phi: f64, height: f64, xyz: *mut f64) -> c_int;
    pub fn eraGc2gde(a: f64, f: f64, xyz: *const f64, elong: *mut f64, phi: *mut f64, height: *mut f64) -> c_int;

    pub fn eraD2dtf(
        scale: *const i8,
        ndp: c_int,
        d1: f64,
        d2: f64,
        iy: *mut c_int,
        im: *mut c_int,
        id: *mut c_int,
        ihmsf: *mut [c_int; 4],
    ) -> c_int;

    pub fn eraDtf2d(
        scale: *const i8,
        iy: c_int,
        im: c_int,
        id: c_int,
        ihr: c_int,
        imn: c_int,
        sec: f64,
        d1: *mut f64,
        d2: *mut f64,
    ) -> c_int;

    pub fn eraGmst00(uta: f64, utb: f64, tta: f64, ttb: f64) -> f64;

    pub fn eraGmst06(uta: f64, utb: f64, tta: f64, ttb: f64) -> f64;

    pub fn eraGmst82(dj1: f64, dj2: f64) -> f64;

    pub fn eraGst00a(uta: f64, utb: f64, tta: f64, ttb: f64) -> f64;

    pub fn eraGst00b(uta: f64, utb: f64) -> f64;

    pub fn eraGst06(uta: f64, utb: f64, tta: f64, ttb: f64, rnpb: *const [f64; 9]) -> f64;

    pub fn eraGst06a(uta: f64, utb: f64, tta: f64, ttb: f64) -> f64;

    pub fn eraGst94(uta: f64, utb: f64) -> f64;

    pub fn eraEe00(date1: f64, date2: f64, epsa: f64, dpsi: f64) -> f64;

    pub fn eraEe00a(date1: f64, date2: f64) -> f64;

    pub fn eraEe00b(date1: f64, date2: f64) -> f64;

    pub fn eraEe06a(date1: f64, date2: f64) -> f64;

    pub fn eraEect00(date1: f64, date2: f64) -> f64;

    pub fn eraEqeq94(date1: f64, date2: f64) -> f64;

    pub fn eraEra00(dj1: f64, dj2: f64) -> f64;

    pub fn eraFad03(t: f64) -> f64;

    pub fn eraFae03(t: f64) -> f64;

    pub fn eraFaf03(t: f64) -> f64;

    pub fn eraFaju03(t: f64) -> f64;

    pub fn eraFal03(t: f64) -> f64;

    pub fn eraFalp03(t: f64) -> f64;

    pub fn eraFama03(t: f64) -> f64;

    pub fn eraFame03(t: f64) -> f64;

    pub fn eraFane03(t: f64) -> f64;

    pub fn eraFaom03(t: f64) -> f64;

    pub fn eraFapa03(t: f64) -> f64;

    pub fn eraFasa03(t: f64) -> f64;

    pub fn eraFaur03(t: f64) -> f64;

    pub fn eraFave03(t: f64) -> f64;

    pub fn eraPn00(
        date1: f64,
        date2: f64,
        dpsi: f64,
        deps: f64,
        epsa: *mut f64,
        rb: *mut [f64; 9],
        rp: *mut [f64; 9],
        rbp: *mut [f64; 9],
        rn: *mut [f64; 9],
        rbpn: *mut [f64; 9],
    );

    pub fn eraPn00a(
        date1: f64,
        date2: f64,
        dpsi: *mut f64,
        deps: *mut f64,
        epsa: *mut f64,
        rb: *mut [f64; 9],
        rp: *mut [f64; 9],
        rbp: *mut [f64; 9],
        rn: *mut [f64; 9],
        rbpn: *mut [f64; 9],
    );

    pub fn eraPn00b(
        date1: f64,
        date2: f64,
        dpsi: *mut f64,
        deps: *mut f64,
        epsa: *mut f64,
        rb: *mut [f64; 9],
        rp: *mut [f64; 9],
        rbp: *mut [f64; 9],
        rn: *mut [f64; 9],
        rbpn: *mut [f64; 9],
    );

    pub fn eraPn06(
        date1: f64,
        date2: f64,
        dpsi: f64,
        deps: f64,
        epsa: *mut f64,
        rb: *mut [f64; 9],
        rp: *mut [f64; 9],
        rbp: *mut [f64; 9],
        rn: *mut [f64; 9],
        rbpn: *mut [f64; 9],
    );

    pub fn eraPn06a(
        date1: f64,
        date2: f64,
        dpsi: *mut f64,
        deps: *mut f64,
        epsa: *mut f64,
        rb: *mut [f64; 9],
        rp: *mut [f64; 9],
        rbp: *mut [f64; 9],
        rn: *mut [f64; 9],
        rbpn: *mut [f64; 9],
    );

    pub fn eraPnm00a(date1: f64, date2: f64, rbpn: *mut [f64; 9]);

    pub fn eraPnm00b(date1: f64, date2: f64, rbpn: *mut [f64; 9]);

    pub fn eraPnm06a(date1: f64, date2: f64, rbpn: *mut [f64; 9]);

    pub fn eraPnm80(date1: f64, date2: f64, rbpn: *mut [f64; 9]);

    pub fn eraPom00(xp: f64, yp: f64, sp: f64, rpom: *mut [f64; 9]);

    pub fn eraS00(date1: f64, date2: f64, x: f64, y: f64) -> f64;

    pub fn eraS00a(date1: f64, date2: f64) -> f64;

    pub fn eraS00b(date1: f64, date2: f64) -> f64;

    pub fn eraS06(date1: f64, date2: f64, x: f64, y: f64) -> f64;

    pub fn eraS06a(date1: f64, date2: f64) -> f64;

    pub fn eraSp00(date1: f64, date2: f64) -> f64;

    pub fn eraXy06(date1: f64, date2: f64, x: *mut f64, y: *mut f64);

    pub fn eraXys00a(date1: f64, date2: f64, x: *mut f64, y: *mut f64, s: *mut f64);

    pub fn eraXys00b(date1: f64, date2: f64, x: *mut f64, y: *mut f64, s: *mut f64);

    pub fn eraXys06a(date1: f64, date2: f64, x: *mut f64, y: *mut f64, s: *mut f64);
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_erfa_cal2jd() {
        let mut djm0 = 0.0;
        let mut djm = 0.0;

        let status = unsafe { eraCal2jd(2000, 1, 1, &mut djm0, &mut djm) };

        assert_eq!(status, 0);
        assert_abs_diff_eq!(djm0, 2400000.5, epsilon = 1e-10);
        assert_abs_diff_eq!(djm, 51544.0, epsilon = 1e-10);
    }

    #[test]
    fn test_erfa_tai_tt() {
        let mut tt1 = 0.0;
        let mut tt2 = 0.0;

        let status = unsafe { eraTaitt(2400000.5, 51544.0, &mut tt1, &mut tt2) };

        assert_eq!(status, 0);
        assert_abs_diff_eq!(tt1, 2400000.5, epsilon = 1e-15);
        // TAI-TT offset is 32.184 seconds = 32.184/86400 days ≈ 0.0003725
        assert_abs_diff_eq!(tt2, 51544.0003725, epsilon = 1e-7);
    }
}
