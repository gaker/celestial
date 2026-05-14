use super::DetectedStar;

/// Drops detections that are within `min_spacing` pixels of a brighter detection.
///
/// Assumes the input is sorted by descending flux (as produced by
/// [`crate::detect::find_bright_stars`]). The brighter of two close detections is
/// kept; the fainter is discarded. Runs in place.
///
/// # Examples
///
/// ```rust,ignore
/// use celestial_solver::detect::{deduplicate, DetectedStar};
///
/// let mut stars: Vec<DetectedStar> = vec![/* sorted by descending flux */];
/// deduplicate(&mut stars, 20.0);
/// ```
pub fn deduplicate(stars: &mut Vec<DetectedStar>, min_spacing: f64) {
    let min_sq = min_spacing * min_spacing;
    let mut keep = vec![true; stars.len()];

    for i in 0..stars.len() {
        if !keep[i] {
            continue;
        }
        for j in (i + 1)..stars.len() {
            if !keep[j] {
                continue;
            }
            let dx = stars[j].x - stars[i].x;
            let dy = stars[j].y - stars[i].y;
            if dx * dx + dy * dy < min_sq {
                keep[j] = false;
            }
        }
    }

    let mut write = 0;
    for (read, &retained) in keep.iter().enumerate() {
        if retained {
            stars.swap(write, read);
            write += 1;
        }
    }
    stars.truncate(write);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn star(x: f64, y: f64, flux: f64) -> DetectedStar {
        DetectedStar {
            x,
            y,
            flux,
            snr: 100.0,
            saturated: false,
            saturated_count: 0,
            background: 0.0,
        }
    }

    #[test]
    fn keeps_isolated_stars() {
        let mut stars = vec![
            star(10.0, 10.0, 1000.0),
            star(50.0, 50.0, 800.0),
            star(90.0, 90.0, 600.0),
        ];
        deduplicate(&mut stars, 20.0);
        assert_eq!(stars.len(), 3);
    }

    #[test]
    fn removes_dimmer_duplicate() {
        let mut stars = vec![
            star(10.0, 10.0, 1000.0),
            star(15.0, 10.0, 500.0),
        ];
        deduplicate(&mut stars, 20.0);
        assert_eq!(stars.len(), 1);
        assert_eq!(stars[0].flux, 1000.0);
    }

    #[test]
    fn keeps_brighter_of_pair() {
        let mut stars = vec![
            star(50.0, 50.0, 2000.0),
            star(55.0, 50.0, 3000.0),
        ];
        stars.sort_by(|a, b| b.flux.partial_cmp(&a.flux).unwrap());
        deduplicate(&mut stars, 20.0);
        assert_eq!(stars.len(), 1);
        assert_eq!(stars[0].flux, 3000.0);
    }

    #[test]
    fn cluster_keeps_brightest() {
        let mut stars = vec![
            star(50.0, 50.0, 5000.0),
            star(52.0, 51.0, 4000.0),
            star(48.0, 49.0, 3000.0),
            star(51.0, 53.0, 2000.0),
        ];
        stars.sort_by(|a, b| b.flux.partial_cmp(&a.flux).unwrap());
        deduplicate(&mut stars, 20.0);
        assert_eq!(stars.len(), 1);
        assert_eq!(stars[0].flux, 5000.0);
    }

    #[test]
    fn respects_spacing_boundary() {
        let mut stars = vec![
            star(0.0, 0.0, 1000.0),
            star(20.0, 0.0, 800.0),
        ];
        deduplicate(&mut stars, 20.0);
        assert_eq!(stars.len(), 2, "stars exactly at min_spacing should both survive");
    }

    #[test]
    fn empty_input() {
        let mut stars: Vec<DetectedStar> = vec![];
        deduplicate(&mut stars, 20.0);
        assert!(stars.is_empty());
    }
}
