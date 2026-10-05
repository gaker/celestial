use super::spk::SpkSegment;
use super::SpkError;
use celestial_core::matrix::Vector3;

const MAX_LINKS: usize = 32;

pub(super) trait Sum: Copy {
    fn zero() -> Self;
    fn plus(self, other: Self) -> Self;
    fn minus(self, other: Self) -> Self;
}

impl Sum for Vector3 {
    fn zero() -> Self {
        Vector3::zeros()
    }

    fn plus(self, other: Self) -> Self {
        self + other
    }

    fn minus(self, other: Self) -> Self {
        self - other
    }
}

impl Sum for (Vector3, Vector3) {
    fn zero() -> Self {
        (Vector3::zeros(), Vector3::zeros())
    }

    fn plus(self, other: Self) -> Self {
        (self.0 + other.0, self.1 + other.1)
    }

    fn minus(self, other: Self) -> Self {
        (self.0 - other.0, self.1 - other.1)
    }
}

struct Path {
    nodes: [i32; MAX_LINKS + 1],
    links: [usize; MAX_LINKS],
    len: usize,
}

impl Path {
    fn new(body: i32) -> Self {
        Self {
            nodes: [body; MAX_LINKS + 1],
            links: [0; MAX_LINKS],
            len: 0,
        }
    }

    // Follows the last segment that covers t from each node to its center,
    // stopping at the first node in `stop` or where no segment continues.
    // Walks in place: returning a Path through Result copied it and cost
    // more than evaluating a short segment.
    fn walk(&mut self, segments: &[SpkSegment], t: f64, stop: &[i32]) -> Result<(), SpkError> {
        while !stop.contains(&self.last()) {
            let node = self.last();
            let Some(k) = segments
                .iter()
                .rposition(|s| s.body() == node && s.covers(t))
            else {
                break;
            };
            self.push(k, segments[k].center(), t)?;
        }
        Ok(())
    }

    fn push(&mut self, link: usize, center: i32, t: f64) -> Result<(), SpkError> {
        if self.nodes().contains(&center) {
            return Err(SpkError::InvalidData(format!(
                "segments covering TDB {} s form a loop through body {}",
                t, center
            )));
        }
        if self.len == MAX_LINKS {
            return Err(SpkError::InvalidData(format!(
                "segments covering TDB {} s chain body {} through more than {} centers",
                t, self.nodes[0], MAX_LINKS
            )));
        }
        self.links[self.len] = link;
        self.len += 1;
        self.nodes[self.len] = center;
        Ok(())
    }

    fn last(&self) -> i32 {
        self.nodes[self.len]
    }

    fn nodes(&self) -> &[i32] {
        &self.nodes[..=self.len]
    }
}

// SPICE's rule: the result is the target's links up to the first node it
// shares with the observer's path to the root, minus the observer's links to
// that node. Summing each side before subtracting keeps a direct pair
// bit-identical to its segment.
pub(super) fn relative<T: Sum>(
    segments: &[SpkSegment],
    (body, center): (i32, i32),
    t: f64,
    link: impl Fn(&SpkSegment) -> Result<T, SpkError>,
) -> Result<Option<T>, SpkError> {
    let mut observer = Path::new(center);
    observer.walk(segments, t, &[])?;
    let mut target = Path::new(body);
    target.walk(segments, t, observer.nodes())?;
    let Some(shared) = observer.nodes().iter().position(|&n| n == target.last()) else {
        return Ok(None);
    };
    let a = sum(segments, &target.links[..target.len], &link)?;
    let b = sum(segments, &observer.links[..shared], &link)?;
    Ok(Some(a.minus(b)))
}

fn sum<T: Sum>(
    segments: &[SpkSegment],
    links: &[usize],
    link: &impl Fn(&SpkSegment) -> Result<T, SpkError>,
) -> Result<T, SpkError> {
    let Some((&first, rest)) = links.split_first() else {
        return Ok(T::zero());
    };
    rest.iter().try_fold(link(&segments[first])?, |acc, &k| {
        Ok(acc.plus(link(&segments[k])?))
    })
}
