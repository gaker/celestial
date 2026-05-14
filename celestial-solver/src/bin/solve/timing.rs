use std::fmt;
use std::time::{Duration, Instant};

pub struct Timer {
    start: Option<Instant>,
    pub elapsed: Duration,
}

impl Timer {
    fn new() -> Self {
        Self { start: None, elapsed: Duration::ZERO }
    }

    pub fn start(&mut self) {
        self.start = Some(Instant::now());
    }

    pub fn stop(&mut self) {
        if let Some(s) = self.start.take() {
            self.elapsed = s.elapsed();
        }
    }
}

impl fmt::Display for Timer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ms = self.elapsed.as_secs_f64() * 1000.0;
        if ms < 1.0 {
            write!(f, "{:.1}us", self.elapsed.as_secs_f64() * 1_000_000.0)
        } else if ms < 1000.0 {
            write!(f, "{ms:.1}ms")
        } else {
            write!(f, "{:.2}s", self.elapsed.as_secs_f64())
        }
    }
}

pub struct Timings {
    pub open: Timer,
    pub solve: Timer,
    pub save: Timer,
}

impl Timings {
    pub fn new() -> Self {
        Self {
            open: Timer::new(),
            solve: Timer::new(),
            save: Timer::new(),
        }
    }

    fn total(&self) -> Timer {
        let elapsed = self.open.elapsed
            + self.solve.elapsed
            + self.save.elapsed;
        Timer { start: None, elapsed }
    }
}

impl fmt::Display for Timings {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "--- timings ---")?;
        writeln!(f, "  open file:   {}", self.open)?;
        writeln!(f, "  solve:       {}", self.solve)?;
        writeln!(f, "  save files:  {}", self.save)?;
        write!(f, "  total:       {}", self.total())
    }
}
