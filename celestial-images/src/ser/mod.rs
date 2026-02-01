mod buffer;
mod error;
mod header;
mod reader;
mod types;
mod writer;

pub use buffer::FrameBuffer;
pub use error::{Result, SerError};
pub use header::SerHeader;
pub use reader::SerReader;
pub use types::{ColorId, SerFile, SerFrame, SerFrameId, SerTimestamp};
pub use writer::SerWriter;
