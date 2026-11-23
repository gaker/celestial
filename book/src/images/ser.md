# SER Video Files

SER files hold the frames of a planetary or lunar capture, with an optional
timestamp for each frame. `SerReader` opens a file and returns frames as
slices of the file's bytes.

## Reading Frames

```rust,ignore
use celestial_images::ser::reader::SerReader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reader = SerReader::open("jupiter.ser")?;
    println!(
        "{} frames of {}x{}",
        reader.frame_count(),
        reader.header().image_width,
        reader.header().image_height
    );
    for frame in reader.iter_frames() {
        let frame = frame?;
        println!("frame {}: {} bytes, {:?}", frame.id, frame.data.len(), frame.timestamp);
    }
    Ok(())
}
```

`open` checks that the file holds every frame the header declares. A file
that is too short fails at `open`, not at the first frame read.

## Files That Change While Open

`SerReader` memory-maps the file, so frame reads do not copy the data.
Another program that changes the file while a `SerReader` has it open
affects the reader:

- If it truncates the file, the next read of a removed page makes the
  kernel stop the process with `SIGBUS` on Linux and macOS. That read can be
  a new frame read or a frame you already hold. Rust cannot turn this into
  an error.
- If it writes over the file in place, nothing reports it. Frames from
  `read_frame`, `frame_slice`, `frame_slices` and `iter_frames` borrow the
  mapped file, so frames already returned change under you.
  `read_frame_data` and `read_frame_u16` return copies, which do not change.
- If it replaces the file by renaming a new one over it, the reader is
  unaffected. It keeps reading the file it opened.

Open a SER file only after the capture program has finished writing it. To
work on a file that may be truncated or rewritten in place, copy it first and
open the copy.
