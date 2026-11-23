use std::ffi::{OsStr, OsString};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

// Tried in turn when a temporary file name is taken, for example by a file
// left behind when an earlier process was killed.
const TEMP_NAME_ATTEMPTS: u32 = 64;

static NEXT_TEMP_ID: AtomicU64 = AtomicU64::new(0);

// Writes `bytes` to a new file in the directory of `path`, flushes it to
// disk, and renames it to `path`. On any error the new file is removed and a
// file already at `path` is left unchanged. The new file takes the
// permissions of the file it replaces. A symbolic link at `path` is
// followed, so the link stays and the file it names is replaced.
pub(super) fn write_atomically(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let path = &follow_links(path)?;
    let (temp_path, file) = create_temp_file(path)?;
    match fill(file, path, bytes).and_then(|()| fs::rename(&temp_path, path)) {
        Ok(()) => Ok(()),
        Err(error) => Err(remove_after_failure(&temp_path, error)),
    }
}

// Same limit as the 40 links Linux follows while resolving a path.
const MAX_LINKS: u32 = 40;

// A dangling link resolves to the missing file it names, so the save creates
// that file.
fn follow_links(path: &Path) -> io::Result<PathBuf> {
    let mut current = path.to_path_buf();
    for _ in 0..MAX_LINKS {
        match fs::symlink_metadata(&current) {
            Ok(meta) if meta.file_type().is_symlink() => {
                let target = fs::read_link(&current)?;
                current = match current.parent() {
                    Some(parent) => parent.join(target),
                    None => target,
                };
            }
            Ok(_) => return Ok(current),
            Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(current),
            Err(error) => return Err(error),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "{} has more than {MAX_LINKS} levels of symbolic links",
            path.display()
        ),
    ))
}

fn create_temp_file(path: &Path) -> io::Result<(PathBuf, File)> {
    let name = path.file_name().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{} does not name a file", path.display()),
        )
    })?;
    let dir = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    for _ in 0..TEMP_NAME_ATTEMPTS {
        let temp_path = dir.join(temp_name(name));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temp_path)
        {
            Ok(file) => return Ok((temp_path, file)),
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        format!(
            "no unused temporary file name in {} after {} attempts",
            dir.display(),
            TEMP_NAME_ATTEMPTS
        ),
    ))
}

// A hidden name beside the target, unique within this process by a counter
// and across processes by the process id.
fn temp_name(name: &OsStr) -> OsString {
    let mut temp_name = OsString::from(".");
    temp_name.push(name);
    let id = NEXT_TEMP_ID.fetch_add(1, Ordering::Relaxed);
    temp_name.push(format!(".{}.{}.tmp", std::process::id(), id));
    temp_name
}

fn fill(mut file: File, path: &Path, bytes: &[u8]) -> io::Result<()> {
    match fs::metadata(path) {
        Ok(existing) => file.set_permissions(existing.permissions())?,
        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
        Err(error) => return Err(error),
    }
    write_contents(&mut file, bytes)?;
    file.sync_all()
}

// The error that stopped the save is returned. A failure to remove the
// temporary file is added to its message, because the file stays on disk.
fn remove_after_failure(temp_path: &Path, error: io::Error) -> io::Error {
    match fs::remove_file(temp_path) {
        Ok(()) => error,
        Err(removal) if removal.kind() == io::ErrorKind::NotFound => error,
        Err(removal) => io::Error::new(
            error.kind(),
            format!(
                "{error}; the temporary file {} was not removed: {removal}",
                temp_path.display()
            ),
        ),
    }
}

fn write_contents(file: &mut File, bytes: &[u8]) -> io::Result<()> {
    #[cfg(test)]
    if let Some(limit) = seam::failing_limit() {
        let written = limit.min(bytes.len());
        file.write_all(&bytes[..written])?;
        return Err(io::Error::other("injected write failure"));
    }
    file.write_all(bytes)
}

// Tests set a byte count on their own thread, and the next writes stop after
// that many bytes and return an error, the way a full disk does.
#[cfg(test)]
mod seam {
    use std::cell::Cell;

    thread_local! {
        static FAIL_AFTER: Cell<Option<usize>> = const { Cell::new(None) };
    }

    pub(super) fn failing_limit() -> Option<usize> {
        FAIL_AFTER.with(Cell::get)
    }

    pub(super) fn with_failing_writes<R>(limit: usize, f: impl FnOnce() -> R) -> R {
        FAIL_AFTER.with(|cell| cell.set(Some(limit)));
        let result = f();
        FAIL_AFTER.with(|cell| cell.set(None));
        result
    }
}

#[cfg(test)]
mod tests {
    use super::seam::with_failing_writes;
    use super::*;
    use crate::fits::compression::CompressionAlgorithm;
    use crate::formats::astro_image::AstroImage;
    use crate::formats::image::Image;
    use crate::formats::pixel_data::PixelData;
    use std::fs;
    use std::path::PathBuf;

    const ORIGINAL: &[u8] = b"original contents of the file";

    #[cfg(unix)]
    #[test]
    fn save_through_a_symlink_replaces_the_target_and_keeps_the_link() {
        let (dir, target) = dir_with_file("target.fits");
        let link = dir.path().join("link.fits");
        std::os::unix::fs::symlink("target.fits", &link).unwrap();
        write_atomically(&link, b"new").unwrap();
        assert!(fs::symlink_metadata(&link)
            .unwrap()
            .file_type()
            .is_symlink());
        assert_eq!(fs::read(&target).unwrap(), b"new");
    }

    #[cfg(unix)]
    #[test]
    fn save_through_a_dangling_symlink_creates_the_target() {
        let dir = tempfile::tempdir().unwrap();
        let link = dir.path().join("link.fits");
        std::os::unix::fs::symlink("missing.fits", &link).unwrap();
        write_atomically(&link, b"new").unwrap();
        assert!(fs::symlink_metadata(&link)
            .unwrap()
            .file_type()
            .is_symlink());
        assert_eq!(fs::read(dir.path().join("missing.fits")).unwrap(), b"new");
    }

    fn dir_with_file(name: &str) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join(name);
        fs::write(&path, ORIGINAL).unwrap();
        (dir, path)
    }

    fn entries(dir: &Path) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }

    fn assert_untouched(dir: &Path, path: &Path) {
        assert_eq!(fs::read(path).unwrap(), ORIGINAL);
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        assert_eq!(entries(dir), vec![name]);
    }

    #[test]
    fn writes_the_bytes_and_replaces_an_existing_file() {
        let (dir, path) = dir_with_file("out.bin");
        write_atomically(&path, b"new").unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"new");
        assert_eq!(entries(dir.path()), vec!["out.bin".to_string()]);
    }

    #[test]
    fn failed_write_keeps_the_existing_file_and_removes_the_temporary_one() {
        let (dir, path) = dir_with_file("out.bin");
        let result = with_failing_writes(3, || write_atomically(&path, &[7u8; 64]));
        assert_eq!(result.unwrap_err().to_string(), "injected write failure");
        assert_untouched(dir.path(), &path);
    }

    #[test]
    fn failed_write_to_a_new_path_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bin");
        let result = with_failing_writes(3, || write_atomically(&path, &[7u8; 64]));
        assert!(result.is_err());
        assert!(entries(dir.path()).is_empty());
    }

    #[test]
    fn failed_rename_over_a_directory_removes_the_temporary_file() {
        let dir = tempfile::tempdir().unwrap();
        let target = dir.path().join("taken");
        fs::create_dir(&target).unwrap();
        fs::write(target.join("inside"), ORIGINAL).unwrap();

        assert!(write_atomically(&target, b"new").is_err());
        assert_eq!(entries(dir.path()), vec!["taken".to_string()]);
        assert_eq!(fs::read(target.join("inside")).unwrap(), ORIGINAL);
    }

    #[test]
    fn path_without_a_file_name_is_rejected() {
        let error = write_atomically(Path::new("/"), b"new").unwrap_err();
        assert_eq!(error.kind(), io::ErrorKind::InvalidInput);
        assert_eq!(error.to_string(), "/ does not name a file");
    }

    #[cfg(unix)]
    #[test]
    fn replacement_keeps_the_permissions_of_the_existing_file() {
        use std::os::unix::fs::PermissionsExt;
        let (_dir, path) = dir_with_file("out.bin");
        fs::set_permissions(&path, fs::Permissions::from_mode(0o600)).unwrap();
        write_atomically(&path, b"new").unwrap();
        let mode = fs::metadata(&path).unwrap().permissions().mode();
        assert_eq!(mode & 0o777, 0o600);
    }

    #[test]
    fn relative_path_writes_in_the_current_directory() {
        let name = format!("atomic-write-test-{}.bin", std::process::id());
        let path = Path::new(&name);
        write_atomically(path, b"new").unwrap();
        let contents = fs::read(path);
        fs::remove_file(path).unwrap();
        assert_eq!(contents.unwrap(), b"new");
    }

    fn rgb_u8() -> Image {
        Image::new(PixelData::U8((0..48).collect()), vec![4, 4, 3])
    }

    #[test]
    fn failed_image_saves_keep_the_existing_file() {
        let mut names = vec!["out.fits", "out.xisf", "packed.fits"];
        if cfg!(feature = "standard-formats") {
            names.extend(["out.png", "out.tiff"]);
        }
        let image = rgb_u8();
        for name in names {
            let (dir, path) = dir_with_file(name);
            let result = with_failing_writes(10, || match name {
                "packed.fits" => {
                    image.save_fits_compressed(&path, CompressionAlgorithm::Rice, None)
                }
                _ => image.save(&path),
            });
            assert!(result.is_err(), "{name} saved");
            assert_untouched(dir.path(), &path);
        }
    }

    #[test]
    fn failed_astro_image_writes_keep_the_existing_file() {
        let data: Vec<f32> = (0..16).map(|i| i as f32).collect();
        let plain = AstroImage::new(&data, vec![4, 4]);
        let compressed = AstroImage::new(&data, vec![4, 4]).compressed(CompressionAlgorithm::Gzip);
        for (name, image) in [
            ("out.fits", &plain),
            ("out.xisf", &plain),
            ("packed.fits", &compressed),
        ] {
            let (dir, path) = dir_with_file(name);
            let result = with_failing_writes(10, || image.write_to(&path));
            assert!(result.is_err(), "{name} saved");
            assert_untouched(dir.path(), &path);
        }
        let (dir, path) = dir_with_file("direct.fits");
        assert!(with_failing_writes(10, || plain.write_fits(&path)).is_err());
        assert_untouched(dir.path(), &path);
        let (dir, path) = dir_with_file("direct.xisf");
        assert!(with_failing_writes(10, || plain.write_xisf(&path)).is_err());
        assert_untouched(dir.path(), &path);
    }
}
