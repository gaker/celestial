// A local HTTP/1.1 server for download tests. It answers one connection per scripted reply,
// in order, and closes each connection after the reply.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

#[derive(Clone)]
pub(crate) struct Reply {
    status: &'static str,
    // The Content-Length header value. It may differ from the bytes sent.
    content_length: Option<usize>,
    // Body pieces, each sent after its pause.
    chunks: Vec<(Duration, Vec<u8>)>,
}

impl Reply {
    pub(crate) fn ok(body: &[u8]) -> Self {
        Reply {
            status: "200 OK",
            content_length: Some(body.len()),
            chunks: vec![(Duration::ZERO, body.to_vec())],
        }
    }

    pub(crate) fn status(status: &'static str) -> Self {
        Reply {
            status,
            content_length: Some(0),
            chunks: Vec::new(),
        }
    }

    // Sends `body` with no Content-Length, so the client reads until the server closes.
    pub(crate) fn without_length(body: &[u8]) -> Self {
        Reply {
            status: "200 OK",
            content_length: None,
            chunks: vec![(Duration::ZERO, body.to_vec())],
        }
    }

    // Sends `body` in `pieces` equal parts with `pause` before each one.
    pub(crate) fn slow(body: &[u8], pieces: usize, pause: Duration) -> Self {
        let size = body.len().div_ceil(pieces);
        Reply {
            status: "200 OK",
            content_length: Some(body.len()),
            chunks: body.chunks(size).map(|c| (pause, c.to_vec())).collect(),
        }
    }

    // Sends the first `sent` bytes of a body whose Content-Length is `declared`, waits `stall`,
    // then closes the connection.
    pub(crate) fn stalled(declared: usize, sent: usize, stall: Duration) -> Self {
        Reply {
            status: "200 OK",
            content_length: Some(declared),
            chunks: vec![(Duration::ZERO, vec![b'x'; sent]), (stall, Vec::new())],
        }
    }
}

pub(crate) struct Server {
    pub(crate) url: String,
    replies: usize,
    requests: Arc<AtomicUsize>,
    handle: Option<JoinHandle<()>>,
}

impl Server {
    pub(crate) fn start(replies: Vec<Reply>) -> Self {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind a local test port");
        let url = format!("http://{}/", listener.local_addr().expect("local address"));
        let replies_len = replies.len();
        let requests = Arc::new(AtomicUsize::new(0));
        let counter = Arc::clone(&requests);
        let handle = thread::spawn(move || {
            for reply in replies {
                let Ok((stream, _)) = listener.accept() else {
                    return;
                };
                counter.fetch_add(1, Ordering::SeqCst);
                // The client may hang up first, which some tests cause on purpose.
                let _ = serve(stream, &reply);
            }
        });
        Server {
            url,
            replies: replies_len,
            requests,
            handle: Some(handle),
        }
    }

    pub(crate) fn requests(&self) -> usize {
        self.requests.load(Ordering::SeqCst)
    }

    // Waits up to 10 s until the server has accepted `count` connections.
    pub(crate) fn wait_for_requests(&self, count: usize) {
        let deadline = Instant::now() + Duration::from_secs(10);
        while self.requests() < count {
            assert!(
                Instant::now() < deadline,
                "only {} requests",
                self.requests()
            );
            thread::sleep(Duration::from_millis(5));
        }
    }

    // Returns the number of client requests, after sending requests of its own for the
    // replies no client asked for, so the server thread can end.
    pub(crate) fn finish(self) -> usize {
        let requests = self.requests();
        let address = self.url.trim_start_matches("http://").trim_end_matches('/');
        for _ in requests..self.replies {
            let mut stream = TcpStream::connect(address).expect("connect to the test server");
            stream
                .write_all(b"GET / HTTP/1.1\r\n\r\n")
                .expect("send a request");
            let mut ignored = Vec::new();
            stream.read_to_end(&mut ignored).expect("read the reply");
        }
        self.join();
        requests
    }

    // Waits until every scripted reply has been sent.
    pub(crate) fn join(mut self) -> usize {
        if let Some(handle) = self.handle.take() {
            handle.join().expect("test server thread panicked");
        }
        self.requests()
    }
}

fn serve(stream: TcpStream, reply: &Reply) -> std::io::Result<()> {
    let mut reader = BufReader::new(stream.try_clone()?);
    let mut line = String::new();
    while reader.read_line(&mut line)? > 2 {
        line.clear();
    }
    let mut stream = stream;
    write!(stream, "HTTP/1.1 {}\r\nConnection: close\r\n", reply.status)?;
    if let Some(len) = reply.content_length {
        write!(stream, "Content-Length: {}\r\n", len)?;
    }
    write!(stream, "\r\n")?;
    stream.flush()?;
    for (pause, chunk) in &reply.chunks {
        thread::sleep(*pause);
        stream.write_all(chunk)?;
        stream.flush()?;
    }
    Ok(())
}
