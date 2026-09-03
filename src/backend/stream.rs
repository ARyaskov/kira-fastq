//! Streaming backend over a user-provided [`std::io::BufRead`].
//!
//! Adapts anything that yields bytes: stdin, a socket, a decoder from another crate. There is no
//! mapping and no random access here, so per-record lines are copied into the reader's scratch
//! buffer. The scan runs over the source's own buffer, so a line costs one `memcpy` and nothing
//! more.

use std::io::BufRead;

use crate::backend::{LineStatus, read_line_bufread};
use crate::error::FastqError;

pub(crate) struct StreamBackend {
    inner: Box<dyn BufRead + Send>,
    logical_offset: u64,
}

impl StreamBackend {
    #[inline]
    pub(crate) fn new(inner: Box<dyn BufRead + Send>) -> Self {
        Self {
            inner,
            logical_offset: 0,
        }
    }

    #[inline]
    pub(crate) fn logical_offset(&self) -> u64 {
        self.logical_offset
    }

    #[inline]
    pub(crate) fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        read_line_bufread(&mut self.inner, out, &mut self.logical_offset)
    }
}
