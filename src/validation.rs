#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValidationMode {
    None,
    Bases,
    Qualities,
    BasesAndQualities,
}
