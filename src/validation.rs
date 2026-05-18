#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ValidationMode {
    #[default]
    None,
    Bases,
    Qualities,
    BasesAndQualities,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Alphabet {
    AcgtnStrict,
    AcgtnCase,
    /// IUPAC ambiguity codes (case-insensitive) plus `.` and `-` gap characters.
    #[default]
    Iupac,
}
