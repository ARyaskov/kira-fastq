#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FastqFormat {
    #[default]
    SingleLine,
    MultiLine,
}
